#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: tracer_diffusion - lateral tracer diffusion
! \label{sec-tracer-diffusion}
!
! !INTERFACE:
   subroutine tracer_diffusion(f,AH_method,AH_const,AH_Prt,AH_stirr_const)
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax,az,au,av,ax
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxc,dyc,dxu,dyu,dxv,dyv,arcd1
#else
   use domain, only: dx,dy,ard1
#endif
   use variables_3d, only: dt,ho,hn
   use variables_3d, only: diffxx,diffxy,diffyx,diffyy
#ifdef _LES_
   use variables_les, only: AmU_3d,AmV_3d
#endif

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in) :: AH_method,AH_const,AH_Prt,AH_stirr_const
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,intent(inout) :: f(I3DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
!  must be saved because of allocation outside a module
   REALTYPE,dimension(:,:),allocatable,save :: dfdxU
#ifndef SLICE_MODEL
   REALTYPE,dimension(:,:),allocatable,save :: dfdyV,dfdxV,dfdyU
#endif
   logical,save                             :: first=.true.
   integer                                  :: rc
   REALTYPE                                 :: h_loc
   integer                                  :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'tracer_diffusion() # ',Ncall
#endif

   if (first) then
      allocate(dfdxU(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'tracer_diffusion: Error allocating memory (dfdxU)'
      dfdxU=_ZERO_

#ifndef SLICE_MODEL
      allocate(dfdyV(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'tracer_diffusion: Error allocating memory (dfdyV)'
      dfdyV=_ZERO_

      allocate(dfdxV(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'tracer_diffusion: Error allocating memory (dfdxV)'
      dfdxV=_ZERO_

      allocate(dfdyU(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'tracer_diffusion: Error allocating memory (dfdyU)'
      dfdyU=_ZERO_
#endif

      first=.false.
   end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!data ranges                                    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!dfdxU (imin-HALO:imax+1   ,jmin-HALO:jmax+HALO)!!!
!!!dfdxV (imin-1   :imax+1   ,jmin-HALO:jmax+1   )!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!dfdyU (imin-HALO:imax+1   ,jmin-1   :jmax+1   )!!!
!!!dfdyV (imin-HALO:imax+HALO,jmin-HALO:jmax+1   )!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do k=1,kmax

!     x-change at U-points
#ifdef SLICE_MODEL
      j=jmax/2
#else
      do j=jmin-HALO,jmax+HALO
#endif
         do i=imin-HALO,imax+1
            if (au(i,j) .ge. 1) then ! zero normal gradient across closed boundary
               dfdxU(i,j) = ( f(i+1,j,k) - f(i,j,k) ) / DXU
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif

#ifndef SLICE_MODEL
      if (AH_method .eq. 3) then
!        interpolation of x-change to V-points
         do j=jmin-HALO,jmax+1
            do i=imax+HALO,imin-1,-1 ! do not change reverse loop order
               if (ax(i-1,j) .eq. 0) then ! zero gradient across closed boundary
                  if (au(i-1,j  ) .gt. 0) dfdxV(i-1,j) = dfdxU(i-1,j  ) !X-points
                  if (au(i-1,j+1) .gt. 0) dfdxV(i-1,j) = dfdxU(i-1,j+1) !X-points
               else
                  dfdxV(i-1,j) = _HALF_ * ( dfdxU(i-1,j) + dfdxU(i-1,j+1) ) !X-points
               end if
!              dfdxV(imax+HALO,:) will be trash
               dfdxV(i,j) = _HALF_ * ( dfdxV(i-1,j) + dfdxV(i,j) ) !V-points
            end do
         end do
      end if

!     y-change at V-points
      do j=jmin-HALO,jmax+1
         do i=imin-HALO,imax+HALO
            if (av(i,j) .ge. 1) then ! zero normal gradient across closed boundary
               dfdyV(i,j) = ( f(i,j+1,k) - f(i,j,k) ) / DYV
            end if
         end do
      end do

      if (AH_method .eq. 3) then
!        interpolation of y-change to U-points
!        this j loop cannot be threaded
         do j=jmax+HALO,jmin-1,-1 ! do not change reverse loop order
            do i=imin-HALO,imax+1
               if (ax(i,j-1) .eq. 0) then ! zero gradient across closed boundary
                  if (av(i  ,j-1) .gt. 0) dfdyU(i,j-1) = dfdyV(i  ,j-1) !X-points
                  if (av(i+1,j-1) .gt. 0) dfdyU(i,j-1) = dfdyV(i+1,j-1) !X-points
               else
                  dfdyU(i,j-1) = _HALF_ * ( dfdyV(i,j-1) + dfdyV(i+1,j-1) ) !X-points
               end if
            end do
            do i=imin-HALO,imax+1
!              dfdyU(:,jmax+HALO) will be trash
               dfdyU(i,j) = _HALF_ * ( dfdyU(i,j-1) + dfdyU(i,j) ) !U-points
            end do
         end do
      end if
#endif

#ifdef SLICE_MODEL
      j = jmax/2
#else
      do j=jmin,jmax
#endif
         do i=imin-1,imax
!           flux at U-points
!           at closed boundaries dfdxU already _ZERO_
!           but dfdyU might be nonzero
            select case (AH_method)
               case(1)
                  dfdxU(i,j) = AH_const * dfdxU(i,j)
#ifdef _LES_
               case(2)
                  dfdxU(i,j) = AmU_3d(i,j,k) / AH_Prt * dfdxU(i,j)
#endif
               case(3)
                  dfdxU(i,j) = AH_stirr_const &
                               * (                           &
                                    diffxx(i,j,k)*dfdxU(i,j) &
#ifndef SLICE_MODEL
                                  + diffxy(i,j,k)*dfdyU(i,j) &
#endif
                                 )
            end select
            if (au(i,j) .eq. 0) then ! can we rely on zero layer heights at land?
               if (az(i  ,j) .gt. 0) h_loc = _HALF_ * hn(i  ,j,k)
               if (az(i+1,j) .gt. 0) h_loc = _HALF_ * hn(i+1,j,k)
            else
               h_loc = _HALF_ * ( hn(i,j,k) + hn(i+1,j,k) )
            end if
            dfdxU(i,j) = h_loc * DYU * dfdxU(i,j)
         end do
#ifndef SLICE_MODEL
      end do
#endif

#ifndef SLICE_MODEL
      do j=jmin-1,jmax
         do i=imin,imax
!           flux at V-points
!           at closed boundaries dfdyV already _ZERO_
!           but dfdxV might be nonzero
            select case (AH_method)
               case(1)
                  dfdyV(i,j) = AH_const * dfdyV(i,j)
#ifdef _LES_
               case(2)
                  dfdyV(i,j) = AmV_3d(i,j,k) / AH_Prt * dfdyV(i,j)
#endif
               case(3)
                  dfdyV(i,j) = AH_stirr_const                &
                               * (                           &
                                    diffyy(i,j,k)*dfdyV(i,j) &
                                  + diffyx(i,j,k)*dfdxV(i,j) &
                                 )
            end select
            if (av(i,j) .eq. 0) then ! can we rely on zero layer heights at land?
               if (az(i,j  ) .gt. 0) h_loc = _HALF_ * hn(i,j  ,k)
               if (az(i,j+1) .gt. 0) h_loc = _HALF_ * hn(i,j+1,k)
            else
               h_loc = _HALF_ * ( hn(i,j,k) + hn(i,j+1,k) )
            end if
            dfdyV(i,j) = h_loc * DXV * dfdyV(i,j)
         end do
      end do
#endif

#ifdef SLICE_MODEL
      j = jmax/2
#else
      do j=jmin,jmax
#endif
         do i=imin,imax
            if (az(i,j) .eq. 1) then
               f(i,j,k)=(  ho(i,j,k)*f(i,j,k)                  &
                         + dt * (                              &
                                   dfdxU(i,j) - dfdxU(i-1,j  ) &
#ifndef SLICE_MODEL
                                 + dfdyV(i,j) - dfdyV(i  ,j-1) &
#endif
                                )                              &
                                * ARCD1                        &
                        )                                      &
                        / hn(i,j,k)
            end if
         end do
#ifndef SLICE_MODEL
      end do
#else
      f(imin:imax,jmax/2+1,k)=f(imin:imax,jmax/2,k)
#endif

   end do

#ifdef DEBUG
   write(debug,*) 'Leaving tracer_diffusion()'
   write(debug,*)
#endif
   return
   end subroutine tracer_diffusion
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
