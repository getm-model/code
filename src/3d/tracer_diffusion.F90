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
   use domain, only: dxu,dyu,dxv,dyv,arcd1
#else
   use domain, only: dx,dy,ard1
#endif
   use variables_3d, only: dt,hn
   use variables_3d, only: diffxx,diffxy,diffyx,diffyy
   use variables_les, only: AmU_3d,AmV_3d
   use variables_2d, only: PP
   use getm_timers, only: tic,toc,TIM_TRACEDIFF

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)  :: AH_method
   REALTYPE,intent(in) :: AH_const,AH_Prt,AH_stirr_const
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
   integer                                  :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'tracer_diffusion() # ',Ncall
#endif
#ifdef SLICE_MODEL
      j=jmax/2
#endif
   call tic(TIM_TRACEDIFF)

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

!  Note (KK): diffusion only within new layer heigts

   do k=1,kmax

!     x-change at U-points
#ifndef SLICE_MODEL
      do j=jmin-HALO,jmax+HALO
#endif
         do i=imin-HALO,imax+HALO-1
            if (au(i,j) .ge. 1) then ! zero gradient across closed boundary
               dfdxU(i,j) = ( f(i+1,j,k) - f(i,j,k) ) / DXU
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif

#ifndef SLICE_MODEL
!     y-change at V-points
      do j=jmin-HALO,jmax+HALO-1
         do i=imin-HALO,imax+HALO
            if (av(i,j) .ge. 1) then ! zero normal gradient across closed boundary
               dfdyV(i,j) = ( f(i,j+1,k) - f(i,j,k) ) / DYV
            end if
         end do
      end do

      if (AH_method .eq. 3) then
!        interpolation of dfdxU to X-points
         do j=jmin-HALO,jmax+HALO-1
            do i=imin-HALO,imax+HALO-1
               if(ax(i,j) .ge. 1) then
                  PP(i,j) = _HALF_ * ( dfdxU(i,j) + dfdxU(i,j+1) )
               else
!                 Note (KK): dfdxX at corners set to 0
!                            dfdxX at closed N/S open bdys not 0, but not needed
                  PP(i,j) = _ZERO_
               end if
            end do
         end do
!        interpolation of dfdxX to V-points
         do j=jmin-HALO,jmax+HALO-1
            do i=imin-HALO+1,imax+HALO-1
!              Note (KK): we only need dfdxV(av=[1|2]), therefore settings for
!                         av=0 (due to zero gradient use of dfdxC) and for
!                         av=3 (zero gradient???) can be skipped
               if (av(i,j).eq.1 .or. av(i,j).eq.2) then
                  dfdxV(i,j) = _HALF_ * ( PP(i-1,j) + PP(i,j) )
               end if
            end do
         end do

!        interpolation of dfdyV to X-points
         do j=jmin-HALO,jmax+HALO-1
            do i=imin-HALO,imax+HALO-1
               if(ax(i,j) .ge. 1) then
                  PP(i,j) = _HALF_ * ( dfdyV(i,j) + dfdyV(i+1,j) )
               else
!                 Note (KK): dfdyX at corners set to 0
!                            dfdyX at closed W/E open bdys not 0, but not needed
                  PP(i,j) = _ZERO_
               end if
            end do
         end do
!        interpolation of dfdyX to U-points
         do j=jmin-HALO+1,jmax+HALO-1
            do i=imin-HALO,imax+HALO-1
!              Note (KK): we only need dfdyU(au=[1|2]), therefore settings for
!                         au=0 (due to zero gradient use of dfdyC) and for
!                         au=3 (zero gradient???) can be skipped
               if (au(i,j).eq.1 .or. au(i,j).eq.2) then
                  dfdyU(i,j) = _HALF_ * ( PP(i,j-1) + PP(i,j) )
               end if
            end do
         end do
      end if
#endif

#ifndef SLICE_MODEL
      do j=jmin,jmax
#endif
         do i=imin-1,imax
!           flux at U-points
!           Note (KK): flux at au=3 is not needed (AmU would be trash anyway),
!                      at closed boundaries dfdxU and diffxy are already 0
!                      therefore only flux calculation at au=[1|2]
            if(au(i,j).eq.1 .or. au(i,j).eq.2) then
               select case (AH_method)
                  case(1)
                     dfdxU(i,j) = AH_const * dfdxU(i,j)
                  case(2)
                     dfdxU(i,j) = AmU_3d(i,j,k) / AH_Prt * dfdxU(i,j)
                  case(3)
                     dfdxU(i,j) = AH_stirr_const &
                                  * (                           &
                                       diffxx(i,j,k)*dfdxU(i,j) &
#ifndef SLICE_MODEL
                                     + diffxy(i,j,k)*dfdyU(i,j) &
#endif
                                    )
               end select
               dfdxU(i,j) = DYU*_HALF_*(hn(i,j,k)+hn(i+1,j,k))*dfdxU(i,j)
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif

#ifndef SLICE_MODEL
      do j=jmin-1,jmax
         do i=imin,imax
!           flux at V-points
!           Note (KK): flux at av=3 is not needed (AmV would be trash anyway),
!                      at closed boundaries dfdyV and diffyx are already 0,
!                      therefore only flux calculation at av=[1|2]
            if(av(i,j).eq.1 .or. av(i,j).eq.2) then
               select case (AH_method)
                  case(1)
                     dfdyV(i,j) = AH_const * dfdyV(i,j)
                  case(2)
                     dfdyV(i,j) = AmV_3d(i,j,k) / AH_Prt * dfdyV(i,j)
                  case(3)
                     dfdyV(i,j) = AH_stirr_const                &
                                  * (                           &
                                       diffyy(i,j,k)*dfdyV(i,j) &
                                     + diffyx(i,j,k)*dfdxV(i,j) &
                                    )
               end select
               dfdyV(i,j) = DXV*_HALF_*(hn(i,j,k)+hn(i,j+1,k))*dfdyV(i,j)
            end if
         end do
      end do
#endif

#ifndef SLICE_MODEL
      do j=jmin,jmax
#endif
         do i=imin,imax
!           Note (KK): tracers at az=2 will be set by external boundary data
!                      therefore dfdx(au=3) and dfdy(av=3) not needed
            if (az(i,j) .eq. 1) then
               f(i,j,k) =  f(i,j,k)                            &
                         + dt * (                              &
                                   dfdxU(i,j) - dfdxU(i-1,j  ) &
#ifndef SLICE_MODEL
                                 + dfdyV(i,j) - dfdyV(i  ,j-1) &
#endif
                                )                              &
                                * ARCD1 / hn(i,j,k)
            end if
         end do
#ifndef SLICE_MODEL
      end do
#else
      f(imin:imax,j+1,k)=f(imin:imax,j,k)
#endif

   end do


   call toc(TIM_TRACEDIFF)
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
