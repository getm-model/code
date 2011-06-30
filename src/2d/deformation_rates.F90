!$Id: deformation_rates.F90,v 1.11 2009-09-30 11:28:45 bjb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: deformation_rates - \label{deformation_rates}
!
! !INTERFACE:
   subroutine deformation_rates(U,V,DU,DV,dudxC,dvdyC,&
                                          dudyX,dvdxX,shearX,&
                                          dudxU,dvdyU,shearU,&
                                          dudxV,dvdyV,shearV )
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,az,au,av,ax
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxc,dyc,dxu,dyu,dxv,dyv,dxx,dyx
#else
   use domain, only: dx,dy
#endif


   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(in)           :: U,V,DU,DV
!
! !OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(out)          :: dudxC,dvdyC
   REALTYPE,dimension(E2DFIELD),intent(out),optional :: dudyX,dvdxX,shearX
   REALTYPE,dimension(E2DFIELD),intent(out),optional :: dudxU,dvdyU,shearU
   REALTYPE,dimension(E2DFIELD),intent(out),optional :: dudxV,dvdyV,shearV
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
!  allocated outside a module, therefore saved
   REALTYPE,dimension(:,:),allocatable,save :: u_vel,v_vel
   integer                                  :: rc
   logical,save                             :: first=.true.
   REALTYPE                                 :: dudy_loc=_ZERO_,dvdx_loc
   integer                                  :: i,j
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'deformation_rates() # ',Ncall
#endif

   if (first) then
      allocate(u_vel(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (u_vel)'
      u_vel=_ZERO_

      allocate(v_vel(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (v_vel)'
      v_vel=_ZERO_

      first = .false.
   end if

!!!!!!!!!!"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Data ranges                                    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!u_vel (imin-HALO:imax+1   ,jmin-HALO:jmax+HALO)!!!
!!!dudxC (imin-1   :imax+1   ,jmin-HALO:jmax+HALO)!!!
!!!dudxU (imin-1   :imax     ,jmin-HALO:jmax+HALO)!!!
!!!dudxV (imin-1   :imax+1   ,jmin-HALO:jmax+1   )!!!
!!! in case of metric correction dudx only jmin-1 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!v_vel (imin-HALO:imax+HALO,jmin-HALO:jmax+1   )!!!
!!!dvdyC (imin-HALO:imax+HALO,jmin-1   :jmax+1   )!!!
!!!dvdyU (imin-HALO:imax+1   ,jmin-1   :jmax+1   )!!!
!!!dvdyV (imin-HALO:imax+HALO,jmin-1   :jmax     )!!!
!!! in case of metric correction dvdy only imin-1 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!shearX(imin-HALO:imax+1   ,jmin-HALO:jmax+1   )!!!
!!!shearU(imin-HALO:imax+1   ,jmin-1   :jmax+1   )!!!
!!!shearV(imin-1   :imax+1   ,jmin-HALO:jmax+1   )!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! for exact horizontal orthogonal coordinates see Kantha and Clayson 2000, page 35

!  zonal velocity
#ifdef SLICE_MODEL
   j=jmax/2
#else
   do j=jmin-HALO,jmax+HALO
#endif
      do i=imin-HALO,imax+1
         if (au(i,j) .gt. 0) then
            u_vel(i,j) = U(i,j)/DU(i,j)
         end if
      end do
#ifndef SLICE_MODEL
   end do
#else
   u_vel(imin-HALO:imax+1,jmax/2+1) = u_vel(imin-HALO:imax+1,jmax/2)
#endif

!  meridional velocity
#ifdef SLICE_MODEL
   j=jmax/2
#else
   do j=jmin-HALO,jmax+1
#endif
      do i=imin-HALO,imax+HALO
         if (av(i,j) .gt. 0) then
            v_vel(i,j) = V(i,j)/DV(i,j)
         end if
      end do
#ifndef SLICE_MODEL
   end do
#else
   v_vel(:,jmax/2-1) = v_vel(:,jmax/2)
   v_vel(:,jmax/2+1) = v_vel(:,jmax/2)
#endif

!!!!!!!!!!!!!!!!!!!!!!!
!!!zonal strain rate!!!
!!!!!!!!!!!!!!!!!!!!!!!

!  zonal strain rate at center points
#ifdef SLICE_MODEL
   j=jmax/2
#else
   do j=jmin-1,jmax+HALO
#endif
      do i=imin-1,imax+1
         if (az(i,j).eq.1 .or. (az(i,j).eq.2 .and. (av(i,j-1).eq.2 .or. av(i,j).eq.2))) then
!           due to mirroring of velocities dudxC can be calculated in open N/S boundary cells
!           KK-TODO: dudxC=0 in open W/E boundary cells ?!
            dudxC(i,j) = (u_vel(i,j) - u_vel(i-1,j))/DXC                      &
!                       metric correction (see Griffies 2004, page 407 or Kantha and Clayson 2000, page 587)
                        +_HALF_*(v_vel(i,j)+v_vel(i,j-1))/DXC*(DXV-DXVJM1)/DYC
         end if
      end do
#ifndef SLICE_MODEL
   end do
#else
   dudxC(imin-1:imax+1,jmax/2+1) = dudxC(imin-1:imax+1,jmax/2)
#endif


   if (present(dudxU)) then
!     interpolation of zonal strain rate to U-points
#ifdef SLICE_MODEL
      j=jmax/2
#else
      do j=jmin-1,jmax+HALO
#endif
         do i=imin-1,imax
            if (au(i,j) .eq. 0) then ! no flow across closed boundary
!              KK-TODO: here we specify Neumann conditions in addition to Dirichlet conditions ?!
               if (az(i  ,j) .eq. 1) dudxU(i,j) = dudxC(i  ,j)
               if (az(i+1,j) .eq. 1) dudxU(i,j) = dudxC(i+1,j)
            else
!              KK-TODO: do we also have to set dudxU(au=2)=0 
!                       (already done for dudxC in neighbouring boundary cell)?!
               dudxU(i,j) = _HALF_*(dudxC(i,j) + dudxC(i+1,j))
            end if
         end do
#ifndef SLICE_MODEL
      end do
#else
      dudxU(imin-1:imax,jmax/2+1) = dudxU(imin-1:imax,jmax/2)
#endif
   end if

   if (present(dudxV)) then
!     interpolation of zonal strain rate to V-points
#ifdef SLICE_MODEL
      j=jmax/2
#else
      do j=jmin-1,jmax+1
#endif
         do i=imin-1,imax+1
            if (av(i,j) .ge. 1) then
!              due to mirroring this also works for av=2
               dudxV(i,j) = _HALF_*(dudxC(i,j) + dudxC(i,j+1))
            end if
         end do
#ifndef SLICE_MODEL
      end do
#else
      dudxV(imin-1:imax+1,jmax/2-1) = dudxV(imin-1:imax+1,jmax/2)
      dudxV(imin-1:imax+1,jmax/2+1) = dudxV(imin-1:imax+1,jmax/2)
#endif
   end if


#ifndef SLICE_MODEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!meridional strain rate!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  meridional strain rate at center points
   do j=jmin-1,jmax+1
      do i=imin-1,imax+HALO
         if (az(i,j).eq.1 .or. (az(i,j).eq.2 .and. (au(i-1,j).eq.2 .or. au(i,j).eq.2))) then
!           due to mirroring of velocities dvdyC can be calculated in open W/E boundary cells
!           KK-TODO: dvdyC=0 in open N/S boundary cells ?!
            dvdyC(i,j) = (v_vel(i,j) - v_vel(i,j-1))/DYC                      &
!                       metric correction (see Griffies 2004, page 407 or Kantha and Clayson 2000, page 587)
                        +_HALF_*(u_vel(i,j)+u_vel(i-1,j))/DYC*(DYU-DYUIM1)/DXC
         end if
      end do
   end do

   if (present(dvdyU)) then
!     interpolation of meriodional strain rate to U-points
      do j=jmin-1,jmax+1
         do i=imin-1,imax+1
            if (au(i,j) .ge. 1) then
!              due to mirroring this also works for au=2
               dvdyU(i,j) = _HALF_*(dvdyC(i,j) + dvdyC(i+1,j))
            end if
         end do
      end do
   end if

   if (present(dvdyV)) then
!     interpolation of meridional strain rate to V-points
      do j=jmin-1,jmax
         do i=imin-1,imax+HALO
            if (av(i,j) .eq. 0) then ! no flow across closed boundary
!              KK-TODO: here we specify Neumann conditions in addition to Dirichlet conditions ?!
               if (az(i,j  ) .eq. 1) dvdyV(i,j) = dvdyC(i,j  )
               if (az(i,j+1) .eq. 1) dvdyV(i,j) = dvdyC(i,j+1)
            else
!              KK-TODO: do we also have to set dvdyV(av=2)=0
!                       (already done for dvdyC in neighbouring boundary cell)?!
               dvdyV(i,j) = _HALF_*(dvdyC(i,j) + dvdyC(i,j+1))
            end if
         end do
      end do
   end if
#endif


!!!!!!!!!!!!!!!!
!!!shear rate!!!
!!!!!!!!!!!!!!!!

!  shear rate at X-points
#ifdef SLICE_MODEL
   j=jmax/2
#else
   do j=jmin-HALO,jmax+1
#endif
      do i=imin-HALO,imax+1
#ifndef SLICE_MODEL
         if (ax(i,j).eq.0 .and. av(i,j).eq.0 .and. av(i+1,j).eq.0) then
!           KK-TODO: do we have to set dudyX=0 for convex corner points as well?
            dudy_loc = _ZERO_ ! slip condition must be set explicitely
         else
!           original (see Kantha and Clayson 2000, page 591)
!           dudy_loc = (u_vel(i,j+1) - u_vel(i,j)) / DYX
!           metric corrected (see Griffies 2004, page 408 and Kantha and Clayson 2000, page 587)
            dudy_loc = DXX * (u_vel(i,j+1)/DXUJP1 - u_vel(i,j)/DXU) / DYX
         end if
         if (present(dudyX)) then
            dudyX(i,j) = dudy_loc
         end if
#endif
         if (ax(i,j).eq.0 .and. au(i,j).eq.0 .and. au(i,j+1).eq.0) then
!           KK-TODO: do we have to include convex corner points here as well?
            dvdx_loc = _ZERO_ ! slip condition must be set explicitely
         else
!           original
!           dvdx_loc = (v_vel(i+1,j) - v_vel(i,j))/DXX
!           metric corrected
            dvdx_loc = DYX * (v_vel(i+1,j)/DYVIP1 - v_vel(i,j)/DYV) / DXX
         end if
         if (present(dvdxX)) then
            dvdxX(i,j) = dvdx_loc
         end if

         if (present(shearX)) then
            shearX(i,j) = dudy_loc + dvdx_loc
         end if
      end do
#ifndef SLICE_MODEL
   end do
#else
   if (present(dvdxX)) then
      dvdxX(imin-HALO:imax+1,jmax/2-1) = dvdxX(imin-HALO:imax+1,jmax/2)
      dvdxX(imin-HALO:imax+1,jmax/2+1) = dvdxX(imin-HALO:imax+1,jmax/2)
   end if
   if (present(shearX)) then
      shearX(imin-HALO:imax+1,jmax/2-1) = shearX(imin-HALO:imax+1,jmax/2)
      shearX(imin-HALO:imax+1,jmax/2+1) = shearX(imin-HALO:imax+1,jmax/2)
   end if
#endif

   if (present(shearU)) then
!     interpolation of shear rate to U-points
#ifdef SLICE_MODEL
      j=jmax/2
#else
      do j=jmin-1,jmax+1
#endif
         do i=imin-HALO,imax+1
            if (au(i,j) .eq. 3)  then
!              KK-TODO: dvdyU(au=3)=0 implies dvdxU(au=3)=dvdxX ?!
!                       here also dudyU(au=3) is assumed to be 0 ?!
!                       on the other hand shearU(au=3) is not used in tracer diffusion
               if (ax(i,j-1) .gt. 0) shearU(i,j) = shearX(i,j-1)
               if (ax(i,j  ) .gt. 0) shearU(i,j) = shearX(i,j  )
            else if (au(i,j).eq.1 .or. au(i,j).eq.2) then
               shearU(i,j) = _HALF_*(shearX(i,j-1) + shearX(i,j))
            end if
         end do
#ifndef SLICE_MODEL
      end do
#else
      shearU(imin-HALO:imax+1,jmax/2+1) = shearU(imin-HALO:imax+1,jmax/2)
#endif
   end if

   if (present(shearV)) then
!     interpolation of shear rate to V-points
#ifdef SLICE_MODEL
      j=jmax/2
#else
      do j=jmin-HALO,jmax+1
#endif
         do i=imin-1,imax+1
            if (av(i,j) .eq. 3)  then
!              KK-TODO: dudxV(av=3)=0 implies dudyV(av=3)=dudyX ?!
!                       here also dvdxV(av=3) is assumed to be 0 ?!
!                       on the other hand shearV(av=3) is not used in tracer diffusion
               if (ax(i-1,j) .gt. 0) shearV(i,j) = shearX(i-1,j)
               if (ax(i  ,j) .gt. 0) shearV(i,j) = shearX(i  ,j)
            else if (av(i,j).eq.1 .or. av(i,j).eq.2) then
               shearV(i,j) = _HALF_*(shearX(i-1,j) + shearX(i,j))
            end if
         end do
#ifndef SLICE_MODEL
      end do
#else
      shearV(imin-1:imax+1,jmax/2-1) = shearV(imin-1:imax+1,jmax/2)
      shearV(imin-1:imax+1,jmax/2+1) = shearV(imin-1:imax+1,jmax/2)
#endif
   end if


#ifdef DEBUG
   write(debug,*) 'Leaving deformation_rates()'
   write(debug,*)
#endif
   return
   end subroutine deformation_rates

!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2011 - Knut Klingbeil                                  !
!-----------------------------------------------------------------------
