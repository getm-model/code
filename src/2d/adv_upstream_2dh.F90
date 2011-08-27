#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:  adv_upstream_2dh - 2D upstream advection 2d
! \label{sec-upstream-2dh-adv}
!
! !INTERFACE:
   subroutine adv_upstream_2dh(dt,f,Di,U,V,Do,Dn,DU,DV, &
                               delxv,delyu,delxu,delyv,area_inv,az,AH)
!  Note (KK): Keep in sync with interface in advection.F90
!
! !DESCRIPTION:
! In this routine, the first-order upstream advection scheme
! is applied for the two horizontal directions in
! one step. The scheme should be positive definite and of
! high resolution. In order to remove truncation errors which might in Wadden
! Sea applications cause non-monotonicity, a truncation of over- and
! undershoots is carried out at the end of this subroutine. Such
! two-dimensional schemes are advantageous in Wadden Sea applications, since
! one-dimensional directional-split schemes might compute negative intermediate
! depths.
!
! The advection terms are calculated according to (\ref{u_discr_advect}) and
! (\ref{v_discr_advect}) and the interface
! fluxes are again calculated according to (\ref{uflux_upstream}) and
! (\ref{vflux_upstream}).
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                        :: dt,AH
   REALTYPE,dimension(E2DFIELD),intent(in)    :: U,V,Do,Dn,DU,DV
   REALTYPE,dimension(E2DFIELD),intent(in)    :: delxv,delyu,delxu,delyv
   REALTYPE,dimension(E2DFIELD),intent(in)    :: area_inv
   integer,dimension(E2DFIELD),intent(in)     :: az
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(inout) :: f,Di
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer         :: rc,i,j,ii,jj
   REALTYPE        :: Dio
#ifdef USE_ALLOCATED_ARRAYS
   REALTYPE, dimension(:,:), allocatable       :: flx,fly
   REALTYPE, dimension(:,:), allocatable       :: cmin,cmax
#else
   REALTYPE        :: flx(I2DFIELD),fly(I2DFIELD)
   REALTYPE        :: cmin(I2DFIELD),cmax(I2DFIELD)
#endif
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'adv_upstream_2dh() # ',Ncall
#endif
#ifdef USE_ALLOCATED_ARRAYS
   allocate(flx(I2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'adv_upstream_2dh: Error allocating memory (flx)'

   allocate(fly(I2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'adv_upstream_2dh: Error allocating memory (fly)'

   allocate(cmax(I2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'adv_upstream_2dh: Error allocating memory (cmax)'

   allocate(cmin(I2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'adv_upstream_2dh: Error allocating memory (cmin)'
#endif
#ifdef SLICE_MODEL
 FATAL 'adv_upstream_2dh(): Do not use adv_upstream_2dh in SLICE_MODEL mode'
 stop
#endif

! NOTE: With the present implementation it is not necessary
!   to initialize flx and fly. Even if they
!   are set to garbage here, then it does not change the result.
!      BJB 2009-09-25.


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rc,i,j,ii,jj,Dio)

! OMP-NOTE: Master thread initializes, while other threads can
!    start on the following two loops.
!$OMP MASTER
   cmax = -100000*_ONE_
   cmin =  100000*_ONE_
!$OMP END MASTER

!  Calculating u-interface low-order fluxes !
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin-1,imax
         if (U(i,j) .gt. _ZERO_) then
            flx(i,j)=U(i,j)*f(i,j)
         else
            flx(i,j)=U(i,j)*f(i+1,j)
         end if
         if ((AH.gt._ZERO_).and.(az(i,j).gt.0).and.(az(i+1,j).gt.0))    &
            flx(i,j)=flx(i,j)-AH*(f(i+1,j)-f(i,j))/delxu(i,j) &
                      *_HALF_*(Dn(i+1,j)+Dn(i,j))
      end do
   end do
!$OMP END DO NOWAIT

!  Calculating v-interface low-order fluxes !
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-1,jmax
      do i=imin,imax
         if (V(i,j) .gt. _ZERO_) then
            fly(i,j)=V(i,j)*f(i,j)
         else
            fly(i,j)=V(i,j)*f(i,j+1)
         end if
         if ((AH.gt._ZERO_).and.(az(i,j).gt.0).and.(az(i,j+1).gt.0))   &
            fly(i,j)=fly(i,j)-AH*(f(i,j+1)-f(i,j))/delyv(i,j)   &
                      *_HALF_*(Dn(i,j+1)+Dn(i,j))
      end do
   end do
!$OMP END DO NOWAIT

!$OMP BARRIER

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1)  then
            do ii=-1,1
               if (az(i+ii,j).eq.1) then
                  if (f(i+ii,j).gt.cmax(i,j)) cmax(i,j)=f(i+ii,j)
                  if (f(i+ii,j).lt.cmin(i,j)) cmin(i,j)=f(i+ii,j)
               end if
            end do
            do jj=-1,1,2
               if (az(i,j+jj).eq.1) then
                  if (f(i,j+jj).gt.cmax(i,j)) cmax(i,j)=f(i,j+jj)
                  if (f(i,j+jj).lt.cmin(i,j)) cmin(i,j)=f(i,j+jj)
               end if
            end do
         end if
      end do
   end do
!$OMP END DO NOWAIT

!$OMP BARRIER


!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1)  then
            Dio=Di(i,j)
            Di(i,j)=Dio                               &
            -dt*((U(i  ,j)*delyu(i  ,j)-U(i-1,j)*delyu(i-1,j)  &
                 +V(i,j  )*delxv(i,j  )-V(i,j-1)*delxv(i,j-1)  &
                 )*area_inv(i,j))
            f(i,j)=(f(i,j)*Dio                               &
            -dt*((flx(i  ,j)*delyu(i  ,j)-flx(i-1,j)*delyu(i-1,j)  &
                 +fly(i,j  )*delxv(i,j  )-fly(i,j-1)*delxv(i,j-1)  &
                 )*area_inv(i,j)))/Di(i,j)
! Force monotonicity, this is needed here for correcting truncations errors:
            if (f(i,j).gt.cmax(i,j)) f(i,j)=cmax(i,j)
            if (f(i,j).lt.cmin(i,j)) f(i,j)=cmin(i,j)
         end if
      end do
   end do
!$OMP END DO NOWAIT


!$OMP END PARALLEL

#ifdef USE_ALLOCATED_ARRAYS
#ifdef FORTRAN90
   deallocate(flx,stat=rc)    ! work array
   if (rc /= 0) stop 'adv_upstream_2dh: Error de-allocating memory (flx)'
   deallocate(fly,stat=rc)    ! work array
   if (rc /= 0) stop 'adv_upstream_2dh: Error de-allocating memory (fly)'
#endif
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving adv_upstream_2dh()'
   write(debug,*)
#endif
   return
   end subroutine adv_upstream_2dh
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
