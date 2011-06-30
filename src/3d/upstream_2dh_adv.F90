#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:  upstream_2dh_adv - 2D upstream advection
! \label{sec-upstream-2dh-adv}
!
! !INTERFACE:
   subroutine upstream_2dh_adv(dt,f,uu,vv,ho,hn,hun,hvn, &
                           delxv,delyu,delxu,delyv,area_inv,az)
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
! Before or after this 2D-horizontal advection is executed,
! a 1D vertical advection
! step, possibly with another scheme, needs to be carried out.
!
! The advection terms are calculated according to (\ref{u_discr_advect}) and
! (\ref{v_discr_advect}) and the interface
! fluxes are again calculated according to (\ref{uflux_upstream}) and
! (\ref{vflux_upstream}).
!
! However, here in contrast to the one-step advection scheme with
! (\ref{adv_one_step}) implemented in {\tt upstream\_adv},
! first the complete horizontal and then the complete vertical
! advection needs to be executed. For consistency and conservation
! reasons, a partial step for the continuity equation needs
! to be executed as well. This is done as follows:
!
! \begin{equation}\label{adv_hor_step}
! \begin{array}{l}
! h^n_{i,j,k} c^n_{i,j,k} =
! h^o_{i,j,k} c^o_{i,j,k} \\ \\
! \displaystyle
! \qquad - \Delta t \Bigg(
! \frac{
! p_{i,j,k}\tilde c^u_{i,j,k}\Delta y^u_{i,j}-
! p_{i-1,j,k}\tilde c^u_{i-1,j,k}\Delta y^u_{i-1,j}
! }{\Delta x^c_{i,j}\Delta y^c_{i,j}}
! +
! \frac{
! q_{i,j,k}\tilde c^v_{i,j,k}\Delta y^v_{i,j}-
! q_{i,j-1,k}\tilde c^v_{i,j-1,k}\Delta y^v_{i,j-1}
! }{\Delta x^c_{i,j}\Delta y^c_{i,j}}\Bigg),
! \end{array}
! \end{equation}
!
! with the layer height changes
!
! \begin{equation}\label{adv_hor_step_h}
! h^n_{i,j,k} =
! h^o_{i,j,k}
! \displaystyle
! - \Delta t \Bigg(
! \frac{
! p_{i,j,k}\Delta y^u_{i,j}-
! p_{i-1,j,k}\Delta y^u_{i-1,j}
! }{\Delta x^c_{i,j}\Delta y^c_{i,j}}
! +
! \frac{
! q_{i,j,k}\Delta y^v_{i,j}-
! q_{i,j-1,k}\Delta y^v_{i,j-1}
! }{\Delta x^c_{i,j}\Delta y^c_{i,j}}\Bigg).
! \end{equation}
!
! Here, $n$ and $o$ denote values before and after this operation,
! respectively, $n$ denote intermediate values when the
! 1D advection step comes after this and $o$ denotes intermediate
! values when the 1D advection step came before this.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax
   use advection_3d, only: hi,hio
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in) :: uu(I3DFIELD),vv(I3DFIELD)
   REALTYPE, intent(in) :: ho(I3DFIELD),hn(I3DFIELD)
   REALTYPE, intent(in) :: hun(I3DFIELD),hvn(I3DFIELD)
   REALTYPE, intent(in) :: delxv(I2DFIELD),delyu(I2DFIELD)
   REALTYPE, intent(in) :: delxu(I2DFIELD),delyv(I2DFIELD)
   REALTYPE, intent(in) :: area_inv(I2DFIELD),dt
   integer, intent(in)  :: az(E2DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout)             :: f(I3DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer         :: rc,i,j,k,ii,jj
#ifdef USE_ALLOCATED_ARRAYS
   REALTYPE, dimension(:,:,:), allocatable       :: flx,fly
   REALTYPE, dimension(:,:,:), allocatable       :: cmin,cmax
#else
   REALTYPE        :: flx(I3DFIELD),fly(I3DFIELD)
   REALTYPE        :: cmin(I3DFIELD),cmax(I3DFIELD)
#endif
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'upstream_2dh_adv() # ',Ncall
#endif
#ifdef USE_ALLOCATED_ARRAYS
   allocate(flx(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'upstream_2dh_adv: Error allocating memory (flx)'

   allocate(fly(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'upstream_2dh_adv: Error allocating memory (fly)'

   allocate(cmax(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'fct_2dh: Error allocating memory (cmax)'

   allocate(cmin(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'fct_2dh: Error allocating memory (cmin)'
#endif
#ifdef SLICE_MODEL
 FATAL 'upstream_2dh_adv(): Do not use upstream_2dh_adv in SLICE_MODEL mode'
 stop
#endif

! NOTE: With the present implementation it is not necessary
!   to initialize flx and fly. Even if they
!   are set to garbage here, then it does not change the result.
!      BJB 2009-09-25.


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rc,i,j,k,ii,jj)

! OMP-NOTE: Master thread initializes, while other threads can
!    start on the following two loops.
!$OMP MASTER
   cmax = -100000*_ONE_
   cmin =  100000*_ONE_
!$OMP END MASTER

   do k=1,kmax   ! Calculating u-interface low-order fluxes !
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin-1,imax
            if (uu(i,j,k) .gt. _ZERO_) then
               flx(i,j,k)=uu(i,j,k)*f(i,j,k)
            else
               flx(i,j,k)=uu(i,j,k)*f(i+1,j,k)
            end if
         end do
      end do
!$OMP END DO NOWAIT
   end do

   do k=1,kmax   ! Calculating v-interface low-order fluxes !
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-1,jmax
         do i=imin,imax
            if (vv(i,j,k) .gt. _ZERO_) then
               fly(i,j,k)=vv(i,j,k)*f(i,j,k)
            else
               fly(i,j,k)=vv(i,j,k)*f(i,j+1,k)
            end if
         end do
      end do
!$OMP END DO NOWAIT
   end do
!$OMP BARRIER

   do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin,imax
            if (az(i,j) .eq. 1)  then
               do ii=-1,1
                  if (az(i+ii,j).eq.1) then
                     if (f(i+ii,j,k).gt.cmax(i,j,k)) cmax(i,j,k)=f(i+ii,j,k)
                     if (f(i+ii,j,k).lt.cmin(i,j,k)) cmin(i,j,k)=f(i+ii,j,k)
                  end if
               end do
               do jj=-1,1,2
                  if (az(i,j+jj).eq.1) then
                     if (f(i,j+jj,k).gt.cmax(i,j,k)) cmax(i,j,k)=f(i,j+jj,k)
                     if (f(i,j+jj,k).lt.cmin(i,j,k)) cmin(i,j,k)=f(i,j+jj,k)
                  end if
               end do
            end if
         end do
      end do
!$OMP END DO NOWAIT
   end do
!$OMP BARRIER


   do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin,imax
            if (az(i,j) .eq. 1)  then
               hio(i,j,k)=hi(i,j,k)
               hi(i,j,k)=hio(i,j,k)                               &
               -dt*((uu(i  ,j,k)*delyu(i  ,j)-uu(i-1,j,k)*delyu(i-1,j)  &
                    +vv(i,j  ,k)*delxv(i,j  )-vv(i,j-1,k)*delxv(i,j-1)  &
                    )*area_inv(i,j))
               f(i,j,k)=(f(i,j,k)*hio(i,j,k)                               &
               -dt*((flx(i  ,j,k)*delyu(i  ,j)-flx(i-1,j,k)*delyu(i-1,j)  &
                    +fly(i,j  ,k)*delxv(i,j  )-fly(i,j-1,k)*delxv(i,j-1)  &
                    )*area_inv(i,j)))/hi(i,j,k)
! Force monotonicity, this is needed here for correcting truncations errors:
               if (f(i,j,k).gt.cmax(i,j,k)) f(i,j,k)=cmax(i,j,k)
               if (f(i,j,k).lt.cmin(i,j,k)) f(i,j,k)=cmin(i,j,k)
            end if
         end do
      end do
!$OMP END DO NOWAIT
   end do

!$OMP END PARALLEL

#ifdef USE_ALLOCATED_ARRAYS
#ifdef FORTRAN90
   deallocate(flx,stat=rc)    ! work array
   if (rc /= 0) stop 'upstream_2dh_adv: Error de-allocating memory (flx)'
   deallocate(fly,stat=rc)    ! work array
   if (rc /= 0) stop 'upstream_2dh_adv: Error de-allocating memory (fly)'
#endif
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving upstream_2dh_adv()'
   write(debug,*)
#endif
   return
   end subroutine upstream_2dh_adv
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
