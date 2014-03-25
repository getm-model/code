#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: vv_momentum_3d - $y$-momentum eq.\ \label{sec-vv-momentum-3d}
!
! !INTERFACE:
   subroutine vv_momentum_3d(n,bdy3d)
!
! !DESCRIPTION:
!
! Here, the budget equation for layer-averaged momentum in eastern direction,
! $q_k$,
! is calculated. The physical equation is given as equation (\ref{vEq}),
! the layer-integrated equation as (\ref{vEqvi}), and after curvilinear
! transformation as (\ref{vEqviCurvi}).
! In this routine, first the Coriolis rotation term, $fp_k$ is calculated,
! either as direct transport averaging, or following \cite{ESPELIDea00}
! by using velocity averages (in case the compiler option {\tt NEW\_CORI}
! is set).
!
! As a next step, explicit forcing terms (advection, diffusion,
! internal pressure gradient, surface stresses) are added up (into the variable
! {\tt ex(k)}), the eddy viscosity is horizontally interpolated to the V-point,
! and the barotropic pressure gradient is calculated (the latter
! includes the pressure gradient correction for drying points, see
! section \ref{Section_dry}).
! Afterwards, the matrix is set up for each water column, and it is solved
! by means of a tri-diagonal matrix solver.
!
! In case that the compiler option {\tt STRUCTURE\_FRICTION} is switched on,
! the frictional effect of structures in the water column is calculated
! by adding the quadratic frictional term $C v \sqrt{u^2+v^2}$ (with a minus
! sign on
! the right hand side) numerically implicitly to the $v$-equation,
! with the friction coefficient $C$. The explicit part of this term,
! $C\sqrt{u^2+v^2}$,
! is calculated in the routine {\tt structure\_friction\_3d.F90}.
!
! Finally, the new velocity profile is shifted such that its vertical
! integral is identical to the time integral of the vertically integrated
! transport.
!
! When GETM is run as a slice model (compiler option {\tt SLICE\_MODEL}
! is activated), the result for $j=2$ is copied to $j=1$ and $j=3$.
! If the compiler option {\tt XZ\_PLUME\_TEST} is set, a slope
! of {\tt yslope} for bottom and isopycnals into the $y$-direction is
! prescribed, which has to be hard-coded as local variable.
!
! !USES:
   use exceptions
   use parameters, only: g,avmmol,rho_0
   use domain, only: rigid_lid
   use domain, only: imin,imax,jmin,jmax,kmax,H,HV,min_depth
   use domain, only: dry_v,corv,au,av,az
#if defined CURVILINEAR || defined SPHERICAL
   use domain, only: dyv,arvd1,dxc,dyx,dyc,dxx
#else
   use domain, only: dx,dy
#endif
   use bdy_3d, only: do_bdy_3d
   use variables_3d, only: Vadv,Dn
   use variables_3d, only: dt,cnpar,kvmin,uu,vv,huo,hvo,hvn,vvEx,ww,hun
   use variables_3d, only: num,nuh,sseo,Dvn,rrv
#ifdef _MOMENTUM_TERMS_
   use variables_3d, only: tdv_v,cor_v,ipg_v,epg_v,vsd_v,hsd_v
#endif
#ifdef XZ_PLUME_TEST
   use variables_3d, only: buoy
#endif
#ifdef STRUCTURE_FRICTION
   use variables_3d, only: sf
#endif
#ifndef NO_BAROCLINIC
   use variables_3d, only: idpdy
#endif
   use halo_zones, only: update_3d_halo,wait_halo,V_TAG
   use meteo, only: tausy,airp
   use m3d, only: ip_fac
   use m3d, only: vel_check,min_vel,max_vel
   use getm_timers, only: tic, toc, TIM_VVMOMENTUM, TIM_VVMOMENTUMH
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: n
   logical, intent(in)                 :: bdy3d
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   logical,save              :: first=.true.
   logical,save              :: no_shift=.false.
   integer                   :: i,j,k,rc
#ifdef NEW_CORI
   REALTYPE,dimension(I3DFIELD) :: work3d
#endif
   REALTYPE, POINTER         :: dif(:)
   REALTYPE, POINTER         :: auxn(:),auxo(:)
   REALTYPE, POINTER         :: a1(:),a2(:)
   REALTYPE, POINTER         :: a3(:),a4(:)
   REALTYPE, POINTER         :: Res(:),ex(:)
   REALTYPE                  :: zp,zm,zy,ResInt,Diff,Uloc
   REALTYPE                  :: gamma=g*rho_0
   REALTYPE                  :: cord_curv=_ZERO_
   REALTYPE                  :: gammai,rho_0i
#ifdef XZ_PLUME_TEST
   REALTYPE                  :: yslope=0.001
#endif
   integer                   :: status
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'vv_momentum_3d() # ',Ncall
#endif
   call tic(TIM_VVMOMENTUM)

   if (first) then
#ifdef NO_BAROTROPIC
      no_shift = .true.
#else
#ifdef SLICE_MODEL
      no_shift = rigid_lid
#endif
#endif
      first = .false.
   end if

   gammai=_ONE_/gamma
   rho_0i=_ONE_/rho_0


!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP    PRIVATE(i,j,k,rc,zp,zm,zy,ResInt,Diff,Uloc,cord_curv)         &
!$OMP    PRIVATE(dif,auxn,auxo,a1,a2,a3,a4,Res,ex)

 ! Each thread allocates its own HEAP storage:
   allocate(dif(1:kmax-1),stat=rc)    ! work array
   if (rc /= 0) stop 'vv_momentum_3d: Error allocating memory (dif)'

   allocate(auxn(1:kmax-1),stat=rc)    ! work array
   if (rc /= 0) stop 'vv_momentum_3d: Error allocating memory (auxn)'

   allocate(auxo(1:kmax-1),stat=rc)    ! work array
   if (rc /= 0) stop 'vv_momentum_3d: Error allocating memory (auxo)'

   allocate(a1(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'vv_momentum_3d: Error allocating memory (a1)'

   allocate(a2(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'vv_momentum_3d: Error allocating memory (a2)'

   allocate(a3(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'vv_momentum_3d: Error allocating memory (a3)'

   allocate(a4(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'vv_momentum_3d: Error allocating memory (a4)'

   allocate(Res(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'vv_momentum_3d: Error allocating memory (Res)'

   allocate(ex(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'vv_momentum_3d: Error allocating memory (ex)'

! Note: We do not need to initialize these work arrays.
!   Tested BJB 2009-09-25.

#ifdef NEW_CORI
!  Espelid et al. [2000], IJNME 49, 1521-1545
!  KK-TODO: check whether h[u|v]o should be replaced by h[u|v]n
   do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax+1
         do i=imin-1,imax
            if (au(i,j) .ne. 0) then
               work3d(i,j,k) = uu(i,j,k)/sqrt(huo(i,j,k))
            else
               work3d(i,j,k) = _ZERO_
            end if
         end do
      end do
!$OMP END DO NOWAIT
   end do
!$OMP BARRIER
#endif

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax

         if ((av(i,j) .eq. 1) .or. (av(i,j) .eq. 2)) then

            if (kmax .gt. kvmin(i,j)) then

               do k=kvmin(i,j),kmax      ! explicit terms
! Espelid et al. [2000], IJNME 49, 1521-1545
#ifdef NEW_CORI
                Uloc=(work3d(i  ,j  ,k)  &
                     +work3d(i-1,j  ,k)  &
                     +work3d(i  ,j+1,k)  &
                     +work3d(i-1,j+1,k)) &
                     *_QUART_*sqrt(hvo(i,j,k))
#else
                Uloc=_QUART_*(uu(i,j,k)+uu(i-1,j,k)+uu(i,j+1,k)+uu(i-1,j+1,k))
#endif
#if defined(SPHERICAL) || defined(CURVILINEAR)
                  cord_curv=(vv(i,j,k)*(DYX-DYXIM1)-Uloc*(DXCJP1-DXC))     &
                        /hvo(i,j,k)*ARVD1
                  ex(k)=(cord_curv-corv(i,j))*Uloc
#else
                  ex(k)=-corv(i,j)*Uloc
#endif
#ifdef _MOMENTUM_TERMS_
                  cor_v(i,j,k)=-dry_v(i,j)*ex(k)
#endif
#ifdef NO_BAROCLINIC
                  ex(k)=dry_v(i,j)*(ex(k)-vvEx(i,j,k))
#else
#ifdef XZ_PLUME_TEST
                  ex(k)=dry_v(i,j)*(ex(k)-vvEx(i,j,k)+idpdy(i,j,k)+yslope*hvn(i,j,k)*(buoy(i,j,kmax)-buoy(i,j,k)))
#else
                  ex(k)=dry_v(i,j)*(ex(k)-vvEx(i,j,k)+ip_fac*idpdy(i,j,k))
#endif
#ifdef _MOMENTUM_TERMS_
                  ipg_v(i,j,k)=-dry_v(i,j)*ip_fac*idpdy(i,j,k)
#endif
#endif
               end do
               ex(kmax)=ex(kmax)                                      &
                       +dry_v(i,j)*_HALF_*(tausy(i,j)+tausy(i,j+1))*rho_0i
!     Eddy viscosity
               do k=kvmin(i,j),kmax-1
                  dif(k)=_HALF_*(num(i,j,k)+num(i,j+1,k)) + avmmol
               end do

!     Auxiliury terms, old and new time level,
!     cnpar: Crank-Nicholson parameter
               do k=kvmin(i,j),kmax-1
                  auxo(k)=2*(1-cnpar)*dt*dif(k)/(hvo(i,j,k+1)+hvo(i,j,k))
                  auxn(k)=2*   cnpar *dt*dif(k)/(hvn(i,j,k+1)+hvn(i,j,k))
               end do

!     Barotropic pressure gradient
#ifndef NO_BAROTROPIC
               zp=max(sseo(i,j+1),-H(i,j  )+min(min_depth,Dn(i,j+1)))
               zm=max(sseo(i,j  ),-H(i,j+1)+min(min_depth,Dn(i,j  )))
               zy=(zp-zm+(airp(i,j+1)-airp(i,j))/gamma)/DYV
#else
               zy=_ZERO_
#endif

!     Matrix elements for surface layer
               k=kmax
               a1(k)=-auxn(k-1)/hvn(i,j,k-1)
               a2(k)=1+auxn(k-1)/hvn(i,j,k)
               a4(k)=vv(i,j,k  )*(1-auxo(k-1)/hvo(i,j,k))              &
                    +vv(i,j,k-1)*auxo(k-1)/hvo(i,j,k-1)                &
                    +dt*ex(k)                                          &
                    -dt*_HALF_*(hvo(i,j,k)+hvn(i,j,k))*g*zy

!     Matrix elements for inner layers
               do k=kvmin(i,j)+1,kmax-1
                  a3(k)=-auxn(k  )/hvn(i,j,k+1)
                  a1(k)=-auxn(k-1)/hvn(i,j,k-1)
                  a2(k)=_ONE_+(auxn(k)+auxn(k-1))/hvn(i,j,k)
                  a4(k)=vv(i,j,k+1)*auxo(k)/hvo(i,j,k+1)               &
                       +vv(i,j,k  )*(1-(auxo(k)+auxo(k-1))/hvo(i,j,k)) &
                       +vv(i,j,k-1)*auxo(k-1)/hvo(i,j,k-1)             &
                       +dt*ex(k)                                       &
                       -dt*_HALF_*(hvo(i,j,k)+hvn(i,j,k))*g*zy
               end do

!     Matrix elements for bottom layer
               k=kvmin(i,j)
               a3(k)=-auxn(k  )/hvn(i,j,k+1)
               a2(k) = _ONE_ + ( auxn(k) + dt*rrv(i,j) )/hvn(i,j,k)
               a4(k)=vv(i,j,k+1)*auxo(k)/hvo(i,j,k+1)                  &
                       +vv(i,j,k  )*(_ONE_-auxo(k)/hvo(i,j,k))         &
                       +dt*ex(k)                                       &
                       -dt*_HALF_*(hvo(i,j,k)+hvn(i,j,k))*g*zy

#ifdef STRUCTURE_FRICTION
               do k=kvmin(i,j),kmax
                  a2(k)=a2(k)+dt*_HALF_*(sf(i,j,k)+sf(i,j+1,k))
               end do
#endif

               call getm_tridiagonal(kmax,kvmin(i,j),kmax,a1,a2,a3,a4,Res)

!     Transport correction: the integral of the new velocities has to
!     be the same than the transport calculated by the external mode, Vadv.

               ResInt= _ZERO_
               do k=kvmin(i,j),kmax
                  ResInt=ResInt+Res(k)
               end do
               Diff=(Vadv(i,j)-ResInt)/Dvn(i,j)

#ifdef _MOMENTUM_TERMS_
               do k=kvmin(i,j),kmax
                  tdv_v(i,j,k)=vv(i,j,k)
                  epg_v(i,j,k)=_HALF_*(hvo(i,j,k)+hvn(i,j,k))*g*zy     &
                              -hvn(i,j,k)*Diff/dt
                  hsd_v(i,j,k)=hsd_v(i,j,k)*dry_v(i,j)
                  if (k .eq. kmax) then
                     vsd_v(i,j,k)=-dt*dry_v(i,j)*_HALF_*               &
                                   (tausy(i,j)+tausy(i,j+1))*rho_0i    &
                                  +auxo(k-1)*(vv(i,j,k)/hvo(i,j,k)     &
                                  -vv(i,j,k-1)/hvo(i,j,k-1))
                  end if
                  if ((k .gt. kvmin(i,j)) .and. (k .lt. kmax)) then
                     vsd_v(i,j,k)=-auxo(k)*(vv(i,j,k+1)/hvo(i,j,k+1)   &
                                  -vv(i,j,k)/hvo(i,j,k))               &
                                  +auxo(k-1)*(vv(i,j,k)/hvo(i,j,k)     &
                                  -vv(i,j,k-1)/hvo(i,j,k-1))
                  end if
                  if (k .eq. kvmin(i,j)) then
                     vsd_v(i,j,k)=-auxo(k)*(vv(i,j,k+1)/hvo(i,j,k+1)   &
                                  -vv(i,j,k)/hvo(i,j,k))
                  end if
               end do
#endif

               if (no_shift) then
                  do k=kvmin(i,j),kmax
                     vv(i,j,k)=Res(k)
                  end do
               else
                  do k=kvmin(i,j),kmax
                     vv(i,j,k)=Res(k)+hvn(i,j,k)*Diff
                  end do
               end if

#ifdef _MOMENTUM_TERMS_
               do k=kvmin(i,j),kmax
                  tdv_v(i,j,k)=(vv(i,j,k)-tdv_v(i,j,k))/dt
                  if (k .eq. kmax) then
                     vsd_v(i,j,k)=(vsd_v(i,j,k)                        &
                                  +auxn(k-1)*(vv(i,j,k)/hvn(i,j,k)     &
                                  -vv(i,j,k-1)/hvn(i,j,k-1)))/dt
                  end if
                  if ((k .gt. kvmin(i,j)) .and. (k .lt. kmax)) then
                     vsd_v(i,j,k)=(vsd_v(i,j,k)                        &
                                 -auxn(k)*(vv(i,j,k+1)/hvn(i,j,k+1)    &
                                  -vv(i,j,k)/hvn(i,j,k))               &
                                 +auxn(k-1)*(vv(i,j,k)/hvn(i,j,k)      &
                                  -vv(i,j,k-1)/hvn(i,j,k-1)))/dt
                  endif
                  if (k .eq. kvmin(i,j)) then
                     vsd_v(i,j,k)=(vsd_v(i,j,k)                        &
                                 -auxn(k)*(vv(i,j,k+1)/hvn(i,j,k+1)    &
                                  -vv(i,j,k)/hvn(i,j,k))               &
                                 +dt*rrv(i,j)*vv(i,j,k)                &
                                  /(_HALF_*(hvn(i,j,k)+hvo(i,j,k))))/dt
                  endif
               end do
#endif

            else ! (kmax .eq. kvmin(i,j))
               vv(i,j,kmax)=Vadv(i,j)
            end if
         end if
      end do
   end do
!$OMP END DO

#ifdef SLICE_MODEL
!$OMP DO SCHEDULE(RUNTIME)
   do i=imin,imax
      do k=kvmin(i,2),kmax
         vv(i,1,k)=vv(i,2,k)
         vv(i,3,k)=vv(i,2,k)
      end do
   end do
!$OMP END DO
#endif

! Each thread must deallocate its own HEAP storage:
   deallocate(dif,stat=rc)
   if (rc /= 0) stop 'vv_momentum_3d: Error deallocating memory (dif)'

   deallocate(auxn,stat=rc)
   if (rc /= 0) stop 'vv_momentum_3d: Error deallocating memory (auxn)'

   deallocate(auxo,stat=rc)
   if (rc /= 0) stop 'vv_momentum_3d: Error deallocating memory (auxo)'

   deallocate(a1,stat=rc)
   if (rc /= 0) stop 'vv_momentum_3d: Error deallocating memory (a1)'

   deallocate(a2,stat=rc)
   if (rc /= 0) stop 'vv_momentum_3d: Error deallocating memory (a2)'

   deallocate(a3,stat=rc)
   if (rc /= 0) stop 'vv_momentum_3d: Error deallocating memory (a3)'

   deallocate(a4,stat=rc)
   if (rc /= 0) stop 'vv_momentum_3d: Error deallocating memory (a4)'

   deallocate(Res,stat=rc)
   if (rc /= 0) stop 'vv_momentum_3d: Error deallocating memory (Res)'

   deallocate(ex,stat=rc)
   if (rc /= 0) stop 'vv_momentum_3d: Error deallocating memory (ex)'

!$OMP END PARALLEL

!  Update the halo zones
   call tic(TIM_VVMOMENTUMH)
   call update_3d_halo(vv,vv,av,imin,jmin,imax,jmax,kmax,V_TAG)

   if (bdy3d) then
!      call do_bdy_3d(2,vv)
   end if

   call wait_halo(V_TAG)
   call toc(TIM_VVMOMENTUMH)
   call mirror_bdy_3d(vv,V_TAG)

   if (vel_check .ne. 0 .and. mod(n,abs(vel_check)) .eq. 0) then
      call check_3d_fields(imin,jmin,imax,jmax,kvmin,kmax,av, &
                           vv,min_vel,max_vel,status)
      if (status .gt. 0) then
         if (vel_check .gt. 0) then
            call getm_error("vv_momentum_3d()", &
                            "out-of-bound values encountered")
         end if
         if (vel_check .lt. 0) then
            LEVEL1 'vv_momentum_3d(): ',status, &
                   ' out-of-bound values encountered'
         end if
      end if
   end if

   call toc(TIM_VVMOMENTUM)
#ifdef DEBUG
   write(debug,*) 'Leaving vv_momentum_3d()'
   write(debug,*)
#endif
   return
   end subroutine vv_momentum_3d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
