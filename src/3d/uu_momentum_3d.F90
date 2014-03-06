#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: uu_momentum_3d - $x$-momentum eq.\ \label{sec-uu-momentum-3d}
!
! !INTERFACE:
   subroutine uu_momentum_3d(n,bdy3d)
!
! !DESCRIPTION:
!
! Here, the budget equation for layer-averaged momentum in eastern direction,
! $p_k$,
! is calculated. The physical equation is given as equation (\ref{uEq}),
! the layer-integrated equation as (\ref{uEqvi}), and after curvilinear
! transformation as (\ref{uEqviCurvi}).
! In this routine, first the Coriolis rotation term, $fq_k$ is calculated,
! either as direct transport averaging, or following \cite{ESPELIDea00}
! by using velocity averages (in case the compiler option {\tt NEW\_CORI}
! is set).
!
! As a next step, explicit forcing terms (advection, diffusion,
! internal pressure gradient, surface stresses) are added up (into the variable
! {\tt ex(k)}), the eddy viscosity is horizontally interpolated to the U-point,
! and the barotropic pressure gradient is calculated (the latter
! includes the pressure gradient correction for drying points, see
! section \ref{Section_dry}).
! Afterwards, the matrix is set up for each water column, and it is solved
! by means of a tri-diagonal matrix solver.
!
! In case that the compiler option {\tt STRUCTURE\_FRICTION} is switched on,
! the frictional effect of structures in the water column is calculated
! by adding the quadratic frictional term $C u \sqrt{u^2+v^2}$ (with a minus sign on
! the right hand side) numerically implicitly to the $u$-equation,
! with the friction coefficient $C$. The explicit part of this term, $C\sqrt{u^2+v^2}$,
! is calculated in the routine {\tt structure\_friction\_3d.F90}.
!
! Finally, the new velocity profile is shifted such that its vertical
! integral is identical to the time integral of the vertically integrated
! transport.
!
! When GETM is run as a slice model (compiler option {\tt SLICE\_MODEL}
! is activated), the result for $j=2$ is copied to $j=3$.
!
! !USES:
   use exceptions
   use parameters, only: g,avmmol,rho_0
   use domain, only: imin,imax,jmin,jmax,kmax,H,min_depth
   use domain, only: dry_u,coru,au,av,az,ax
#if defined CURVILINEAR || defined SPHERICAL
   use domain, only: dxu,arud1,dxx,dyc,dyx,dxc
#else
   use domain, only: dx,dy
#endif
   use bdy_3d, only: do_bdy_3d
   use variables_3d, only: uuEuler,UEulerAdv,Dn
   use variables_3d, only: dt,cnpar,kumin,uu,vv,huo,hun,hvo,uuEx,ww,hvn
   use variables_3d, only: num,nuh,sseo,Dun,rru
#ifdef _MOMENTUM_TERMS_
   use variables_3d, only: tdv_u,cor_u,ipg_u,epg_u,vsd_u,hsd_u
#endif
#ifdef STRUCTURE_FRICTION
   use variables_3d, only: sf
#endif
   use variables_3d, only: idpdx
   use halo_zones, only: update_3d_halo,wait_halo,U_TAG
   use meteo, only: tausx,airp
   use waves, only: waves_method,NO_WAVES
   use variables_waves, only: uuStokes
   use m3d, only: calc_ip,ip_fac
   use m3d, only: vel_check,min_vel,max_vel
   use getm_timers, only: tic, toc, TIM_UUMOMENTUM, TIM_UUMOMENTUMH
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
   integer                   :: i,j,k,rc
   REALTYPE,dimension(I3DFIELD) :: work3d
   REALTYPE, POINTER         :: dif(:)
   REALTYPE, POINTER         :: auxn(:),auxo(:)
   REALTYPE, POINTER         :: a1(:),a2(:)
   REALTYPE, POINTER         :: a3(:),a4(:)
   REALTYPE, POINTER         :: Res(:),ex(:)
   REALTYPE                  :: zp,zm,zx,ResInt,Diff,Vloc
   REALTYPE                  :: gamma=g*rho_0
   REALTYPE                  :: cord_curv=_ZERO_
   REALTYPE                  :: gammai,rho_0i
   integer                   :: status
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'uu_momentum_3d() # ',Ncall
#endif
   call tic(TIM_UUMOMENTUM)

   gammai=_ONE_/gamma
   rho_0i=_ONE_/rho_0

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP    PRIVATE(i,j,k,rc,zp,zm,zx,ResInt,Diff,Vloc,cord_curv)         &
!$OMP    PRIVATE(dif,auxn,auxo,a1,a2,a3,a4,Res,ex)

! Each thread allocates its own HEAP storage:
   allocate(dif(1:kmax-1),stat=rc)    ! work array
   if (rc /= 0) stop 'uu_momentum_3d: Error allocating memory (dif)'

   allocate(auxn(1:kmax-1),stat=rc)    ! work array
   if (rc /= 0) stop 'uu_momentum_3d: Error allocating memory (auxn)'

   allocate(auxo(1:kmax-1),stat=rc)    ! work array
   if (rc /= 0) stop 'uu_momentum_3d: Error allocating memory (auxo)'

   allocate(a1(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'uu_momentum_3d: Error allocating memory (a1)'

   allocate(a2(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'uu_momentum_3d: Error allocating memory (a2)'

   allocate(a3(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'uu_momentum_3d: Error allocating memory (a3)'

   allocate(a4(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'uu_momentum_3d: Error allocating memory (a4)'

   allocate(Res(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'uu_momentum_3d: Error allocating memory (Res)'

   allocate(ex(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'uu_momentum_3d: Error allocating memory (ex)'

! Note: We do not need to initialize these work arrays.
!   Tested BJB 2009-09-25.

!  calculate vv at X-points (needed for Coriolis)
!  KK-TODO: check whether h[u|v]o should be replaced by h[u|v]n
   do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-1,jmax
         do i=imin,imax
            if (ax(i,j) .ne. 0) then
#ifdef NEW_CORI
! Espelid et al. [2000], IJNME 49, 1521-1545
               work3d(i,j,k) = _HALF_ * (  vv(i  ,j,k)/sqrt(hvo(i  ,j,k)) &
                                         + vv(i+1,j,k)/sqrt(hvo(i+1,j,k)) )
#else
               work3d(i,j,k) = _HALF_ * ( vv(i,j,k) + vv(i+1,j,k) )
#endif
            else
               work3d(i,j,k) = _ZERO_
            end if
         end do
      end do
!$OMP END DO NOWAIT
   end do
!$OMP BARRIER

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax

         if (au(i,j) .eq. 1 .or. au(i,j) .eq. 2) then
            if (kmax .gt. kumin(i,j)) then
!              explicit terms
               do k=kumin(i,j),kmax
!                 Coriolis
                  Vloc = _HALF_ * ( work3d(i,j-1,k) + work3d(i,j,k) )
#ifdef NEW_CORI
!                 Espelid et al. [2000], IJNME 49, 1521-1545
                  Vloc = Vloc * sqrt(huo(i,j,k))
#endif
#if defined(SPHERICAL) || defined(CURVILINEAR)
                  cord_curv=(Vloc*(DYCIP1-DYC)-uu(i,j,k)*(DXX-DXXJM1))   &
                        /huo(i,j,k)*ARUD1
                  ex(k)=(cord_curv+coru(i,j))*Vloc
#else
                  ex(k)=coru(i,j)*Vloc
#endif
#ifdef _MOMENTUM_TERMS_
                  cor_u(i,j,k)=-dry_u(i,j)*ex(k)
#endif
!                 advection / diffusion
                  ex(k) = ex(k) - uuEx(i,j,k)
               end do
!              internal pressure gradient
               if (calc_ip) then
                  do k=kumin(i,j),kmax
                     ex(k) = ex(k) + ip_fac*idpdx(i,j,k)
#ifdef _MOMENTUM_TERMS_
                     ipg_u(i,j,k)=-dry_u(i,j)*ip_fac*idpdx(i,j,k)
#endif
                  end do
               end if
!              surface stress
               ex(kmax) = ex(kmax) + _HALF_*(tausx(i,j)+tausx(i+1,j))*rho_0i
!              finalise explicit terms
               do k=kumin(i,j),kmax
                  ex(k) = dry_u(i,j) * ex(k)
               end do
!     Eddy viscosity
               do k=kumin(i,j),kmax-1
                  dif(k)=0.5*(num(i,j,k)+num(i+1,j,k)) + avmmol
               end do

!     Auxilury terms, old and new time level,
               do k=kumin(i,j),kmax-1
                  auxo(k)=_TWO_*(_ONE_-cnpar)*dt*dif(k)/(huo(i,j,k+1)+huo(i,j,k))
                  auxn(k)=_TWO_*       cnpar *dt*dif(k)/(hun(i,j,k+1)+hun(i,j,k))
               end do

!     Barotropic pressure gradient
#ifndef NO_BAROTROPIC
               zp=max(sseo(i+1,j),-H(i  ,j)+min(min_depth,Dn(i+1,j)))
               zm=max(sseo(i  ,j),-H(i+1,j)+min(min_depth,Dn(i  ,j)))
               zx=(zp-zm+(airp(i+1,j)-airp(i,j))*gammai)/DXU
#else
               zx=_ZERO_
#endif

!     Matrix elements for surface layer
               k=kmax
               a1(k)=-auxn(k-1)/hun(i,j,k-1)
               a2(k)=1.+auxn(k-1)/hun(i,j,k)
               a4(k)=uuEuler(i,j,k  )*(1-auxo(k-1)/huo(i,j,k))         &
                    +uuEuler(i,j,k-1)*auxo(k-1)/huo(i,j,k-1)           &
                    +dt*ex(k)                                          &
                    -dt*_HALF_*(huo(i,j,k)+hun(i,j,k))*g*zx

!     Matrix elements for inner layers
               do k=kumin(i,j)+1,kmax-1
                  a3(k)=-auxn(k  )/hun(i,j,k+1)
                  a1(k)=-auxn(k-1)/hun(i,j,k-1)
                  a2(k)=_ONE_+(auxn(k)+auxn(k-1))/hun(i,j,k)
                  a4(k)=uuEuler(i,j,k+1)*auxo(k)/huo(i,j,k+1)               &
                       +uuEuler(i,j,k  )*(1-(auxo(k)+auxo(k-1))/huo(i,j,k)) &
                       +uuEuler(i,j,k-1)*auxo(k-1)/huo(i,j,k-1)             &
                       +dt*ex(k)                                            &
                       -dt*_HALF_*(huo(i,j,k)+hun(i,j,k))*g*zx
               end do

!     Matrix elements for bottom layer
               k=kumin(i,j)
               a3(k)=-auxn(k  )/hun(i,j,k+1)
               a2(k) = _ONE_ + ( auxn(k) + dt*rru(i,j) )/hun(i,j,k)
               a4(k)=uuEuler(i,j,k+1)*auxo(k)/huo(i,j,k+1)             &
                    +uuEuler(i,j,k  )*(1-auxo(k)/huo(i,j,k))           &
                    +dt*ex(k)                                          &
                    -dt*_HALF_*(huo(i,j,k)+hun(i,j,k))*g*zx

#ifdef STRUCTURE_FRICTION
               do k=kumin(i,j),kmax
                  a2(k)=a2(k)+dt*_HALF_*(sf(i,j,k)+sf(i+1,j,k))
               end do
#endif

               call getm_tridiagonal(kmax,kumin(i,j),kmax,a1,a2,a3,a4,Res)

!     Transport correction: the integral of the new velocities has to
!     be the same than the transport calculated by the external mode, Uadv.

               ResInt= _ZERO_
               do k=kumin(i,j),kmax
                  ResInt=ResInt+Res(k)
               end do
               Diff=(UEulerAdv(i,j)-ResInt)/Dun(i,j)

               do k=kumin(i,j),kmax
#ifdef _MOMENTUM_TERMS_
                  tdv_u(i,j,k)=uu(i,j,k)
                  epg_u(i,j,k)=_HALF_*(huo(i,j,k)+hun(i,j,k))*g*zx     &
                              -hun(i,j,k)*Diff/dt
                  hsd_u(i,j,k)=hsd_u(i,j,k)*dry_u(i,j)
                  if (k .eq. kmax) then
                     vsd_u(i,j,k)=-dt*dry_u(i,j)*_HALF_*               &
                                   (tausx(i,j)+tausx(i+1,j))*rho_0i    &
                                  +auxo(k-1)*(uu(i,j,k)/huo(i,j,k)     &
                                  -uu(i,j,k-1)/huo(i,j,k-1))
                  end if
                  if ((k .gt. kumin(i,j)) .and. (k .lt. kmax)) then
                     vsd_u(i,j,k)=-auxo(k)*(uu(i,j,k+1)/huo(i,j,k+1)   &
                                  -uu(i,j,k)/huo(i,j,k))               &
                                  +auxo(k-1)*(uu(i,j,k)/huo(i,j,k)     &
                                  -uu(i,j,k-1)/huo(i,j,k-1))
                  end if
                  if (k .eq. kumin(i,j)) then
                     vsd_u(i,j,k)=-auxo(k)*(uu(i,j,k+1)/huo(i,j,k+1)   &
                                  -uu(i,j,k)/huo(i,j,k))
                  end if
#endif
#ifndef NO_BAROTROPIC
                  uuEuler(i,j,k)=Res(k) +hun(i,j,k)*Diff
#else
                  uuEuler(i,j,k)=Res(k)
#endif
#ifdef _MOMENTUM_TERMS_
                  tdv_u(i,j,k)=(uu(i,j,k)-tdv_u(i,j,k))/dt
                  if (k .eq. kmax) then
                     vsd_u(i,j,k)=(vsd_u(i,j,k)                        &
                                  +auxn(k-1)*(uu(i,j,k)/hun(i,j,k)     &
                                  -uu(i,j,k-1)/hun(i,j,k-1)))/dt
                  end if
                  if ((k .gt. kumin(i,j)) .and. (k .lt. kmax)) then
                     vsd_u(i,j,k)=(vsd_u(i,j,k)                        &
                                 -auxn(k)*(uu(i,j,k+1)/hun(i,j,k+1)    &
                                  -uu(i,j,k)/hun(i,j,k))               &
                                 +auxn(k-1)*(uu(i,j,k)/hun(i,j,k)      &
                                  -uu(i,j,k-1)/hun(i,j,k-1)))/dt
                  endif
                  if (k .eq. kumin(i,j)) then
                     vsd_u(i,j,k)=(vsd_u(i,j,k)                        &
                                 -auxn(k)*(uu(i,j,k+1)/hun(i,j,k+1)    &
                                  -uu(i,j,k)/hun(i,j,k))               &
                                 +dt*rru(i,j)*uu(i,j,k)                &
                                  /(_HALF_*(hun(i,j,k)+huo(i,j,k))))/dt
                  endif
#endif
               end do
            else  ! if (kmax .eq. kumin(i,j))
                  uuEuler(i,j,kmax)=UEulerAdv(i,j)
            end if
         end if
      end do
   end do
!$OMP END DO

#ifdef SLICE_MODEL
!$OMP DO SCHEDULE(RUNTIME)
   do i=imin,imax
      do k=kumin(i,2),kmax
         uuEuler(i,3,k)=uuEuler(i,2,k)
      end do
   end do
!$OMP END DO NOWAIT
#endif

! Each thread must deallocate its own HEAP storage:
   deallocate(dif,stat=rc)
   if (rc /= 0) stop 'uu_momentum_3d: Error deallocating memory (dif)'

   deallocate(auxn,stat=rc)
   if (rc /= 0) stop 'uu_momentum_3d: Error deallocating memory (auxn)'

   deallocate(auxo,stat=rc)
   if (rc /= 0) stop 'uu_momentum_3d: Error deallocating memory (auxo)'

   deallocate(a1,stat=rc)
   if (rc /= 0) stop 'uu_momentum_3d: Error deallocating memory (a1)'

   deallocate(a2,stat=rc)
   if (rc /= 0) stop 'uu_momentum_3d: Error deallocating memory (a2)'

   deallocate(a3,stat=rc)
   if (rc /= 0) stop 'uu_momentum_3d: Error deallocating memory (a3)'

   deallocate(a4,stat=rc)
   if (rc /= 0) stop 'uu_momentum_3d: Error deallocating memory (a4)'

   deallocate(Res,stat=rc)
   if (rc /= 0) stop 'uu_momentum_3d: Error deallocating memory (Res)'

   deallocate(ex,stat=rc)
   if (rc /= 0) stop 'uu_momentum_3d: Error deallocating memory (ex)'

!$OMP END PARALLEL

!  Update the halo zones
   call tic(TIM_UUMOMENTUMH)
   call update_3d_halo(uuEuler,uuEuler,au,imin,jmin,imax,jmax,kmax,U_TAG)

   if (bdy3d) then
!      call do_bdy_3d(1,uu)
!     Note (KK): modification of uu AND uuEuler necessary for waves!!!
   end if

   call wait_halo(U_TAG)
   call toc(TIM_UUMOMENTUMH)
   call mirror_bdy_3d(uuEuler,U_TAG)

   if (waves_method .ne. NO_WAVES) then
      uu = uuEuler + uuStokes
   end if

   if (vel_check .ne. 0 .and. mod(n,abs(vel_check)) .eq. 0) then
      call check_3d_fields(imin,jmin,imax,jmax,kumin,kmax,au, &
                           uu,min_vel,max_vel,status)
      if (status .gt. 0) then
         if (vel_check .gt. 0) then
            call getm_error("uu_momentum_3d()", &
                            "out-of-bound values encountered")
         end if
         if (vel_check .lt. 0) then
            LEVEL1 'uu_momentum_3d(): ',status, &
                   ' out-of-bound values encountered'
         end if
      end if
   end if

   call toc(TIM_UUMOMENTUM)
#ifdef DEBUG
   write(debug,*) 'Leaving uu_momentum_3d()'
   write(debug,*)
#endif
   return
   end subroutine uu_momentum_3d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
