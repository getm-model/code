!$Id: uu_momentum_3d.F90,v 1.12 2006-12-15 09:57:50 kbk Exp $
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
! Finally, the new velocity profile is shifted such that its vertical
! integral is identical to the time integral of the vertically integrated
! transport.
! If the compiler option {\tt MUDFLAT} is defined, this fitting of profiles 
! is made with
! respect to the new surface elevation, otherwise to the 
! old surface elevation.
!
! When GETM is run as a slice model (compiler option {\tt SLICE\_MODEL}
! is activated), the result for $j=2$ is copied to $j=3$. 
!
! !USES:
   use exceptions
   use parameters, only: g,avmmol,rho_0
   use domain, only: imin,jmin,imax,jmax
   use domain, only: iimin,iimax,jjmin,jjmax,kmax,H,HU,min_depth
   use domain, only: dry_u,coru,au,av,az,ax
#if defined CURVILINEAR || defined SPHERICAL
   use domain, only: dxu,arud1,dxx,dyc,dyx,dxc
#else
   use domain, only: dx,dy
#endif
   use variables_2d, only: Uint,D
   use bdy_3d, only: do_bdy_3d
   use variables_3d, only: dt,cnpar,kumin,uu,vv,huo,hun,hvo,uuEx,ww,hvn
   use variables_3d, only: num,nuh,sseo,ssun,rru
   use variables_3d, only: ssuo
#ifndef NO_BAROCLINIC
   use variables_3d, only: idpdx
#endif
#ifdef UV_TVD
   use variables_3d, only: uadv,vadv,wadv,huadv,hvadv,hoadv,hnadv
#endif
   use halo_zones, only: update_3d_halo,wait_halo,U_TAG
   use meteo, only: tausx,airp
   use m3d, only: vel_check,min_vel,max_vel
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: n
   logical, intent(in)                 :: bdy3d
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,rc
   REALTYPE                  :: dif(1:kmax-1)
   REALTYPE                  :: auxn(1:kmax-1),auxo(1:kmax-1)
   REALTYPE                  :: a1(0:kmax),a2(0:kmax)
   REALTYPE                  :: a3(0:kmax),a4(0:kmax)
   REALTYPE                  :: Res(0:kmax),ex(0:kmax)
   REALTYPE                  :: zp,zm,zx,ResInt,Diff,Vloc
   REALTYPE                  :: gamma=g*rho_0
   REALTYPE                  :: cord_curv=_ZERO_
   integer                   :: status 
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'uu_momentum_3d() # ',Ncall
#endif

   do j=jjmin,jjmax
      do i=iimin,iimax

         if (au(i,j) .eq. 1 .or. au(i,j) .eq. 2) then
            if (kmax .gt. kumin(i,j)) then
               do k=kumin(i,j),kmax ! explicit terms
! Espelid et al. [2000], IJNME 49, 1521-1545
#ifdef NEW_CORI
                  Vloc=(vv(i  ,j  ,k)/sqrt(hvo(i  ,j  ,k))   &
                       +vv(i+1,j  ,k)/sqrt(hvo(i+1,j  ,k))   &
                       +vv(i  ,j-1,k)/sqrt(hvo(i  ,j-1,k))   &
                       +vv(i+1,j-1,k)/sqrt(hvo(i+1,j-1,k)))  &
                       *0.25*sqrt(huo(i,j,k))
#else
                  Vloc=0.25*(vv(i,j,k)+vv(i+1,j,k)+vv(i,j-1,k)+vv(i+1,j-1,k))
#endif
#if defined(SPHERICAL) || defined(CURVILINEAR)
                  cord_curv=(Vloc*(DYCIP1-DYC)-uu(i,j,k)*(DXX-DXXJM1))   &
                        /huo(i,j,k)*ARUD1
                  ex(k)=(cord_curv+coru(i,j))*Vloc
#else
                  ex(k)=coru(i,j)*Vloc
#endif
#ifdef NO_BAROCLINIC
                  ex(k)=dry_u(i,j)*(ex(k)-uuEx(i,j,k))
#else
                  ex(k)=dry_u(i,j)*(ex(k)-uuEx(i,j,k)+idpdx(i,j,k))
#endif
               end do
               ex(kmax)=ex(kmax)                                         &
                       +dry_u(i,j)*0.5*(tausx(i,j)+tausx(i+1,j))/rho_0
!     Eddy viscosity
               do k=kumin(i,j),kmax-1
                  dif(k)=0.5*(num(i,j,k)+num(i+1,j,k)) + avmmol
               end do

!     Auxilury terms, old and new time level,
               do k=kumin(i,j),kmax-1
                  auxo(k)=2*(1-cnpar)*dt*dif(k)/(huo(i,j,k+1)+huo(i,j,k))
                  auxn(k)=2*   cnpar *dt*dif(k)/(hun(i,j,k+1)+hun(i,j,k))
               end do

!     Barotropic pressure gradient
               zp=max(sseo(i+1,j),-H(i  ,j)+min(min_depth,D(i+1,j)))
               zm=max(sseo(i  ,j),-H(i+1,j)+min(min_depth,D(i  ,j)))
               zx=(zp-zm+(airp(i+1,j)-airp(i,j))/gamma)/DXU

!     Matrix elements for surface layer
               k=kmax
               a1(k)=-auxn(k-1)/hun(i,j,k-1)
               a2(k)=1.+auxn(k-1)/hun(i,j,k)
               a4(k)=uu(i,j,k  )*(1-auxo(k-1)/huo(i,j,k))              &
                    +uu(i,j,k-1)*auxo(k-1)/huo(i,j,k-1)                &
                    +dt*ex(k)                                          &
                    -dt*0.5*(huo(i,j,k)+hun(i,j,k))*g*zx

!     Matrix elements for inner layers
               do k=kumin(i,j)+1,kmax-1
                  a3(k)=-auxn(k  )/hun(i,j,k+1)
                  a1(k)=-auxn(k-1)/hun(i,j,k-1)
                  a2(k)=1.+(auxn(k)+auxn(k-1))/hun(i,j,k)
                  a4(k)=uu(i,j,k+1)*auxo(k)/huo(i,j,k+1)               &
                       +uu(i,j,k  )*(1-(auxo(k)+auxo(k-1))/huo(i,j,k)) &
                       +uu(i,j,k-1)*auxo(k-1)/huo(i,j,k-1)             &
                       +dt*ex(k)                                       &
                       -dt*0.5*(huo(i,j,k)+hun(i,j,k))*g*zx
               end do

!     Matrix elements for bottom layer
               k=kumin(i,j)
               a3(k)=-auxn(k  )/hun(i,j,k+1)
               a2(k)=1.+auxn(k)/hun(i,j,k)                             &
                     +dt*rru(i,j)/(0.5*(hun(i,j,k)+huo(i,j,k)))
               a4(k)=uu(i,j,k+1)*auxo(k)/huo(i,j,k+1)                  &
                    +uu(i,j,k  )*(1-auxo(k)/huo(i,j,k))                &
                    +dt*ex(k)                                          &
                    -dt*0.5*(huo(i,j,k)+hun(i,j,k))*g*zx

               call getm_tridiagonal(kmax,kumin(i,j),kmax,a1,a2,a3,a4,Res)

!     Transport correction: the integral of the new velocities has to
!     be the same than the transport calculated by the external mode, Uint.

               ResInt= _ZERO_
               do k=kumin(i,j),kmax
                  ResInt=ResInt+Res(k)
               end do
#ifdef MUDFLAT
               Diff=(Uint(i,j)-ResInt)/(ssun(i,j)+HU(i,j))
#else
               Diff=(Uint(i,j)-ResInt)/(ssuo(i,j)+HU(i,j))
#endif

               do k=kumin(i,j),kmax
                  uu(i,j,k)=Res(k) +hun(i,j,k)*Diff
               end do
            else  ! if (kmax .eq. kumin(i,j))
                  uu(i,j,kmax)=Uint(i,j)
            end if
         end if
      end do
end do

#ifdef SLICE_MODEL
   do i=iimin,iimax
      do k=kumin(i,2),kmax
         uu(i,3,k)=uu(i,2,k)
      end do
   end do
#endif


!  Update the halo zones
   call update_3d_halo(uu,uu,au,iimin,jjmin,iimax,jjmax,kmax,U_TAG)

   if (bdy3d) then
!      call do_bdy_3d(1,uu)
   end if

   call wait_halo(U_TAG)
   call mirror_bdy_3d(uu,U_TAG)

   if (vel_check .ne. 0 .and. mod(n,abs(vel_check)) .eq. 0) then
      call check_3d_fields(imin,jmin,imax,jmax,au,       &
                           iimin,jjmin,iimax,jjmax,kmax, &
                           kumin,uu,min_vel,max_vel,status)
      if (status .gt. 0) then
         if (vel_check .gt. 0) then
            call getm_error("uu_momentum_3d()", &
                            "out-of-bound values encountered")
         end if
         if (vel_check .lt. 0) then
            LEVEL1 'do_salinity(): ',status, &
                   ' out-of-bound values encountered'
         end if
      end if
   end if



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
