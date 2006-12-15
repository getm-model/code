!$Id: vv_momentum_3d.F90,v 1.16 2006-12-15 10:25:42 kbk Exp $
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
! Finally, the new velocity profile is shifted such that its vertical
! integral is identical to the time integral of the vertically integrated
! transport.
! If the compiler option {\tt MUDFLAT} is defined, this fitting of profiles
! is made with
! respect to the new surface elevation, otherwise to the
! old surface elevation.
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
   use domain, only: imin,jmin,imax,jmax
   use domain, only: iimin,iimax,jjmin,jjmax,kmax,H,HV,min_depth
   use domain, only: dry_v,corv,au,av,az,ax
#if defined CURVILINEAR || defined SPHERICAL
   use domain, only: dyv,arvd1,dxc,dyx,dyc,dxx
#else
   use domain, only: dx,dy
#endif
   use variables_2d, only: Vint,D
   use bdy_3d, only: do_bdy_3d
   use variables_3d, only: dt,cnpar,kvmin,uu,vv,huo,hvo,hvn,vvEx,ww,hun
   use variables_3d, only: num,nuh,sseo,ssvn,rrv
   use variables_3d, only: ssvo
#ifdef XZ_PLUME_TEST
   use variables_3d, only: rho
#endif
#ifndef NO_BAROCLINIC
   use variables_3d, only: idpdy
#endif
#ifdef UV_TVD
   use variables_3d, only: uadv,vadv,wadv,huadv,hvadv,hoadv,hnadv
#endif
   use halo_zones, only: update_3d_halo,wait_halo,V_TAG
   use meteo, only: tausy,airp
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
   REALTYPE                  :: auxn(1:kmax-1),auxo(1:kmax-1),fuu(1:kmax)
   REALTYPE                  :: a1(0:kmax),a2(0:kmax)
   REALTYPE                  :: a3(0:kmax),a4(0:kmax)
   REALTYPE                  :: Res(0:kmax),ex(0:kmax)
   REALTYPE                  :: zp,zm,zy,ResInt,Diff,Uloc
   REALTYPE                  :: gamma=g*rho_0
   REALTYPE                  :: cord_curv=_ZERO_
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
   do j=jjmin,jjmax
      do i=iimin,iimax

         if ((av(i,j) .eq. 1) .or. (av(i,j) .eq. 2)) then

            if (kmax .gt. kvmin(i,j)) then

               do k=kvmin(i,j),kmax      ! explicit terms
! Espelid et al. [2000], IJNME 49, 1521-1545
#ifdef NEW_CORI
                Uloc=(uu(i  ,j  ,k)/sqrt(huo(i  ,j  ,k))  &
                     +uu(i-1,j  ,k)/sqrt(huo(i-1,j  ,k))  &
                     +uu(i  ,j+1,k)/sqrt(huo(i  ,j+1,k))  &
                     +uu(i-1,j+1,k)/sqrt(huo(i-1,j+1,k))) &
                     *0.25*sqrt(hvo(i,j,k))
#else
                Uloc=0.25*(uu(i,j,k)+uu(i-1,j,k)+uu(i,j+1,k)+uu(i-1,j+1,k))
#endif
#if defined(SPHERICAL) || defined(CURVILINEAR)
                  cord_curv=(vv(i,j,k)*(DYX-DYXIM1)-Uloc*(DXCJP1-DXC))     &
                        /hvo(i,j,k)*ARVD1
                  ex(k)=(cord_curv-corv(i,j))*Uloc
#else
                  ex(k)=-corv(i,j)*Uloc
#endif
#ifdef NO_BAROCLINIC
                  ex(k)=dry_v(i,j)*(ex(k)-vvEx(i,j,k))
#else
#ifdef XZ_PLUME_TEST
                  ex(k)=dry_v(i,j)*(ex(k)-vvEx(i,j,k)+idpdy(i,j,k)+yslope*hvn(i,j,k)*(rho(i,j,kmax)-rho(i,j,k)))
#else
                  ex(k)=dry_v(i,j)*(ex(k)-vvEx(i,j,k)+idpdy(i,j,k))
#endif
#endif
               end do
               ex(kmax)=ex(kmax)                                      &
                       +dry_v(i,j)*.5*(tausy(i,j)+tausy(i,j+1))/rho_0
!     Eddy viscosity
               do k=kvmin(i,j),kmax-1
                  dif(k)=0.5*(num(i,j,k)+num(i,j+1,k)) + avmmol
               end do

!     Auxiliury terms, old and new time level,
!     cnpar: Crank-Nicholson parameter
               do k=kvmin(i,j),kmax-1
                  auxo(k)=2*(1-cnpar)*dt*dif(k)/(hvo(i,j,k+1)+hvo(i,j,k))
                  auxn(k)=2*   cnpar *dt*dif(k)/(hvn(i,j,k+1)+hvn(i,j,k))
               end do

!     Barotropic pressure gradient
               zp=max(sseo(i,j+1),-H(i,j  )+min(min_depth,D(i,j+1)))
               zm=max(sseo(i,j  ),-H(i,j+1)+min(min_depth,D(i,j  )))
               zy=(zp-zm+(airp(i,j+1)-airp(i,j))/gamma)/DYV

!     Matrix elements for surface layer
               k=kmax
               a1(k)=-auxn(k-1)/hvn(i,j,k-1)
               a2(k)=1+auxn(k-1)/hvn(i,j,k)
               a4(k)=vv(i,j,k  )*(1-auxo(k-1)/hvo(i,j,k))              &
                    +vv(i,j,k-1)*auxo(k-1)/hvo(i,j,k-1)                &
                    +dt*ex(k)                                          &
                    -dt*0.5*(hvo(i,j,k)+hvn(i,j,k))*g*zy

!     Matrix elements for inner layers
               do k=kvmin(i,j)+1,kmax-1
                  a3(k)=-auxn(k  )/hvn(i,j,k+1)
                  a1(k)=-auxn(k-1)/hvn(i,j,k-1)
                  a2(k)=1+(auxn(k)+auxn(k-1))/hvn(i,j,k)
                  a4(k)=vv(i,j,k+1)*auxo(k)/hvo(i,j,k+1)               &
                       +vv(i,j,k  )*(1-(auxo(k)+auxo(k-1))/hvo(i,j,k)) &
                       +vv(i,j,k-1)*auxo(k-1)/hvo(i,j,k-1)             &
                       +dt*ex(k)                                       &
                       -dt*0.5*(hvo(i,j,k)+hvn(i,j,k))*g*zy
               end do

!     Matrix elements for bottom layer
               k=kvmin(i,j)
               a3(k)=-auxn(k  )/hvn(i,j,k+1)
               a2(k)=1+auxn(k)/hvn(i,j,k)                              &
                        +dt*rrv(i,j)/(0.5*(hvn(i,j,k)+hvo(i,j,k)))
               a4(k)=vv(i,j,k+1)*auxo(k)/hvo(i,j,k+1)                  &
                       +vv(i,j,k  )*(1-auxo(k)/hvo(i,j,k))             &
                       +dt*ex(k)                                       &
                       -dt*0.5*(hvo(i,j,k)+hvn(i,j,k))*g*zy

               call getm_tridiagonal(kmax,kvmin(i,j),kmax,a1,a2,a3,a4,Res)

!     Transport correction: the integral of the new velocities has to
!     be the same than the transport calculated by the external mode, Vint.

               ResInt= _ZERO_
               do k=kvmin(i,j),kmax
                  ResInt=ResInt+Res(k)
               end do

#ifdef MUDFLAT
               Diff=(Vint(i,j)-ResInt)/(ssvn(i,j)+HV(i,j))
#else
               Diff=(Vint(i,j)-ResInt)/(ssvo(i,j)+HV(i,j))
#endif


               do k=kvmin(i,j),kmax
                  vv(i,j,k)=Res(k)+hvn(i,j,k)*Diff
               end do
            else ! (kmax .eq. kvmin(i,j))
               vv(i,j,kmax)=Vint(i,j)
            end if
         end if
      end do
   end do

#ifdef SLICE_MODEL
   do i=iimin,iimax
      do k=kvmin(i,2),kmax
         vv(i,1,k)=vv(i,2,k)
         vv(i,3,k)=vv(i,2,k)
      end do
   end do
#endif

!  Update the halo zones
   call update_3d_halo(vv,vv,av,iimin,jjmin,iimax,jjmax,kmax,V_TAG)

   if (bdy3d) then
!      call do_bdy_3d(2,vv)
   end if

   call wait_halo(V_TAG)
   call mirror_bdy_3d(vv,V_TAG)

   if (vel_check .ne. 0 .and. mod(n,abs(vel_check)) .eq. 0) then
      call check_3d_fields(imin,jmin,imax,jmax,av,       &
                           iimin,jjmin,iimax,jjmax,kmax, &
                           kvmin,vv,min_vel,max_vel,status)
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
