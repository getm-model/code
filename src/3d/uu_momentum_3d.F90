!$Id: uu_momentum_3d.F90,v 1.8 2006-02-10 22:41:56 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: uu_momentum_3d - $x$-momentum eq.\ \label{sec-uu-momentum-3d}
!
! !INTERFACE:
   subroutine uu_momentum_3d(bdy3d)
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
   use parameters, only: g,avmmol,rho_0
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
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical, intent(in)                 :: bdy3d
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: uu_momentum_3d.F90,v $
!  Revision 1.8  2006-02-10 22:41:56  hb
!  Source code documentation extended
!
!  Revision 1.7  2006-01-29 20:32:33  hb
!  Small LaTeX corrections to source code documentation
!
!  Revision 1.6  2006-01-28 20:07:54  hb
!  Extensions to compiler option SLICE_MODEL for better representation of zero gradients in y-direction
!
!  Revision 1.5  2004-07-28 14:58:18  hb
!  Changing subroutine calling order via MUDFLAT
!
!  Revision 1.4  2004/04/20 16:49:37  hb
!  call to coordinates moved for better consistency (see JMB)
!
!  Revision 1.3  2003/04/23 12:16:34  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 13:36:38  kbk
!  parallel support, cleaned code + NO_3D, NO_BAROCLINIC
!
!  Revision 1.1.1.1  2002/05/02 14:00:56  gotm
!  recovering after CVS crash
!
!  Revision 1.17  2001/10/26 09:11:28  bbh
!  Stresses in meteo.F90 are in N/m2 - divide by rho_0 where necessary
!
!  Revision 1.16  2001/10/12 11:43:04  bbh
!  Fixed conflicts
!
!  Revision 1.15  2001/10/12 11:39:20  bbh
!  TVD moved out of ??_momentum_3d.F90 and into uv_advect_3d.F90
!
!  Revision 1.14  2001/09/19 13:07:00  bbh
!  Moved advection related 3D fields to global allocation
!
!  Revision 1.13  2001/09/04 07:25:28  bbh
!  Coriolis based on summing up transports
!
!  Revision 1.12  2001/09/03 20:04:21  bbh
!  Allow individual advection settings for momentum, salinity and temperature
!
!  Revision 1.11  2001/09/03 12:56:58  bbh
!  Advection can now be split into different schemes for each direction
!
!  Revision 1.10  2001/08/30 08:56:26  bbh
!  Preparing for 3D boundary conditions
!
!  Revision 1.9  2001/08/27 11:51:45  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.8  2001/08/01 08:31:22  bbh
!  CURVILINEAR now implemented
!
!  Revision 1.7  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.6  2001/05/21 13:07:19  bbh
!  dt and cnpar is in variables_3d.F90
!
!  Revision 1.5  2001/05/20 07:51:40  bbh
!  Internal pressure included
!
!  Revision 1.4  2001/05/18 10:00:50  bbh
!  Added masks to calls to update_3d_halo()
!
!  Revision 1.3  2001/05/03 20:12:31  bbh
!  Use of variables_3d
!
!  Revision 1.2  2001/04/24 08:38:15  bbh
!  cosmetics - :-)
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
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
