!$Id: vv_momentum_3d.F90,v 1.5 2004-04-20 16:49:37 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: vv_momentum_3d() - 3D-vmomentum equation.
!
! !INTERFACE:
   subroutine vv_momentum_3d(bdy3d)
!
! !DESCRIPTION:
!
! !USES:
   use parameters, only: g,avmmol,rho_0
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
#ifndef NO_BAROCLINIC
   use variables_3d, only: idpdy
#endif
#ifdef UV_TVD
   use variables_3d, only: uadv,vadv,wadv,huadv,hvadv,hoadv,hnadv
#endif
   use halo_zones, only: update_3d_halo,wait_halo,V_TAG
   use meteo, only: tausy,airp
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
!  $Log: vv_momentum_3d.F90,v $
!  Revision 1.5  2004-04-20 16:49:37  hb
!  call to coordinates moved for better consistency (see JMB)
!
!  Revision 1.4  2003/06/29 17:06:23  kbk
!  Corv --> corv
!
!  Revision 1.3  2003/04/23 12:16:34  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 13:36:38  kbk
!  parallel support, cleaned code + NO_3D, NO_BAROCLINIC
!
!  Revision 1.1.1.1  2002/05/02 14:00:57  gotm
!  recovering after CVS crash
!
!  Revision 1.16  2001/10/26 09:11:28  bbh
!  Stresses in meteo.F90 are in N/m2 - divide by rho_0 where necessary
!
!  Revision 1.15  2001/10/12 11:43:04  bbh
!  Fixed conflicts
!
!  Revision 1.14  2001/10/12 11:39:20  bbh
!  TVD moved out of ??_momentum_3d.F90 and into uv_advect_3d.F90
!
!  Revision 1.13  2001/09/19 13:07:00  bbh
!  Moved advection related 3D fields to global allocation
!
!  Revision 1.12  2001/09/04 07:25:28  bbh
!  Coriolis based on summing up transports
!
!  Revision 1.11  2001/09/03 20:04:21  bbh
!  Allow individual advection settings for momentum, salinity and temperature
!
!  Revision 1.10  2001/09/03 12:56:58  bbh
!  Advection can now be split into different schemes for each direction
!
!  Revision 1.9  2001/08/30 08:56:26  bbh
!  Preparing for 3D boundary conditions
!
!  Revision 1.8  2001/08/27 11:51:45  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.7  2001/08/01 08:31:22  bbh
!  CURVILINEAR now implemented
!
!  Revision 1.6  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.5  2001/05/21 13:07:19  bbh
!  dt and cnpar is in variables_3d.F90
!
!  Revision 1.4  2001/05/20 07:51:40  bbh
!  Internal pressure included
!
!  Revision 1.3  2001/05/18 10:00:50  bbh
!  Added masks to calls to update_3d_halo()
!
!  Revision 1.2  2001/05/03 20:12:31  bbh
!  Use of variables_3d
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
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
                  ex(k)=dry_v(i,j)*(ex(k)-vvEx(i,j,k)+idpdy(i,j,k))
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

               Diff=(Vint(i,j)-ResInt)/(ssvo(i,j)+HV(i,j))


               do k=kvmin(i,j),kmax
                  vv(i,j,k)=Res(k)+hvn(i,j,k)*Diff
               end do
            else ! (kmax .eq. kvmin(i,j))
               vv(i,j,kmax)=Vint(i,j)
            end if
         end if
      end do
   end do

!  Update the halo zones
   call update_3d_halo(vv,vv,av,iimin,jjmin,iimax,jjmax,kmax,V_TAG)

   if (bdy3d) then
!      call do_bdy_3d(2,vv)
   end if

   call wait_halo(V_TAG)

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
