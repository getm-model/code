!$Id: momentum.F90,v 1.5 2003-04-07 15:54:16 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: momentum() - 2D-momentum for all interior points.
!
! !INTERFACE:
   subroutine momentum(n,tausx,tausy,airp)
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)	:: n
   REALTYPE, intent(in)	:: tausx(E2DFIELD)
   REALTYPE, intent(in)	:: tausy(E2DFIELD)
   REALTYPE, intent(in)	:: airp(E2DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: momentum.F90,v $
!  Revision 1.5  2003-04-07 15:54:16  kbk
!  parallel support
!
!  Revision 1.1.1.1  2002/05/02 14:00:44  gotm
!  recovering after CVS crash
!
!  Revision 1.8  2001/10/26 09:11:28  bbh
!  Stresses in meteo.F90 are in N/m2 - divide by rho_0 where necessary
!
!  Revision 1.7  2001/10/12 09:19:14  bbh
!  SALTWEDGE_TEST - should not be here
!
!  Revision 1.6  2001/09/04 07:25:28  bbh
!  Coriolis based on summing up transports
!
!  Revision 1.5  2001/08/01 08:25:52  bbh
!  CURVILINEAR now implemented
!
!  Revision 1.4  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.3  2001/05/18 12:55:13  bbh
!  Included masks in calls to update_2d_halo()
!
!  Revision 1.2  2001/05/03 19:35:01  bbh
!  Use of variables_2d
!
!  Revision 1.1.1.1  2001/04/17 08:43:07  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
  logical	:: ufirst=.false.  
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'Momentum() # ',Ncall
#endif

   if(ufirst) then
      call umomentum(tausx,airp)
      call vmomentum(tausy,airp)
      ufirst = .false.
   else
      call vmomentum(tausy,airp)
      call umomentum(tausx,airp)
      ufirst = .true.
   end if

   return
   end subroutine momentum
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: umomentum() - 2D-momentum for all interior points.
!
! !INTERFACE:
   subroutine umomentum(tausx,airp)
!
! !DESCRIPTION:
!
! !USES:
   use parameters, only: g,rho_0
   use domain,     only: kmax,imin,imax,jmin,jmax,H,au,min_depth,Cori,dry_u
   use domain,     only: av,corv
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain,       only: dxu,arvd1,dxc,dyx
   use variables_2d, only: V
#else
   use domain,       only: dx
#endif
   use m2d, only: dtm
   use variables_2d, only: D,z,UEx,U,DU,fV,SlUx,SlRu,ru,fU,DV
   use halo_zones, only : update_2d_halo,wait_halo,U_TAG
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)	:: tausx(E2DFIELD),airp(E2DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
   integer	:: i,j
   REALTYPE	:: zx(E2DFIELD),Slr(E2DFIELD),tausu(E2DFIELD)
   REALTYPE	:: zp,zm,Uloc,Uold
   REALTYPE	:: gamma=rho_0*g
   REALTYPE	:: cord_curv=_ZERO_
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'umomentum() # ',Ncall
#endif

   do j=jmin,jmax
      do i=imin,imax
         if (au(i,j) .gt. 0) then
            zp=max(z(i+1,j),-H(i  ,j)+min(min_depth,D(i+1,j)))
            zm=max(z(i  ,j),-H(i+1,j)+min(min_depth,D(i  ,j)))
            zx(i,j)=(zp-zm+(airp(i+1,j)-airp(i,j))/gamma)/DXU
            tausu(i,j)=0.5*(tausx(i,j)+tausx(i+1,j))
         end if
      end do
   end do

   where (U .gt. 0)
      Slr=max(Slru, _ZERO_ )
   else where
      Slr=min(Slru, _ZERO_ )
   end where

   where ((au .eq. 1) .or. (au .eq. 2))
      U=(U-dtm*(g*DU*zx+dry_u*(-tausu/rho_0-fV+UEx+SlUx+Slr)))/(1+dtm*ru/DU)
   end where

!  now u is calculated
   call update_2d_halo(U,U,au,imin,jmin,imax,jmax,U_TAG)
   call wait_halo(U_TAG)
   call mirror_bdy_2d(U,U_TAG)

! Semi-implicit treatment of Coriolis force for V-momentum eq.
   do j=jmin,jmax
      do i=imin,imax
         if(av(i,j) .ge. 1) then
! Espelid et al. [2000], IJNME 49, 1521-1545
#ifdef NEW_CORI
            Uloc= &
               ( U(i,j  )/sqrt(DU(i,j  ))+ U(i-1,j  )/sqrt(DU(i-1,j  ))   &
               + U(i,j+1)/sqrt(DU(i,j+1))+ U(i-1,j+1)/sqrt(DU(i-1,j+1)))  &
               *0.25*sqrt(DV(i,j))
#else
            Uloc=0.25*( U(i,j)+ U(i-1,j)+ U(i,j+1)+ U(i-1,j+1))
#endif
#if defined(SPHERICAL) || defined(CURVILINEAR)
            cord_curv=(V(i,j)*(DYX-DYXIM1)-Uloc*(DXCJP1-DXC))/DV(i,j)*ARVD1
            fU(i,j)=(cord_curv+corv(i,j))*Uloc
#else
            fU(i,j)=corv(i,j)*Uloc
#endif
         else
            fU(i,j)= _ZERO_
         end if
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving umomentum()'
   write(debug,*)
#endif
   return
   end subroutine umomentum
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: vmomentum() - 2D-momentum for all interior points.
!
! !INTERFACE:
   subroutine vmomentum(tausy,airp)
!
! !DESCRIPTION:
!
! !USES:
   use parameters, only: g,rho_0
   use domain,     only: imin,imax,jmin,jmax,H,av,min_depth,Cori
   use domain,     only: dry_v,au,coru
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain,     only: dyv,arud1,dxx,dyc
   use m2d,        only: U
#else
   use domain,     only: dy
#endif
   use m2d,        only: dtm
   use variables_2d, only: D,z,VEx,V,DV,fU,SlVx,SlRv,rv,fV,DU
   use halo_zones, only : update_2d_halo,wait_halo,V_TAG
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)	:: tausy(E2DFIELD),airp(E2DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
   integer	:: i,j
   REALTYPE	:: zy(E2DFIELD),Slr(E2DFIELD),tausv(E2DFIELD)
   REALTYPE	:: zp,zm,Vloc
   REALTYPE	:: gamma=rho_0*g
   REALTYPE	:: cord_curv=_ZERO_
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'vmomentum() # ',Ncall
#endif

   do j=jmin,jmax
      do i=imin,imax
         if (av(i,j) .gt. 0) then
            zp=max(z(i,j+1),-H(i,j  )+min(min_depth,D(i,j+1)))
            zm=max(z(i,j  ),-H(i,j+1)+min(min_depth,D(i,j  )))
            zy(i,j)=(zp-zm+(airp(i,j+1)-airp(i,j))/gamma)/DYV
            tausv(i,j)=0.5*(tausy(i,j)+tausy(i,j+1))
         end if
      end do
   end do

   where (V .gt. 0)
      Slr=max(Slrv, _ZERO_ )
   else where
      Slr=min(Slrv, _ZERO_ )
   end where

   where ((av .eq. 1) .or. (av .eq. 2))
      V=(V-dtm*(g*DV*zy+dry_v*(-tausv/rho_0+fU+VEx+SlVx+Slr)))/(1+dtm*rv/DV)
   end where

!  now v is calculated
   call update_2d_halo(V,V,av,imin,jmin,imax,jmax,V_TAG)
   call wait_halo(V_TAG)
   call mirror_bdy_2d(V,V_TAG)

!  Semi-implicit treatment of Coriolis force for U-momentum eq.
   do j=jmin,jmax
      do i=imin,imax
         if(au(i,j) .ge. 1) then
! Espelid et al. [2000], IJNME 49, 1521-1545
#ifdef NEW_CORI
            Vloc = &
	       ( V(i,j  )/sqrt(DV(i,j  ))+ V(i+1,j  )/sqrt(DV(i+1,j  )) + &
                 V(i,j-1)/sqrt(DV(i,j-1))+ V(i+1,j-1)/sqrt(DV(i+1,j-1)))  &
                    *0.25*sqrt(DU(i,j))
#else
            Vloc = 0.25*( V(i,j)+ V(i+1,j)+ V(i,j-1)+ V(i+1,j-1))
#endif
#if defined(SPHERICAL) || defined(CURVILINEAR)
         cord_curv=(Vloc*(DYCIP1-DYC)-U(i,j)*(DXX-DXXJM1))/DU(i,j)*ARUD1
         fV(i,j)=(cord_curv+coru(i,j))*Vloc
#else
         fV(i,j)=coru(i,j)*Vloc
#endif
         else
            fV(i,j) = _ZERO_
         end if
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving vmomentum()'
   write(debug,*)
#endif
   return
   end subroutine vmomentum
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
