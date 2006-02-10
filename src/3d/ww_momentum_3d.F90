!$Id: ww_momentum_3d.F90,v 1.6 2006-02-10 22:41:56 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ww_momentum_3d - continuity eq.\ \label{sec-ww-momentum-3d}
!
! !INTERFACE:
   subroutine ww_momentum_3d()
!
! !DESCRIPTION:
!
! Here, the local continuity equation is calculated in order to obtain
! the grid-related vertical velocity $\bar w_k$. An layer-integrated equation
! for this quantity is given as equation (\ref{ContiLayerInt}) which 
! has been derived
! from the differential formulation (\ref{Konti}).
!
! Since the kinematic boundary condition must hold (and is used for the
! derivation of (\ref{ContiLayerInt})), the grid-related vertical
! velocity at the surface muzst be zero, i.e.\ $\bar w_{k_{\max}}=0$.
! This is a good consistence check for the mode splitting, since this is
! only fulfilled if the vertically integrated continuity equation
! (which is the sea surface elevation equation calculated on the
! micro time step) and this local continuity equation are compatible.
!
! The physical vertical velocity is then recalculated from the grid-related
! vertical velocity by means of (\ref{conservative_w}), ... which should
! soon be coded in the routine {\tt tow} in the directory {\tt futils}.
!
! !USES:
   use domain, only: iimin,iimax,jjmin,jjmax,kmax
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: arcd1,dxv,dyu
#else
   use domain, only: dx,dy,ard1
#endif
   use variables_3d, only: dt,kmin,uu,vv,ww,ho,hn
!  #define CALC_HALO_WW
#ifndef CALC_HALO_WW
   use domain, only: az
   use halo_zones, only: update_3d_halo,wait_halo,z_TAG
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: ww_momentum_3d.F90,v $
!  Revision 1.6  2006-02-10 22:41:56  hb
!  Source code documentation extended
!
!  Revision 1.5  2003-08-14 13:04:37  kbk
!  forgot to include az
!
!  Revision 1.4  2003/08/14 12:58:00  kbk
!  update halo zones using - update_3d_halo() - or calculating
!
!  Revision 1.3  2003/04/23 12:16:34  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 13:05:11  kbk
!  cleaned code
!
!  Revision 1.1.1.1  2002/05/02 14:00:57  gotm
!  recovering after CVS crash
!
!  Revision 1.6  2001/08/27 11:50:17  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.5  2001/08/01 08:31:22  bbh
!  CURVILINEAR now implemented
!
!  Revision 1.4  2001/05/21 13:07:19  bbh
!  dt and cnpar is in variables_3d.F90
!
!  Revision 1.3  2001/05/03 20:12:31  bbh
!  Use of variables_3d
!
!  Revision 1.2  2001/05/01 07:13:27  bbh
!  use: kmax from m3d to domain
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
   REALTYPE                  :: dtm1
   integer                   :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'ww_momentum_3d() # ',Ncall
#endif

   dtm1=_ONE_/dt

   do k=1,kmax
#ifdef CALC_HALO_WW
      do j=jjmin-1,jjmax+1
         do i=iimin-1,iimax+1
#else
      do j=jjmin,jjmax
         do i=iimin,iimax
#endif
            if (k .lt. kmin(i,j)) then
               ww(i,j,kmin(i,j))= _ZERO_
            else
               ww(i,j,k)=ww(i,j,k-1)                                   &
                        -(hn(i,j,k)-ho(i  ,j  ,k))*dtm1 &
                        -((uu(i,j,k)*DYU-uu(i-1,j  ,k)*DYUIM1) &
                        +(vv(i,j,k)*DXV-vv(i  ,j-1,k)*DXVJM1))*ARCD1
           end if
         end do
      end do
   end do

#ifndef CALC_HALO_WW
   call update_3d_halo(ww,ww,az,iimin,jjmin,iimax,jjmax,kmax,z_TAG)
   call wait_halo(z_TAG)
#endif

!     Consistency test: ww(i,j,kmax) must always be zero !

#ifdef DEBUG
   write(debug,*) 'Leaving ww_momentum_3d()'
   write(debug,*)
#endif
   return
   end subroutine ww_momentum_3d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
