!$Id: bottom_friction_3d.F90,v 1.4 2003-04-23 12:16:34 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: bottom_friction_3d() - 3D-bottom friction.
!
! !INTERFACE:
   subroutine bottom_friction_3d
!
! !DESCRIPTION:
!
! !USES:
   use parameters, only: kappa,avmmol
   use domain, only: iimin,iimax,jjmin,jjmax,kmax,au,av,min_depth
   use variables_2d, only: zub,zvb,zub0,zvb0
   use variables_3d, only: kumin,kvmin,uu,vv,huo,hun,hvo,hvn,rru,rrv
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
!  $Log: bottom_friction_3d.F90,v $
!  Revision 1.4  2003-04-23 12:16:34  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.3  2003/04/07 16:29:48  kbk
!  parallel support
!
!  Revision 1.1.1.1  2002/05/02 14:00:53  gotm
!  recovering after CVS crash
!
!  Revision 1.8  2001/10/22 07:52:04  bbh
!  0. -> _ZERO_
!
!  Revision 1.7  2001/07/26 13:00:41  bbh
!  *** empty log message ***
!
!  Revision 1.6  2001/07/26 12:57:47  bbh
!  Testing advection schems - using ifdef HAIDVOGEL_TEST
!
!  Revision 1.5  2001/06/25 13:15:33  bbh
!  Fixed a few typos found by the DECFOR compiler
!
!  Revision 1.4  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.3  2001/05/18 12:59:36  bbh
!  Dont update halo zones for rru and rrv
!
!  Revision 1.2  2001/05/03 20:12:31  bbh
!  Use of variables_3d
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
   integer                   :: i,j,kk
   REALTYPE                  :: r,hh,fricvel
   logical, save             :: first=.true.
   REALTYPE, save            :: uuloc(I2DFIELD)
   REALTYPE, save            :: uvloc(I2DFIELD)
   REALTYPE, save            :: vuloc(I2DFIELD)
   REALTYPE, save            :: vvloc(I2DFIELD)
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'bottom_friction_3d() # ',Ncall
#endif

   if(first) then
      uuloc = _ZERO_
      uvloc = _ZERO_
      vuloc = _ZERO_
      vvloc = _ZERO_
      first = .false.
   end if

   do j=jjmin,jjmax
      do i=iimin,iimax
         if (au(i,j) .ge. 1) then
            kk = kumin(i,j)
            uuloc(i,j)=uu(i,j,kk)/(huo(i,j,kk))
#if 0
            uvloc(i,j)=( vv(i,j  ,kk)+vv(i+1,j  ,kk)        &
                        +vv(i,j-1,kk)+vv(i+1,j-1,kk) )      &
                       /( hvo(i,j  ,kk)+hvo(i+1,j  ,kk)     &
                         +hvo(i,j-1,kk)+hvo(i+1,j-1,kk) )
#else
            uvloc(i,j)=0.25*( vv(i  ,j  ,kk)/hvo(i  ,j  ,kk)         &
                             +vv(i+1,j  ,kk)/hvo(i+1,j  ,kk)         &
                             +vv(i  ,j-1,kk)/hvo(i  ,j-1,kk)         &
                             +vv(i+1,j-1,kk)/hvo(i+1,j-1,kk) )
#endif
         else
            uuloc(i,j) = _ZERO_
            uvloc(i,j) = _ZERO_
         end if
      end do
   end do

   do j=jjmin,jjmax
      do i=iimin,iimax
         if (av(i,j) .ge. 1) then
            kk = kvmin(i,j)
#if 0
            vuloc(i,j)=( uu(i  ,j  ,kk) + uu(i-1,j  ,kk)       &
                      +  uu(i  ,j+1,kk) + uu(i-1,j+1,kk) )     &
                         /(huo(i,j  ,kk)+huo(i-1,j  ,kk)       &
                          +huo(i,j+1,kk)+huo(i-1,j+1,kk) )
#else
            vuloc(i,j)=0.25*( uu(i  ,j  ,kk)/huo(i  ,j  ,kk)    &
                            + uu(i-1,j  ,kk)/huo(i-1,j  ,kk)    &
                            + uu(i  ,j+1,kk)/huo(i  ,j+1,kk)    &
                            + uu(i-1,j+1,kk)/huo(i-1,j+1,kk) )
#endif
            vvloc(i,j)=vv(i,j,kk)/(hvo(i,j,kk))
         else
            vuloc(i,j) = _ZERO_
            vvloc(i,j) = _ZERO_
         end if
      end do
   end do

#if 1
   do j=jjmin,jjmax
      do i=iimin,iimax
         if (au(i,j) .ge. 1) then
            hh=max(min_depth/kmax,hun(i,j,kumin(i,j)))
            r=(zub(i,j)+0.5*hh)/zub(i,j)
            r=(kappa/log(r))**2
!            fricvel=sqrt(r*(uuloc(i,j)**2+uvloc(i,j)**2))
!            zub(i,j)=min(hh,zub0(i,j)+0.1*avmmol/max(avmmol,fricvel))
!            r=(zub(i,j)+0.5*hh)/zub(i,j)
!            r=(kappa/log(r))**2
            rru(i,j)=r*sqrt(uuloc(i,j)**2+uvloc(i,j)**2)
         end if
      end do
   end do
#endif

#if 1
   do j=jjmin,jjmax
      do i=iimin,iimax
         if (av(i,j) .ge. 1) then
            hh=max(min_depth/kmax,hvn(i,j,kvmin(i,j)))
            r=(zvb(i,j)+0.5*hh)/zvb(i,j)
            r=(kappa/log(r))**2
!            fricvel=sqrt(r*(vuloc(i,j)**2+vvloc(i,j)**2))
!            zvb(i,j)=min(hh,zvb0(i,j)+0.1*avmmol/max(avmmol,fricvel))
!            r=(zvb(i,j)+0.5*hh)/zvb(i,j)
!            r=(kappa/log(r))**2
            rrv(i,j)=r*sqrt(vuloc(i,j)**2+vvloc(i,j)**2)
         end if
      end do
   end do
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving bottom_friction_3d()'
   write(debug,*)
#endif
   return
   end subroutine bottom_friction_3d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
