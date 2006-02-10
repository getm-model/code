!$Id: slow_bottom_friction.F90,v 1.6 2006-02-10 22:41:56 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: slow_bottom_friction() - slow bottom friction
! \label{sec-slow-bottom-friction}
!
! !INTERFACE:
   subroutine slow_bottom_friction
!
! !DESCRIPTION:
!
! !USES:
   use parameters, only: kappa
   use domain, only: iimin,iimax,jjmin,jjmax,HU,HV,min_depth,au,av
   use variables_2d, only: zub,zvb,ru,rv,Uinto,Vinto
   use variables_3d, only: ssuo,ssun,ssvo,ssvn
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
!  $Log: slow_bottom_friction.F90,v $
!  Revision 1.6  2006-02-10 22:41:56  hb
!  Source code documentation extended
!
!  Revision 1.5  2003-09-12 16:27:27  kbk
!  removed save attributes for local variables - now compiles using PGF
!
!  Revision 1.4  2003/08/14 10:53:23  kbk
!  need temporary velocities in some halo zones
!
!  Revision 1.3  2003/04/23 12:16:34  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 13:36:38  kbk
!  parallel support, cleaned code + NO_3D, NO_BAROCLINIC
!
!  Revision 1.1.1.1  2002/05/02 14:00:55  gotm
!  recovering after CVS crash
!
!  Revision 1.3  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.2  2001/05/03 20:12:31  bbh
!  Use of variables_3d
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
   integer                   :: i,j
   REALTYPE                  :: uloc,vloc,HH
   logical,save              :: first=.true.
   REALTYPE                  :: Ui(I2DFIELD)
   REALTYPE                  :: Vi(I2DFIELD)
   REALTYPE                  :: ruu(I2DFIELD)
   REALTYPE                  :: rvv(I2DFIELD)
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'slow_bottom_friction() # ',Ncall
#endif

   do j=jjmin,jjmax
      do i=iimin,iimax
         if(au(i,j) .ge. 1) then
            Ui(i,j)=Uinto(i,j)/(ssuo(i,j)+HU(i,j))
         else
            Ui(i,j)=_ZERO_
         end if
      end do
   end do
!  need velocity in some halo points as well
   Ui(:,jjmax+1) = Uinto(:,jjmax+1)/(ssuo(:,jjmax+1)+HU(:,jjmax+1))
   Ui(iimin-1,:) = Uinto(iimin-1,:)/(ssuo(iimin-1,:)+HU(iimin-1,:))

   do j=jjmin,jjmax
      do i=iimin,iimax
         if(av(i,j) .ge. 1) then
            Vi(i,j)=Vinto(i,j)/(ssvo(i,j)+HV(i,j))
         else
            Vi(i,j)=_ZERO_
         end if
      end do
   end do
!  need velocity in some halo points as well
   Vi(:,jjmin-1) = Vinto(:,jjmin-1)/(ssvo(:,jjmin-1)+HV(:,jjmin-1))
   Vi(iimax+1,:) = Vinto(iimax+1,:)/(ssvo(iimax+1,:)+HV(iimax+1,:))

   do j=jjmin,jjmax
      do i=iimin,iimax
         if (au(i,j) .ge. 1) then
            HH=max(min_depth,ssun(i,j)+HU(i,j))
            ruu(i,j)=(zub(i,j)+0.5*HH)/zub(i,j)
            if (ruu(i,j) .le. _ONE_) then
               STDERR i,j,ssuo(i,j),' Bottom xfriction coefficient infinite.'
               stop 'slow_bottom_friction()'
            end if
            ruu(i,j)=(kappa/log(ruu(i,j)))**2
         end if
      end do
   end do

   do j=jjmin,jjmax
      do i=iimin,iimax
         if (av(i,j) .ge. 1) then
            HH=max(min_depth,ssvn(i,j)+HV(i,j))
            rvv(i,j)=(zvb(i,j)+0.5*HH)/zvb(i,j)
            if (rvv(i,j) .le. _ONE_) then
               STDERR i,j,ssvo(i,j),' Bottom yfriction coefficient infinite.'
               stop 'slow_bottom_friction()'
            end if
            rvv(i,j)=(kappa/log(rvv(i,j)))**2
         end if
      end do
   end do

   do j=jjmin,jjmax
      do i=iimin,iimax
         if (au(i,j) .ge. 1) then
            uloc=Ui(i,j)
            vloc=0.25*( Vi(i  ,j  )    &
                       +Vi(i+1,j  )    &
                       +Vi(i  ,j-1)    &
                       +Vi(i+1,j-1) )
            ru(i,j)=ruu(i,j)*sqrt(uloc**2+vloc**2)
         else
            ru(i,j)=_ZERO_
         end if
      end do
   end do

   do j=jjmin,jjmax
      do i=iimin,iimax
         if (av(i,j) .ge. 1) then
            uloc=0.25*( Ui(i  ,j  )    &
                       +Ui(i-1,j  )    &
                       +Ui(i  ,j+1)    &
                       +Ui(i-1,j+1) )
            vloc=Vi(i,j)
            rv(i,j)=rvv(i,j)*sqrt(uloc**2+vloc**2)
         else
            rv(i,j)=_ZERO_
         end if
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving slow_bottom_friction()'
   write(debug,*)
#endif
   return
   end subroutine slow_bottom_friction
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
