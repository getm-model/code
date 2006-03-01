!$Id: slow_bottom_friction.F90,v 1.7 2006-03-01 14:45:12 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: slow_bottom_friction - slow bed friction
! \label{sec-slow-bottom-friction}
!
! !INTERFACE:
   subroutine slow_bottom_friction
!
! !DESCRIPTION:
!
! This routine basically calculates the bed friction, as it would come out
! if the vertically and macro timestep averaged velocity would be used.
! The output of this subroutine is thus $R\sqrt{u^2+v^2}$ on the U-points
! (see variable {\tt ruu}) and on the V-points (see {\tt rvv}) with the
! vertically and macro timestep averaged velocity components on the
! old time step, $u$ and $v$,
! which are in the code denoted by {\tt Ui} and {\tt Vi}, respectively.
! The drag coefficient $R$ is given by eq.\
! (\ref{bottom_vert}) on page \pageref{bottom_vert}.
! The results for the variables {\tt ruu} and {\tt rvv} will then be used
! in the routine {\tt slow\_terms} described on page \pageref{sec-slow-terms}
! for the calculation of the slow terms $S^x_F$ and $S^y_F$, see
! section \ref{SectionVerticalIntegrated}.
!
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
