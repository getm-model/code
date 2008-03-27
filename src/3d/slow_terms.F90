!$Id: slow_terms.F90,v 1.10 2008-03-27 08:51:23 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: slow_terms - calculation of slow terms \label{sec-slow-terms}
!
! !INTERFACE:
   subroutine slow_terms
!
! !DESCRIPTION:
!
! Here, the calculation of the so-called slow terms (which are the 
! interaction terms between the barotropic and the baroclinic mode) is
! completed. The mathematical form of these slow terms is given by 
! equations (\ref{Slowfirst}) - (\ref{Slowlast}), see section
! \ref{SectionVerticalIntegrated}.
! These calculations have been prepared in the routines 
! {\tt slow\_bottom\_friction}, {\tt slow\_advection} and
! {\tt slow\_diffusion}.

!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax,HU,HV,au,av
   use variables_2d, only: Uint,Vint,UEx,VEx,Slru,Slrv,SlUx,SlVx,ru,rv
   use variables_3d, only: kumin,kvmin,uu,vv,huo,hun,hvo,hvn
   use variables_3d, only: ssuo,ssun,ssvo,ssvn,uuEx,vvEx,rru,rrv
#ifndef NO_BAROCLINIC
   use variables_3d, only: idpdx,idpdy
#endif
#ifdef STRUCTURE_FRICTION
   use variables_3d, only: sf
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
! !LOCAL VARIABLES:
   integer                   :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'slow_terms() # ',Ncall
#endif

   if (kmax .gt. 1) then

      do j=jmin,jmax
         do i=imin,imax
            if (au(i,j) .ge. 1) then
               SlUx(i,j)=-UEx(i,j)
            end if
         end do
      end do

      do k=1,kmax
         do j=jmin,jmax
            do i=imin,imax
               if (au(i,j) .ge. 1) then
                  if (k .ge. kumin(i,j)) then
#ifdef NO_BAROCLINIC
                     SlUx(i,j)=SlUx(i,j)+uuEx(i,j,k)
#else
                     SlUx(i,j)=SlUx(i,j)+uuEx(i,j,k)-idpdx(i,j,k)
#endif
                  end if
               end if
            end do
         end do
      end do

      do j=jmin,jmax
         do i=imin,imax
            if (av(i,j) .ge. 1) then
               SlVx(i,j)=-VEx(i,j)
            end if
         end do
      end do

      do k=1,kmax
         do j=jmin,jmax
            do i=imin,imax
               if (av(i,j) .ge. 1) then
                  if (k .ge. kvmin(i,j)) then
#ifdef NO_BAROCLINIC
                     SlVx(i,j)=SlVx(i,j)+vvEx(i,j,k)
#else
                     SlVx(i,j)=SlVx(i,j)+vvEx(i,j,k)-idpdy(i,j,k)
#endif
                  end if
               end if
            end do
         end do
      end do

   else

      do k=1,kmax
         do j=jmin,jmax
            do i=imin,imax
               if (au(i,j) .ge. 1) then
                  if (k .ge. kumin(i,j)) then
#ifdef NO_BAROCLINIC
                     SlUx(i,j)= _ZERO_
#else
                     SlUx(i,j)=-idpdx(i,j,k)
#endif
                  end if
               end if
            end do
         end do
      end do

      do k=1,kmax
         do j=jmin,jmax
            do i=imin,imax
               if (av(i,j) .ge. 1) then
                  if (k.ge.kvmin(i,j)) then
#ifdef NO_BAROCLINIC
                     SlVx(i,j)= _ZERO_
#else
                     SlVx(i,j)=-idpdy(i,j,k)
#endif
                  end if
               end if
            end do
         end do
      end do
   endif

   do j=jmin,jmax
      do i=imin,imax
         if (au(i,j) .ge. 1) then
            k=kumin(i,j)
            if (kmax .gt. 1) then
#ifdef NO_SLR
               STDERR 'NO_SLR U'
               Slru(i,j)= _ZERO_
#else
               Slru(i,j)=-Uint(i,j)/(0.5*(ssuo(i,j)+ssun(i,j))         &
                                       +HU(i,j))*ru(i,j)               &
                     +uu(i,j,k)/(0.5*(huo(i,j,k)+hun(i,j,k)))*rru(i,j)
#endif
#ifdef STRUCTURE_FRICTION
               do k=1,kmax
                  Slru(i,j)=Slru(i,j)+uu(i,j,k)*0.5*(sf(i,j,k)+sf(i+1,j,k))
               end do
#endif
            else
               Slru(i,j)= _ZERO_
            end if
         end if
      end do
   end do

   do j=jmin,jmax
      do i=imin,imax
         if (av(i,j) .ge. 1) then
            k=kvmin(i,j)
            if (kmax .gt. 1) then
#ifdef NO_SLR
               STDERR 'NO_SLR V'
               Slrv(i,j)= _ZERO_
#else
               Slrv(i,j)=-Vint(i,j)/(0.5*(ssvo(i,j)+ssvn(i,j))         &
                                    +HV(i,j))*rv(i,j)                  &
                  +vv(i,j,k)/(0.5*(hvo(i,j,k)+hvn(i,j,k)))*rrv(i,j)
#endif
#ifdef STRUCTURE_FRICTION
               do k=1,kmax
                  Slrv(i,j)=Slrv(i,j)+vv(i,j,k)*0.5*(sf(i,j,k)+sf(i,j+1,k))
               end do
#endif
            else
               Slrv(i,j)=_ZERO_
            end if
         end if
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving slow_terms()'
   write(debug,*)
#endif
   return
   end subroutine slow_terms
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
