!$Id: slow_terms.F90,v 1.13 2010-03-22 05:02:58 hb Exp $
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
   use m3d, only: ip_fac
   use getm_timers, only: tic, toc, TIM_SLOWTERMS
#ifndef NO_BAROCLINIC
   use variables_3d, only: idpdx,idpdy
#endif
#ifdef STRUCTURE_FRICTION
   use variables_3d, only: sf
#endif

!$ use omp_lib
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
   REALTYPE                  :: vertsum
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'slow_terms() # ',Ncall
#endif
   call tic(TIM_SLOWTERMS)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,vertsum)

   if (kmax .gt. 1) then

!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin,imax
            if (au(i,j) .ge. 1) then
               vertsum=-UEx(i,j)
               do k=1,kmax
                  if (k .ge. kumin(i,j)) then
#ifdef NO_BAROCLINIC
                     vertsum=vertsum+uuEx(i,j,k)
#else
                     vertsum=vertsum+uuEx(i,j,k)-ip_fac*idpdx(i,j,k)
#endif
                  end if
               end do
               SlUx(i,j)=vertsum
            end if
         end do
      end do
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin,imax
            if (av(i,j) .ge. 1) then
               vertsum=-VEx(i,j)
               do k=1,kmax
                  if (k .ge. kvmin(i,j)) then
#ifdef NO_BAROCLINIC
                     vertsum=vertsum+vvEx(i,j,k)
#else
                     vertsum=vertsum+vvEx(i,j,k)-ip_fac*idpdy(i,j,k)
#endif
                  end if
               end do
               SlVx(i,j)=vertsum
            end if
         end do
      end do
!$OMP END DO

   else
!
! Here kmax=1, so the loops degenerate and there is no need 
! to test for k .ge. kumin(i,j).
      k=1
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin,imax
            if (au(i,j) .ge. 1) then
#ifdef NO_BAROCLINIC
               SlUx(i,j)= _ZERO_
#else
               SlUx(i,j)=-ip_fac*idpdx(i,j,k)
#endif
            end if
         end do
      end do
!$OMP END DO

!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin,imax
            if (av(i,j) .ge. 1) then
#ifdef NO_BAROCLINIC
               SlVx(i,j)= _ZERO_
#else
               SlVx(i,j)=-ip_fac*idpdy(i,j,k)
#endif
            end if
         end do
      end do
!$OMP END DO
   endif

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (au(i,j) .ge. 1) then
            k=kumin(i,j)
            if (kmax .gt. 1) then
#ifdef NO_SLR
               STDERR 'NO_SLR U'
               Slru(i,j)= _ZERO_
#else
               Slru(i,j)=-Uint(i,j)/(_HALF_*(ssuo(i,j)+ssun(i,j))         &
                                       +HU(i,j))*ru(i,j)               &
                     +uu(i,j,k)/(_HALF_*(huo(i,j,k)+hun(i,j,k)))*rru(i,j)
#endif
#ifdef STRUCTURE_FRICTION
               do k=1,kmax
                  Slru(i,j)=Slru(i,j)+uu(i,j,k)*_HALF_*(sf(i,j,k)+sf(i+1,j,k))
               end do
#endif
            else
               Slru(i,j)= _ZERO_
            end if
         end if
      end do
   end do
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (av(i,j) .ge. 1) then
            k=kvmin(i,j)
            if (kmax .gt. 1) then
#ifdef NO_SLR
               STDERR 'NO_SLR V'
               Slrv(i,j)= _ZERO_
#else
               Slrv(i,j)=-Vint(i,j)/(_HALF_*(ssvo(i,j)+ssvn(i,j))         &
                                    +HV(i,j))*rv(i,j)                  &
                  +vv(i,j,k)/(_HALF_*(hvo(i,j,k)+hvn(i,j,k)))*rrv(i,j)
#endif
#ifdef STRUCTURE_FRICTION
               do k=1,kmax
                  Slrv(i,j)=Slrv(i,j)+vv(i,j,k)*_HALF_*(sf(i,j,k)+sf(i,j+1,k))
               end do
#endif
            else
               Slrv(i,j)=_ZERO_
            end if
         end if
      end do
   end do
!$OMP END DO
!$OMP END PARALLEL

   call toc(TIM_SLOWTERMS)
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
