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
   use domain, only: imin,imax,jmin,jmax,kmax,au,av,z0_method
   use variables_2d, only: Uint,Vint,UEx,VEx,Slru,Slrv,SlUx,SlVx,ru,rv
   use variables_3d, only: kumin,kvmin,uu,vv,hun,hvn,Dn,Dun,Dvn
   use variables_3d, only: uuEx,vvEx,rru,rrv
#ifndef NO_BAROCLINIC
   use variables_3d, only: idpdx,idpdy
#endif
#ifdef STRUCTURE_FRICTION
   use variables_3d, only: sf
#endif
   use m3d, only: ip_fac
   use m2d_general, only: bottom_friction,calc_uvex
   use getm_timers, only: tic, toc, TIM_SLOWTERMS
!$ use omp_lib
   IMPLICIT NONE
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

      if (z0_method .ne. 0) then
         call bottom_friction(Uint,Vint,Dun,Dvn,ru,rv)
      end if
      call calc_uvex(0,Uint,Vint,Dn,Dun,Dvn)

!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin,imax

            if (au(i,j) .ge. 1) then

               vertsum=-UEx(i,j)
               do k=kumin(i,j),kmax
#ifdef NO_BAROCLINIC
                     vertsum=vertsum+uuEx(i,j,k)
#else
                     vertsum=vertsum+uuEx(i,j,k)-ip_fac*idpdx(i,j,k)
#endif
               end do
               SlUx(i,j)=vertsum

#ifdef NO_SLR
               STDERR 'NO_SLR U'
               Slru(i,j)= _ZERO_
#else
               k=kumin(i,j)
               Slru(i,j) =   rru(i,j)*uu(i,j,k)/hun(i,j,k) &
                           - ru(i,j)*Uint(i,j)/Dun(i,j)
#endif

#ifdef STRUCTURE_FRICTION
               do k=kumin(i,j),kmax
                  Slru(i,j)=Slru(i,j)+uu(i,j,k)*_HALF_*(sf(i,j,k)+sf(i+1,j,k))
               end do
#endif

            end if

            if (av(i,j) .ge. 1) then

               vertsum=-VEx(i,j)
               do k=kvmin(i,j),kmax
#ifdef NO_BAROCLINIC
                     vertsum=vertsum+vvEx(i,j,k)
#else
                     vertsum=vertsum+vvEx(i,j,k)-ip_fac*idpdy(i,j,k)
#endif
               end do
               SlVx(i,j)=vertsum

#ifdef NO_SLR
               STDERR 'NO_SLR V'
               Slrv(i,j)= _ZERO_
#else
               k=kvmin(i,j)
               Slrv(i,j) =   rrv(i,j)*vv(i,j,k)/hvn(i,j,k) &
                           - rv(i,j)*Vint(i,j)/Dvn(i,j)
#endif

#ifdef STRUCTURE_FRICTION
               do k=kvmin(i,j),kmax
                  Slrv(i,j)=Slrv(i,j)+vv(i,j,k)*_HALF_*(sf(i,j,k)+sf(i,j+1,k))
               end do
#endif

            end if

         end do
      end do
!$OMP END DO

#ifndef NO_BAROCLINIC

   else
!
! Here kmax=1, so the loops degenerate and there is no need
! to test for k .ge. kumin(i,j).
      k=1
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin,imax
            if (au(i,j) .ge. 1) then
               SlUx(i,j)=-ip_fac*idpdx(i,j,k)
            end if
            if (av(i,j) .ge. 1) then
               SlVx(i,j)=-ip_fac*idpdy(i,j,k)
            end if
         end do
      end do
!$OMP END DO

#endif

   end if

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
