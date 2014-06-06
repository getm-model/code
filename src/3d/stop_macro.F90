#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: stop_macro - terminates the macro loop \label{sec-stop-macro}
!
! !INTERFACE:
   subroutine stop_macro(reset)
!
! !DESCRIPTION:
!
! This routine should be called from {\tt m3d} at the end of each macro
! time step in order to calculate the so-called slow terms and to
! reinitialise the transports {\tt Uint} and {\tt Vint} to zero.
! The mathematical form of the interaction terms between the barotropic
! and the baroclinic mode is given by
! equations (\ref{Slowfirst}) - (\ref{Slowlast}), see section
! \ref{SectionVerticalIntegrated}.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax,au,av
   use variables_2d, only: Uint,Vint,UEx,VEx,Slru,Slrv,SlUx,SlVx,ru,rv
   use variables_3d, only: kumin,kvmin,uu,vv,hun,hvn,Dn,Dveln,Dun,Dvn
   use variables_3d, only: dt,Uadv,Vadv,uuEx,vvEx,rru,rrv
#ifndef NO_BAROCLINIC
   use variables_3d, only: idpdx,idpdy
#endif
#ifdef STRUCTURE_FRICTION
   use variables_3d, only: sf
#endif
   use m3d, only: ip_fac
   use m2d, only: uv_advect,uv_diffusion,bottom_friction
   use getm_timers, only: tic, toc, TIM_STOPMCR
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical,intent(in) :: reset
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   REALTYPE                  :: vertsum,rrdraghm1,rdragDm1
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'stop_macro() # ',Ncall
#endif
   call tic(TIM_STOPMCR)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,vertsum,rrdraghm1,rdragDm1)

   if (kmax .gt. 1) then

      call bottom_friction(Uadv,Vadv,Dun,Dvn,ru,rv)
      call uv_advect(Uadv,Vadv,Dn,Dveln,Dun,Dvn)
      call uv_diffusion(0,Uadv,Vadv,Dn,Dun,Dvn) ! Has to be called after uv_advect.

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
!              KK-TODO: adapt taubx calculation in stresses_3d and use taubx here.
!                       consider uuEx/UEx as well for "new" transport?
!                       if no special treatment of Slru is needed anymore, it can
!                       be combined with SlUx.
               k=kumin(i,j)
               rrdraghm1 = rru(i,j) / hun(i,j,k)
               rdragDm1  = ru (i,j) / Dun(i,j)
               Slru(i,j) =   rrdraghm1*uu  (i,j,k) / (_ONE_+dt*rrdraghm1) &
                           - rdragDm1 *Uadv(i,j)   / (_ONE_+dt*rdragDm1 )
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
!              KK-TODO: adapt tauby calculation in stresses_3d and use tauby here.
!                       consider vvEx/VEx as well for "new" transport?
!                       if no special treatment of Slrv is needed anymore, it can
!                       be combined with SlVx.
               k=kvmin(i,j)
               rrdraghm1 = rrv(i,j) / hvn(i,j,k)
               rdragDm1  = rv (i,j) / Dvn(i,j)
               Slrv(i,j) =   rrdraghm1*vv  (i,j,k) / (_ONE_+dt*rrdraghm1) &
                           - rdragDm1 *Vadv(i,j)   / (_ONE_+dt*rdragDm1 )
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

   if (reset) then
      Uint= _ZERO_
      Vint= _ZERO_
   end if

   call toc(TIM_STOPMCR)
#ifdef DEBUG
   write(debug,*) 'Leaving stop_macro()'
   write(debug,*)
#endif
   return
   end subroutine stop_macro
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
