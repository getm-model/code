#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: bottom_friction - calculates the 2D bottom friction.
!
! !INTERFACE:
   subroutine bottom_friction(runtype)
!
! !DESCRIPTION:
!
! In this routine the bottom friction for the external (vertically integrated)
! mode is calculated. This is done separately for the $U$-equation in the
! U-points and for the $V$-equation in the V-points.
! The drag coefficient $R$ for the external mode is given in eq.\
! (\ref{bottom_vert}) on page \pageref{bottom_vert}.
! For {\tt runtype=1} (only vertically integrated calculations), the
! bottom roughness length is depending on the bed friction
! velocity $u_*^b$ and the molecular viscosity $\nu$:
!
! \begin{equation}\label{Defz0b}
! z_0^b = 0.1 \frac{\nu}{u_*^b} + \left(z^b_0\right)_{\min},
! \end{equation}
!
! see e.g.\ \cite{KAGAN95}, i.e.\ the given roughness may be increased
! by viscous effects.
! After this, the drag coefficient is multiplied by the absolute value of the
! local velocity, which is alculated by dividing the local transports by the
! local water depths and by properly interpolating these velocities
! to the U- and V-points. The resulting fields are {\tt ru}, representing
! $R\sqrt{u^2+v^2}$ on the U-points and {\tt rv}, representing
! this quantity on the V-points.
!
! !USES:
   use parameters, only: kappa
   use m2d, only: avmmol
   use domain, only: imin,imax,jmin,jmax,au,av,min_depth
   use variables_2d
   use getm_timers,  only: tic, toc, TIM_BOTTFRICT
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: runtype
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  !LOCAL VARIABLES:
   integer                   :: i,j
   REALTYPE                  :: uloc(E2DFIELD),vloc(E2DFIELD)
   REALTYPE                  :: HH(E2DFIELD),fricvel(E2DFIELD)
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'bottom_friction() # ',Ncall
#endif
   CALL tic(TIM_BOTTFRICT)

#ifdef DEBUG
   if(Ncall .eq. 1) then
      STDERR 'bottom_friction(): checking for 0 depth values'
      do j=jmin,jmax
         do i=imin,imax
            if (av(i,j) .eq. 1) then
               if(DU(i,j) .eq. _ZERO_) then
                  STDERR 'DU=0 uv_advect ',i,j,i,j,av(i,j)
               end if
               if(DU(i-1,j) .eq. _ZERO_) then
                  STDERR 'DU=0 uv_advect ',i-1,j,i,j,av(i,j)
               end if
               if(DU(i,j+1) .eq. _ZERO_) then
                  STDERR 'DU=0 uv_advect ',i,j+1,i,j,av(i,j)
               end if
               if(DU(i-1,j+1) .eq. _ZERO_) then
                  STDERR 'DU=0 uv_advect ',i-1,j+1,i,j,av(i,j)
               end if
            end if
         end do
      end do
      do j=jmin,jmax
         do i=imin,imax
            if (au(i,j) .eq. 1) then
               if(DV(i,j) .eq. _ZERO_) then
                  STDERR 'DV=0 uv_advect ',i,j,i,j,au(i,j)
               end if
               if(DV(i+1,j) .eq. _ZERO_) then
                  STDERR 'DV=0 uv_advect ',i+1,j,i,j,au(i,j)
               end if
               if(DV(i,j-1) .eq. _ZERO_) then
                  STDERR 'DV=0 uv_advect ',i,j-1,i,j,au(i,j)
               end if
               if(DV(i+1,j-1) .eq. _ZERO_) then
                  STDERR 'DV=0 uv_advect ',i+1,j-1,i,j,au(i,j)
               end if
            end if
         end do
      end do
   end if
#endif


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)

!  The x-direction

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (au(i,j) .gt. 0) then
            vloc(i,j)=_QUART_* ( V(i  ,j  )/DV(i  ,j  )   &
                                +V(i+1,j  )/DV(i+1,j  )   &
                                +V(i  ,j-1)/DV(i  ,j-1)   &
                                +V(i+1,j-1)/DV(i+1,j-1) )
         else
            vloc(i,j) = _ZERO_
         end if
      end do
   end do
!OMP END DO NOWAIT

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO
         if (au(i,j) .gt. 0) then
            uloc(i,j) = U(i,j)/DU(i,j)
            HH(i,j)=max(min_depth,DU(i,j))
#ifndef DEBUG
            ruu(i,j)=(kappa/log((zub(i,j)+_HALF_*HH(i,j))/zub(i,j)))**2
#else
            ruu(i,j)=(zub(i,j)+_HALF_*HH(i,j))/zub(i,j)
            if (ruu(i,j) .le. _ONE_) then
!$OMP CRITICAL
               STDERR 'bottom_friction friction coefficient infinite.'
               STDERR 'i = ',i,' j = ',j
               STDERR 'min_depth = ',min_depth,' DU = ',DU(i,j)
!$OMP END CRITICAL
               stop
            end if
            ruu(i,j)=(kappa/log(ruu(i,j)))**2
#endif
         end if
      end do
   end do
!$OMP END DO

   if (runtype .eq. 1) then
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO
            if (au(i,j) .gt. 0) then
               fricvel(i,j)=sqrt(ruu(i,j)*(uloc(i,j)**2+vloc(i,j)**2))
               zub(i,j)=min(HH(i,j),zub0(i,j)+_TENTH_*avmmol/max(avmmol,fricvel(i,j)))
               ruu(i,j)=(zub(i,j)+_HALF_*HH(i,j))/zub(i,j)
               ruu(i,j)=(kappa/log(ruu(i,j)))**2
               ru(i,j)=ruu(i,j)*sqrt(uloc(i,j)**2+vloc(i,j)**2)
            end if
         end do
      end do
!$OMP END DO
   else
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO
            if (au(i,j) .gt. 0) then
               ru(i,j)=ruu(i,j)*sqrt(uloc(i,j)**2+vloc(i,j)**2)
            end if
         end do
      end do
!$OMP END DO
   end if

!  The y-direction
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (av(i,j) .gt. 0) then
            uloc(i,j)=_QUART_* ( U(i  ,j  )/DU(i  ,j  )   &
                                +U(i-1,j  )/DU(i-1,j  )   &
                                +U(i  ,j+1)/DU(i  ,j+1)   &
                                +U(i-1,j+1)/DU(i-1,j+1) )
         else
            uloc(i,j) = _ZERO_
         end if
      end do
   end do
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO
         if (av(i,j) .gt. 0) then
            vloc(i,j)=V(i,j)/DV(i,j)
            HH(i,j)=max(min_depth,DV(i,j))
#ifndef DEBUG
            rvv(i,j)=(kappa/log((zvb(i,j)+_HALF_*HH(i,j))/zvb(i,j)))**2
#else
            rvv(i,j)=(zvb(i,j)+_HALF_*HH(i,j))/zvb(i,j)
            if (rvv(i,j) .le. _ONE_) then
!$OMP CRITICAL
               STDERR 'bottom_friction friction coefficient infinite.'
               STDERR 'i = ',i,' j = ',j
               STDERR 'min_depth = ',min_depth,' DV = ',DU(i,j)
!$OMP END CRITICAL
               stop
            end if
            rvv(i,j)=(kappa/log(rvv(i,j)))**2
#endif
         end if
      end do
   end do
!$OMP END DO

   if (runtype .eq. 1) then
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO
            if (av(i,j) .gt. 0) then
               fricvel(i,j)=sqrt(rvv(i,j)*(uloc(i,j)**2+vloc(i,j)**2))
               zvb(i,j)=min(HH(i,j),zvb0(i,j)+_TENTH_*avmmol/max(avmmol,fricvel(i,j)))
               rvv(i,j)=(zvb(i,j)+_HALF_*HH(i,j))/zvb(i,j)
               rvv(i,j)=(kappa/log(rvv(i,j)))**2
               rv(i,j)=rvv(i,j)*sqrt(uloc(i,j)**2+vloc(i,j)**2)
            end if
         end do
      end do
!$OMP END DO
   else
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO
            if (av(i,j) .gt. 0) then
               rv(i,j)=rvv(i,j)*sqrt(uloc(i,j)**2+vloc(i,j)**2)
            end if
         end do
      end do
!$OMP END DO
   end if

!$OMP END PARALLEL

   CALL toc(TIM_BOTTFRICT)
#ifdef DEBUG
   write(debug,*) 'Leaving bottom_friction()'
   write(debug,*)
#endif
   return
   end subroutine bottom_friction
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
