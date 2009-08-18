!$Id: bottom_friction.F90,v 1.8 2009-08-18 10:24:43 bjb Exp $
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
   use parameters, only: kappa,avmmol
   use domain, only: imin,imax,jmin,jmax,au,av,min_depth
   use variables_2d
   use getm_timers, only: tic, toc, TIM_BOTTFRICT
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: runtype
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
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

   do j=jmin,jmax
      do i=imin,imax
         if (au(i,j) .gt. 0) then
            vloc(i,j)=0.25* ( V(i  ,j  )/DV(i  ,j  )   &
                             +V(i+1,j  )/DV(i+1,j  )   &
                             +V(i  ,j-1)/DV(i  ,j-1)   &
                             +V(i+1,j-1)/DV(i+1,j-1) )
         else
            vloc(i,j) = _ZERO_
         end if
      end do
   end do

!  The x-direction

#ifndef DEBUG
   where (au .gt. 0)
      uloc=U/DU
      HH=max(min_depth,DU)
      ruu=(kappa/log((zub+0.5*HH)/zub))**2
   end where
#else
   uloc=U/DU
   HH=max(min_depth,DU)
   ruu=(zub+0.5*HH)/zub

   do j=jmin,jmax
      do i=imin,imax
         if (ruu(i,j) .le. _ONE_) then
            STDERR 'bottom_friction friction coefficient infinite.'
            STDERR 'i = ',i,' j = ',j
            STDERR 'min_depth = ',min_depth,' DU = ',DU(i,j)
            stop
         end if
      end do
   end do
   where (au .gt. 0)
      ruu=(kappa/log(ruu))**2
   end where
#endif

   if (runtype .eq. 1) then
      where (au .gt. 0)
         fricvel=sqrt(ruu*(uloc**2+vloc**2))
         zub=min(HH,zub0+0.1*avmmol/max(avmmol,fricvel))
         ruu=(zub+0.5*HH)/zub
         ruu=(kappa/log(ruu))**2
      end where
   end if

   where (au .gt. 0)
      ru=ruu*sqrt(uloc**2+vloc**2)
   end where

!  The y-direction
   do j=jmin,jmax
      do i=imin,imax
         if (av(i,j) .gt. 0) then
            uloc(i,j)=0.25* ( U(i  ,j  )/DU(i  ,j  )   &
                             +U(i-1,j  )/DU(i-1,j  )   &
                             +U(i  ,j+1)/DU(i  ,j+1)   &
                             +U(i-1,j+1)/DU(i-1,j+1) )
         else
            uloc(i,j) = _ZERO_
         end if
      end do
   end do

#ifndef DEBUG
   where (av .gt. 0)
      vloc=V/DV
      HH=max(min_depth,DV)
      rvv=(kappa/log((zvb+0.5*HH)/zvb))**2
   end where
#else
   vloc=V/DV
   HH=max(min_depth,DV)
   rvv=(zvb+0.5*HH)/zvb

   do j=jmin,jmax
      do i=imin,imax
         if (rvv(i,j) .le. _ONE_) then
            STDERR 'bottom_friction friction coefficient infinite.'
            STDERR 'i = ',i,' j = ',j
            STDERR 'min_depth = ',min_depth,' DV = ',DU(i,j)
            stop
         end if
      end do
   end do

   where (av .gt. 0)
      rvv=(kappa/log(rvv))**2
   end where
#endif

   if (runtype .eq. 1) then
      where (av .gt. 0)
         fricvel=sqrt(rvv*(uloc**2+vloc**2))
         zvb=min(HH,zvb0+0.1*avmmol/max(avmmol,fricvel))
         rvv=(zvb+0.5*HH)/zvb
         rvv=(kappa/log(rvv))**2
      end where
   end if

   where (av .gt. 0)
      rv=rvv*sqrt(uloc**2+vloc**2)
   end where

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
