#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: bottom_friction_3d - bottom friction
! \label{sec-bottom-friction-3d}
!
! !INTERFACE:
   subroutine bottom_friction_3d
!
! !DESCRIPTION:
!
! Based on the assumption that the velocity distribution in the bottom
! layer is logarithmic,
! the product of the drag coefficient with the
! absolute value of the current speed in the bottom layer,
!
! \begin{equation}
! r \sqrt{u_b^2+v_b^2}
! \end{equation}
!
! with the velocity components of the bottom layer, $u_b$ and $v_b$,
! and the drag coefficient
!
! \begin{equation}\label{r}
! r = \left(\frac{\kappa}{\ln \left(\frac{0.5h_1+z_0^b}{z_0^b}\right)}
! \right)^2,
! \end{equation}
!
! is calculated and
! provided as output parameters {\tt rru} (for U-points) and
! {\tt rrv} (for V-points). The layer height $h_1$ in (\ref{r}) is set to
! the thickness of the bottom layer in the respective U- or V-point.
!
! There are some experimental options for the interested user included
! here. It is possible to change the interpolation of $u$ to V-points
! and of $v$ to U-points from velocity-based interpolation (as done
! presently) to transport-based averaging (commented out). Furthermore,
! the user may activate some outcommented lines which allow the
! consideration of flow-depending bottom roughness length $z_0^b$
! according to (\ref{Defz0b}), see page \pageref{Defz0b}.
!
! For a derivation of (\ref{r}), see section \ref{SectionBedFric} on
! page \pageref{SectionBedFric}.
!
! !USES:
   use parameters, only: kappa
   use m2d, only: avmmol
   use domain, only: imin,imax,jmin,jmax,kmax,au,av,min_depth
   use variables_2d, only: zub,zvb,zub0,zvb0
   use variables_3d, only: kumin,kvmin,uu,vv,huo,hun,hvo,hvn,rru,rrv
   use getm_timers, only: tic, toc, TIM_BOTTFRICT3D
!$ use omp_lib
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: i,j,kk
   REALTYPE                  :: r,hh,fricvel
   logical, save             :: first=.true.
   REALTYPE                  :: uuloc(I2DFIELD)
   REALTYPE                  :: uvloc(I2DFIELD)
   REALTYPE                  :: vuloc(I2DFIELD)
   REALTYPE                  :: vvloc(I2DFIELD)
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'bottom_friction_3d() # ',Ncall
#endif
   call tic(TIM_BOTTFRICT3D)

#if 0
   if(first) then
      uuloc = _ZERO_
      uvloc = _ZERO_
      vuloc = _ZERO_
      vvloc = _ZERO_
      first = .false.
   end if
#endif

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,kk,r,hh,fricvel)
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
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
!$OMP END DO

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
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
!$OMP END DO

#if 1
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
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
!$OMP END DO
#endif

#if 1
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
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
!$OMP END DO
#endif

!$OMP END PARALLEL

#ifdef SLICE_MODEL
   do i=imin,imax
      rrv(i,1)=rrv(i,2)
      rrv(i,3)=rrv(i,2)
   enddo
#endif

   call toc(TIM_BOTTFRICT3D)
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
