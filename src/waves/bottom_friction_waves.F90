#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: bottom_friction_waves - calculates the 2D bottom friction.
!
! !INTERFACE:
   subroutine bottom_friction_waves(U1,V1,DU1,DV1,Dvel,u_vel,v_vel,velU,velV,ru,rv,zub,zvb,taubmax)
!
! !DESCRIPTION:
!
! !USES:
   use parameters, only: kappa,avmmol
   use domain, only: imin,imax,jmin,jmax,az,au,av
   use domain, only: z0,zub0,zvb0
   use variables_waves, only: waveH,waveT,waveK,coswavedir,sinwavedir
   use waves, only: waves_bbl_method,NO_WBBL,WBBL_DATA2,WBBL_SOULSBY05
   use waves, only: wbbl_rdrag
   use getm_timers, only: tic,toc,TIM_WAVES
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT VARIABLES:
   REALTYPE,dimension(E2DFIELD),intent(in)    :: U1,V1,DU1,DV1,Dvel,u_vel,v_vel,velU,velV
!
! !INPUT/OUTPUT VARIABLES:
   REALTYPE,dimension(E2DFIELD),intent(inout) :: ru,rv,zub,zvb
   REALTYPE,dimension(:,:),pointer,intent(inout),optional :: taubmax
!
! !OUTPUT VARIABLES:
!
! !REVISION HISTORY:
!  Original author(s): Ulf Gräwe
!                      Saeed Moghimi
!                      Knut Klingbeil
!
!  !LOCAL VARIABLES:
   REALTYPE,dimension(E2DFIELD)             :: taubw,wbbl,taubmx,taubmy
   REALTYPE                                 :: Hrms,omegam1,uorb,aorb
   REALTYPE                                 :: Rew,fwr,fws,fwl,fw
   REALTYPE                                 :: tauwr,tauws,tauwl,tauw
   REALTYPE                                 :: wbl
   REALTYPE                                 :: cdm1,cosangle
   REALTYPE                                 :: ttransx,ttransy,ttrans
   REALTYPE                                 :: coscurrdir,sincurrdir
   REALTYPE                                 :: tauc,taubm,taubc,taube,taubp
   integer                                  :: i,j,rc
   logical                                  :: calc_taubmax
   logical,save                             :: first=.true.
   REALTYPE,save      :: avmmolm1
   REALTYPE,parameter :: sqrthalf=sqrt(_HALF_)
   REALTYPE,parameter :: pi=3.1415926535897932384626433832795029d0
   REALTYPE,parameter :: oneovertwopi=_HALF_/pi
   REALTYPE,parameter :: ar = 0.24d0 ! 0.26d0
   REALTYPE,parameter :: as = 0.24d0 ! 0.22d0
   REALTYPE,parameter :: Rew_crit = 5.0d5 ! (Stanev et al., 2009)
!   REALTYPE,parameter :: Rew_crit = 1.5d5 ! (Soulsby & Clarke, 2005)
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'bottom_friction_waves() # ',Ncall
#endif
#ifdef SLICE_MODEL
   j = jmax/2 ! this MUST NOT be changed!!!
#endif

   if (waves_bbl_method .eq. NO_WBBL) return

   call tic(TIM_WAVES)

   if (first) then
      avmmolm1 = _ONE_ / avmmol
      first = .false.
   end if

   calc_taubmax = .false.
   if (present(taubmax)) then
      calc_taubmax = associated(taubmax)
   end if

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j)                                         &
!$OMP          PRIVATE(i,Hrms,omegam1,uorb,aorb,Rew,fwr,fws,fwl,fw)    &
!$OMP          PRIVATE(tauw,tauwr,tauws,tauwl)                         &
!$OMP          PRIVATE(wbl,cdm1,cosangle,coscurrdir,sincurrdir)        &
!$OMP          PRIVATE(ttransx,ttransy,ttrans)                         &
!$OMP          PRIVATE(tauc,tauw,taubm,taubc,taube,taubp)

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
   do j=jmin-HALO,jmax+HALO
#endif
      do i=imin-HALO,imax+HALO
         if (az(i,j) .ne. 0) then
            Hrms = sqrthalf * waveH(i,j)
            omegam1 = oneovertwopi * waveT(i,j)
!           wave orbital velocity amplitude at bottom (peak orbital velocity, ubot in SWAN)
            uorb = _HALF_ * Hrms / ( omegam1*sinh(waveK(i,j)*Dvel(i,j)) )
!           wave orbital excursion
            aorb = omegam1 * uorb
!           wave Reynolds number
            Rew = aorb * uorb * avmmolm1

!           Note (KK): We do not calculate fw alone, because for small
!                      uorb this can become infinite.

!           KK-TODO: For combined wave-current flow, the decision on
!                    turbulent or laminar flow depends on Rew AND Rec!
!                    (Soulsby & Clarke, 2005)
!                    However, here we decide according to Lettmann et al. (2009).
!                    (Or do we want to assume always turbulent currents?)
            if ( Rew .gt. Rew_crit ) then
!              wave friction factor for rough turbulent flow
!               fwr = 1.39d0 * (aorb/z0(i,j))**(-0.52d0)
               tauwr = _HALF_ * 1.39d0 * (omegam1/z0(i,j))**(-0.52d0) * uorb**(2-0.52d0)
!              wave friction factor for smooth turbulent flow
!               fws = 0.0521d0 * Rew**(-0.187d0)
               tauws = _HALF_ * (omegam1*avmmolm1)**(-0.187d0) * uorb**(2-2*0.187d0)
!              wave friction factor
!              Note (KK): For combined wave-current flow, the decision on
!                         rough or smooth flow depends on the final taubmax.
!                         (Soulsby & Clarke, 2005)
!                         However, here we decide according to Stanev et al. (2009).
!                         (as for wave-only flow)
!               fw = max( fwr , fws )
               tauw = max( tauwr , tauws )
            else
!              wave friction factor for laminar flow
!               fwl = _TWO_ * Rew**(-_HALF_)
!               fw = fwl
               tauwl = (omegam1*avmmolm1)**(-_HALF_) * uorb
               tauw = tauwl
            end if

!           wave-only bottom stress
!            taubw(i,j) = _HALF_ * fw * uorb**2
            taubw(i,j) = tauw

!           bbl thickness (Soulsby & Clarke, 2005)
            wbbl(i,j) = max( 12.0d0*z0(i,j) , ar*omegam1*sqrt(taubw(i,j)) )
         end if
      end do
#ifndef SLICE_MODEL
   end do
#endif
!$OMP END DO

#ifdef SLICE_MODEL
!$OMP SINGLE
   taubw(:,j+1) = taubw(:,j)
   wbbl (:,j+1) = wbbl (:,j)
!$OMP END SINGLE
#endif

!  The x-direction

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
   do j=jmin-HALO+1,jmax+HALO-1
#endif
      do i=imin-HALO,imax+HALO-1
!        velU must be zero at land!!!
         if (velU(i,j) .gt. _ZERO_) then
            tauc = ru(i,j) * velU(i,j)
            tauw = _HALF_ * ( taubw(i,j) + taubw(i+1,j) )
            wbl = _HALF_ * ( wbbl(i,j) + wbbl(i+1,j) )
            ru(i,j) = wbbl_rdrag(tauc,tauw,ru(i,j),velU(i,j),DU1(i,j),wbl,zub(i,j))
            cdm1 = velU(i,j) / ru(i,j)
            zub(i,j) = _HALF_*DU1(i,j) / (exp(kappa*sqrt(cdm1)) - _ONE_)
         end if
      end do
#ifndef SLICE_MODEL
   end do
#endif
!$OMP END DO


!  The y-direction

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
   do j=jmin-HALO,jmax+HALO-1
#endif
      do i=imin-HALO+1,imax+HALO-1
!        velV must be zero at land!!!
         if (velV(i,j) .gt. _ZERO_) then
            tauc = rv(i,j) * velV(i,j)
            tauw = _HALF_ * ( taubw(i,j) + taubw(i,j+1) )
            wbl = _HALF_ * ( wbbl(i,j) + wbbl(i,j+1) )
            rv(i,j) = wbbl_rdrag(tauc,tauw,rv(i,j),velV(i,j),DV1(i,j),wbl,zvb(i,j))
            cdm1 = velV(i,j) / rv(i,j)
            zvb(i,j) = _HALF_*DV1(i,j) / (exp(kappa*sqrt(cdm1)) - _ONE_)
         end if
      end do
#ifndef SLICE_MODEL
   end do
#endif
!$OMP END DO

#ifdef SLICE_MODEL
!$OMP SINGLE
   ru (imin-HALO  :imax+HALO-1,j+1) = ru (imin-HALO  :imax+HALO-1,j)
   rv (imin-HALO+1:imax+HALO-1,j-1) = rv (imin-HALO+1:imax+HALO-1,j)
   rv (imin-HALO+1:imax+HALO-1,j+1) = rv (imin-HALO+1:imax+HALO-1,j)
   zub(imin-HALO  :imax+HALO-1,j+1) = zub(imin-HALO  :imax+HALO-1,j)
   zvb(imin-HALO+1:imax+HALO-1,j-1) = zvb(imin-HALO+1:imax+HALO-1,j)
   zvb(imin-HALO+1:imax+HALO-1,j+1) = zvb(imin-HALO+1:imax+HALO-1,j)
!$OMP END SINGLE
#endif


   if (calc_taubmax) then

!$OMP WORKSHARE
!     velocities must be zero at land!!!
      taubmx = ru * u_vel
      taubmy = rv * v_vel
!$OMP END WORKSHARE

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin-HALO+1,jmax+HALO-1
#endif
         do i=imin-HALO+1,imax+HALO-1
            if (az(i,j) .ne. 0) then
               taubm = _HALF_*sqrt(  ( taubmx(i-1,j  ) + taubmx(i,j) )**2 &
                                   + ( taubmy(i  ,j-1) + taubmy(i,j) )**2 )

               select case(waves_bbl_method)
                  case (WBBL_DATA2)
                     taubp = taubw(i,j)
                  case (WBBL_SOULSBY05)
                     taubc = taubmax(i,j)
                     taube = sqrt( taubc**2 + taubw(i,j)**2 )
                     taubp = sqrt( taube * taubw(i,j))
               end select

!              we need the cosine of the angle between waves and currents
               ttransx = U1(i-1,j) + U1(i,j)
               ttransy = V1(i,j-1) + V1(i,j)
               ttrans  = sqrt( ttransx**2 + ttransy**2 )
               if (ttrans .gt. _ZERO_) then
                  coscurrdir = ttransx / ttrans
                  sincurrdir = ttransy / ttrans
                  cosangle = coscurrdir*coswavedir(i,j) + sincurrdir*sinwavedir(i,j)
               else
                  cosangle = _ZERO_
               end if

               taubmax(i,j) = sqrt( taubm**2 + taubp**2 + _TWO_*taubm*taubp*cosangle )
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO

#ifdef SLICE_MODEL
!$OMP SINGLE
      taubmax(imin-HALO+1:imax+HALO-1,j+1) = taubmax(imin-HALO+1:imax+HALO-1,j)
!$OMP END SINGLE
#endif

   end if


!$OMP END PARALLEL

   call toc(TIM_WAVES)

#ifdef DEBUG
   write(debug,*) 'Leaving bottom_friction_waves()'
   write(debug,*)
#endif
   return
   end subroutine bottom_friction_waves
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2014 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
