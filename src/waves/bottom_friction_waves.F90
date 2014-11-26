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
   use parameters, only: kappa
   use domain, only: imin,imax,jmin,jmax,az,au,av
   use domain, only: z0,zub0,zvb0
   use variables_waves, only: waveH,waveT,waveK,coswavedir,sinwavedir
   use waves, only: waves_bbl_method,NO_WBBL,WBBL_DATA2,WBBL_SOULSBY05
   use waves, only: wbbl_tauw,wbbl_rdrag
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
   REALTYPE                                 :: tauw,wbl
   REALTYPE                                 :: cdm1,cosangle
   REALTYPE                                 :: ttransx,ttransy,ttrans
   REALTYPE                                 :: tauc,taubm,taubc,taube,taubp
   integer                                  :: i,j,rc
   logical                                  :: calc_taubmax

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

   calc_taubmax = .false.
   if (present(taubmax)) then
      calc_taubmax = associated(taubmax)
   end if

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j)                                         &
!$OMP          PRIVATE(i,tauwwbl,cdm1,cosangle,ttransx,ttransy,ttrans) &
!$OMP          PRIVATE(tauc,tauw,taubm,taubc,taube,taubp)

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
   do j=jmin-HALO,jmax+HALO
#endif
      do i=imin-HALO,imax+HALO
         if (az(i,j) .ne. 0) then
            taubw(i,j) = wbbl_tauw(waveT(i,j),waveH(i,j),waveK(i,j), &
                                   Dvel(i,j),z0(i,j),wbbl(i,j))
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
                  cosangle = (coswavedir(i,j)*ttransx + sinwavedir(i,j)*ttransy ) / ttrans
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
