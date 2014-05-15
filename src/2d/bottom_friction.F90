#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: bottom_friction - calculates the 2D bottom friction.
!
! !INTERFACE:
   subroutine bottom_friction(U1,V1,DU1,DV1,Dvel,ru,rv,kwe,zub,zvb,taubmax)
!  Note (KK): keep in sync with interface in m2d.F90
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
   use domain, only: imin,imax,jmin,jmax,az,au,av
   use domain, only: bottfric_method,cd_min,z0d_iters,zub0,zvb0
   use waves, only: waves_method,NO_WAVES,bottom_friction_waves
   use getm_timers, only: tic,toc,TIM_BOTTFRIC
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(in) :: U1,V1,DU1,DV1,Dvel
   logical,intent(in),optional             :: kwe !keyword-enforcer
!
! !OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(out)                 :: ru,rv
   REALTYPE,dimension(E2DFIELD),intent(out),target,optional :: zub,zvb
   REALTYPE,dimension(:,:),pointer,intent(out),optional     :: taubmax
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  !LOCAL VARIABLES:
   REALTYPE,dimension(E2DFIELD)             :: work2d,velU,velV
   REALTYPE,dimension(E2DFIELD),target      :: t_zub,t_zvb
   REALTYPE,dimension(:,:),allocatable,save :: u_vel,v_vel
   REALTYPE,dimension(:,:),pointer          :: p_zub,p_zvb
   REALTYPE                                 :: vel,cd,cdsqrt,z0d
   integer                                  :: i,j,it,rc
   logical,save                             :: first=.true.
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'bottom_friction() # ',Ncall
#endif
#ifdef SLICE_MODEL
   j = jmax/2 ! this MUST NOT be changed!!!
#endif
   CALL tic(TIM_BOTTFRIC)

   if (bottfric_method.eq.2 .or. bottfric_method.eq.3) then

      if (first) then
         allocate(u_vel(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_2d: Error allocating memory (u_vel)'
         u_vel=_ZERO_

         allocate(v_vel(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_2d: Error allocating memory (v_vel)'
         v_vel=_ZERO_

         first = .false.
      end if

      if (present(zub)) then
         p_zub => zub
      else
         p_zub => t_zub
      end if
      if (present(zvb)) then
         p_zvb => zvb
      else
         p_zvb => t_zvb
      end if

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j)                                         &
!$OMP          PRIVATE(i,vel,cd,cdsqrt,it,z0d)


!     KK-TODO: the present implementation sets normal velocity outside open
!              bdy cell to zero (we need proper mirror)

!     zonal velocity
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin-HALO,jmax+HALO
#endif
         do i=imin-HALO,imax+HALO-1
            if (au(i,j) .ge. 1) then
               u_vel(i,j) = U1(i,j)/DU1(i,j)
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO NOWAIT

!     meridional velocity
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin-HALO,jmax+HALO-1
#endif
         do i=imin-HALO,imax+HALO
            if (av(i,j) .ge. 1) then
               v_vel(i,j) = V1(i,j)/DV1(i,j)
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO

#ifdef SLICE_MODEL
!$OMP SINGLE
      v_vel(:,j-1) = v_vel(:,j)
!$OMP END SINGLE
#endif

!     The x-direction

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin-HALO+1,jmax+HALO-1
#endif
!        calculate v_velC
         do i=imin-HALO,imax+HALO
            if (az(i,j) .ge. 1) then
               work2d(i,j) = _HALF_ * ( v_vel(i,j-1) + v_vel(i,j) )
            end if
         end do
#ifdef SLICE_MODEL
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME)
#endif
         do i=imin-HALO,imax+HALO-1
            if ( au(i,j) .ge. 1 ) then
               velU(i,j) = sqrt( u_vel(i,j)**2 + (_HALF_*(work2d(i,j)+work2d(i+1,j)))**2 )
               z0d = zub0(i,j)
!              Note (KK): note shifting of log profile so that U(-H)=0
               cdsqrt = kappa / log( _ONE_ + _HALF_*DU1(i,j)/z0d )
               if (avmmol.gt._ZERO_ .and. velU(i,j).gt._ZERO_) then
                  do it=1,z0d_iters
                     z0d = zub0(i,j) + _TENTH_*avmmol/(cdsqrt*velU(i,j))
!                    KK-TODO: clipping of z0d at DU as in the old code?
                     cdsqrt = kappa / log( _ONE_ + _HALF_*DU1(i,j)/z0d )
                  end do
               end if
               cd = max( cd_min , cdsqrt**2) ! see Blumberg and Mellor (1987)
               ru(i,j) = cd * velU(i,j)
               p_zub(i,j) = z0d
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO


!     The y-direction

!     calculate u_velC
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin-HALO,jmax+HALO
#endif
         do i=imin-HALO+1,imax+HALO-1
            if (az(i,j) .ge. 1) then
               work2d(i,j) = _HALF_ * ( u_vel(i-1,j) + u_vel(i,j) )
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO

#ifdef SLICE_MODEL
!$OMP SINGLE
      work2d(imin-HALO+1:imax+HALO-1,j+1) = work2d(imin-HALO+1:imax+HALO-1,j)
!$OMP END SINGLE
#endif

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin-HALO,jmax+HALO-1
#endif
         do i=imin-HALO+1,imax+HALO-1
            if ( av(i,j) .ge. 1 ) then
               velV(i,j) = sqrt( (_HALF_*(work2d(i,j)+work2d(i,j+1)))**2 + v_vel(i,j)**2 )
               z0d = zvb0(i,j)
!              Note (KK): note shifting of log profile so that V(-H)=0
               cdsqrt = kappa / log( _ONE_ + _HALF_*DV1(i,j)/z0d )
               if (avmmol.gt._ZERO_ .and. velV(i,j).gt._ZERO_) then
                  do it=1,z0d_iters
                     z0d = zvb0(i,j) + _TENTH_*avmmol/(cdsqrt*velV(i,j))
!                    KK-TODO: clipping of z0d at DV as in the old code?
                     cdsqrt = kappa / log( _ONE_ + _HALF_*DV1(i,j)/z0d )
                  end do
               end if
               cd = max( cd_min , cdsqrt**2) ! see Blumberg and Mellor (1987)
               rv(i,j) = cd * velV(i,j)
               p_zvb(i,j) = z0d
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO

!$OMP END PARALLEL

#ifdef SLICE_MODEL
      ru   (imin-HALO  :imax+HALO-1,j+1) = ru   (imin-HALO  :imax+HALO-1,j)
      rv   (imin-HALO+1:imax+HALO-1,j-1) = rv   (imin-HALO+1:imax+HALO-1,j)
      rv   (imin-HALO+1:imax+HALO-1,j+1) = rv   (imin-HALO+1:imax+HALO-1,j)
      p_zub(imin-HALO  :imax+HALO-1,j+1) = p_zub(imin-HALO  :imax+HALO-1,j)
      p_zvb(imin-HALO+1:imax+HALO-1,j-1) = p_zvb(imin-HALO+1:imax+HALO-1,j)
      p_zvb(imin-HALO+1:imax+HALO-1,j+1) = p_zvb(imin-HALO+1:imax+HALO-1,j)
#endif

      if (waves_method .ne. NO_WAVES) then
         call toc(TIM_BOTTFRIC)
         if (present(taubmax)) then
            call bottom_friction_waves(U1,V1,DU1,DV1,Dvel,velU,velV,ru,rv,p_zub,p_zvb,taubmax)
         else
            call bottom_friction_waves(U1,V1,DU1,DV1,Dvel,velU,velV,ru,rv,p_zub,p_zvb)
         end if
         call tic(TIM_BOTTFRIC)
      end if

   end if

   CALL toc(TIM_BOTTFRIC)
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
