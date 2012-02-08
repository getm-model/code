#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: bottom_friction - calculates the 2D bottom friction.
!
! !INTERFACE:
   subroutine bottom_friction(U,V,DU,DV,ru,rv,zub,zvb)
!  Note (KK): keep in sync with interface in m2d_general.F90
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
   use parameters,only: kappa,avmmol
   use domain,only: imin,imax,jmin,jmax,az,au,av,zub0,zvb0
   use variables_2d,only: PP
   use getm_timers,only: tic,toc,TIM_BOTTFRIC
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(in)           :: U,V,DU,DV
!
! !OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(out)          :: ru,rv
   REALTYPE,dimension(E2DFIELD),intent(out),optional :: zub,zvb
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  !LOCAL VARIABLES:
   REALTYPE,dimension(:,:),allocatable,save :: u_vel,v_vel
   REALTYPE                                 :: vel,cd,z0b
   integer                                  :: i,j,it,rc
   logical,save                             :: first=.true.
   integer,parameter                        :: it_max=0
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'bottom_friction() # ',Ncall
#endif
#ifdef SLICE_MODEL
   j = jmax/2
#endif
   CALL tic(TIM_BOTTFRIC)

   if (first) then
      allocate(u_vel(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (u_vel)'
      u_vel=_ZERO_

      allocate(v_vel(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (v_vel)'
      v_vel=_ZERO_

      first = .false.
   end if

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,vel,cd,it,z0b)

!  zonal velocity
#ifndef SLICE_MODEL
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO
#endif
      do i=imin-HALO,imax+HALO-1
         if (au(i,j) .ge. 1) then
            u_vel(i,j) = U(i,j)/DU(i,j)
         end if
      end do
#ifndef SLICE_MODEL
   end do
!$OMP END DO NOWAIT
#endif

!  meridional velocity
#ifndef SLICE_MODEL
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO-1
#endif
      do i=imin-HALO,imax+HALO
         if (av(i,j) .ge. 1) then
            v_vel(i,j) = V(i,j)/DV(i,j)
         end if
      end do
#ifndef SLICE_MODEL
   end do
!$OMP END DO
#else
   v_vel(:,j-1) = v_vel(:,j)
#endif

!  The x-direction

#ifndef SLICE_MODEL
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO+1,jmax+HALO-1
#endif
!     calculate v_velC
      do i=imin-HALO,imax+HALO
         if (az(i,j) .ge. 1) then
            PP(i,j) = _HALF_ * ( v_vel(i,j-1) + v_vel(i,j) )
         end if
      end do

      do i=imin-HALO,imax+HALO-1
         if (au(i,j).eq.1 .or. au(i,j).eq.2) then
            vel = sqrt( u_vel(i,j)**2 + (_HALF_*(PP(i,j)+PP(i+1,j)))**2 )
            z0b = zub0(i,j)
!           Note (KK): note shifting of log profile so that U(-H)=0
            cd = kappa / log( _ONE_ + _HALF_*DU(i,j)/z0b )
            if (avmmol.gt._ZERO_ .and. vel.gt._ZERO_) then
               do it=1,it_max
                  z0b = zub0(i,j) + _TENTH_*avmmol/(cd*vel)
                  cd = kappa / log( _ONE_ + _HALF_*DU(i,j)/z0b )
               end do
            end if
            !cd = max( 0.0025d0 , cd) ! see Blumberg and Mellor (1987)
            ru(i,j) = (cd**2) * vel
            if (present(zub)) then
               zub(i,j) = z0b
            end if
         end if
      end do
#ifndef SLICE_MODEL
   end do
!$OMP END DO
#else
   ru(imin-HALO:imax+HALO-1,j+1) = ru(imin-HALO:imax+HALO-1,j)
   if (present(zub)) then
      zub(imin-HALO:imax+HALO-1,j+1) = zub(imin-HALO:imax+HALO-1,j)
   end if
#endif

!  The y-direction

!  calculate u_velC
#ifndef SLICE_MODEL
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO
#endif
      do i=imin-HALO+1,imax+HALO-1
         if (az(i,j) .ge. 1) then
            PP(i,j) = _HALF_ * ( u_vel(i-1,j) + u_vel(i,j) )
         end if
      end do
#ifndef SLICE_MODEL
!$OMP END DO
   end do
#else
   PP(imin-HALO+1:imax+HALO-1,j+1) = PP(imin-HALO+1:imax+HALO-1,j)
#endif

#ifndef SLICE_MODEL
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO-1
#endif
      do i=imin-HALO+1,imax+HALO-1
         if (av(i,j).eq.1 .or. av(i,j).eq.2) then
            vel = sqrt( (_HALF_*(PP(i,j)+PP(i,j+1)))**2 + v_vel(i,j)**2 )
            z0b = zvb0(i,j)
!           Note (KK): note shifting of log profile so that V(-H)=0
            cd = kappa / log( _ONE_ + _HALF_*DV(i,j)/z0b )
            if (avmmol.gt._ZERO_ .and. vel.gt._ZERO_) then
               do it=1,it_max
                  z0b = zvb0(i,j) + _TENTH_*avmmol/(cd*vel)
                  cd = kappa / log( _ONE_ + _HALF_*DV(i,j)/z0b )
               end do
            end if
            !cd = max( 0.0025d0 , cd) ! see Blumberg and Mellor (1987)
            rv(i,j) = (cd**2) * vel
            if (present(zvb)) then
               zvb(i,j) = z0b
            end if
         end if
      end do
#ifndef SLICE_MODEL
   end do
!$OMP END DO
#else
   rv(imin-HALO+1:imax+HALO-1,j-1) = rv(imin-HALO+1:imax+HALO-1,j)
   rv(imin-HALO+1:imax+HALO-1,j+1) = rv(imin-HALO+1:imax+HALO-1,j)
   if (present(zvb)) then
      zvb(imin-HALO+1:imax+HALO-1,j-1) = zvb(imin-HALO+1:imax+HALO-1,j)
      zvb(imin-HALO+1:imax+HALO-1,j+1) = zvb(imin-HALO+1:imax+HALO-1,j)
   end if
#endif

!$OMP END PARALLEL

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
