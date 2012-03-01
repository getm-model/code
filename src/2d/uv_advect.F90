#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: uv_advect - 2D advection of momentum \label{sec-uv-advect}
!
! !INTERFACE:
   subroutine uv_advect(U,V,DU,DV)
!
! !DESCRIPTION:
!
! The advective terms in the vertically integrated
! momentum equation are discretised in
! a momentum-conservative form. This is carried out here for the
! advective terms in the $U$-equation (\ref{UMOM}) and the
! $V$-equation (\ref{VMOM}) (after applying the curvilinear
! coordinate transformationand multiplying these
! equations with $mn$).
!
! First advection term in (\ref{UMOM}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! \left(mn\,\partial_{\cal X}\left(\frac{U^2}{Dn}\right)\right)_{i,j}\approx \\ \\
! \quad
! \displaystyle
! \frac{
! \frac12(U_{i+1,j}+U_{i,j})\tilde u_{i+1,j}\Delta y^c_{i+1,j}-
! \frac12(U_{i,j}+U_{i-1,j})\tilde u_{i,j}\Delta y^c_{i,j}
! }{\Delta x^u_{i,j}\Delta y^u_{i,j}}
! \end{array}
! \end{equation}
!
! For the upwind scheme used here, the inter-facial velocities which are defined
! on T-points are here
! calculated as:
!
! \begin{equation}
! \tilde u_{i,j}=
! \left\{
! \begin{array}{ll}
! \displaystyle
! \frac{U_{i-1,j}}{D^u_{i-1,j}} & \mbox{ for } \frac12(U_{i,j}+U_{i-1,j})>0\\ \\
! \displaystyle
! \frac{U_{i,j}}{D^u_{i,j}} & \mbox{ else. }
! \end{array}
! \right.
! \end{equation}
!
! Second advection term in (\ref{UMOM}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! \left(mn\,\partial_{\cal Y}y\left(\frac{UV}{Dm}\right)\right)_{i,j,k}\approx \\ \\
! \displaystyle
! \quad
! \frac{
! \frac12(V_{i+1,j}+V_{i,j})\tilde u_{i,j}\Delta x^+_{i,j}-
! \frac12(V_{i+1,j-1}+V_{i,j-1})\tilde u_{i,j-1}\Delta x^+_{i,j-1}
! }{\Delta x^u_{i,j}\Delta y^u_{i,j}}
! \end{array}
! \end{equation}
!
! For the upwind scheme used here, the inter-facial
! velocities which are defined on
! X-points are here
! calculated as:
!
! \begin{equation}
! \tilde u_{i,j}=
! \left\{
! \begin{array}{ll}
! \displaystyle
! \frac{U_{i,j}}{D^u_{i,j}} & \mbox{ for } \frac12(V_{i+1,j,k}+V_{i,j,k})>0\\ \\
! \displaystyle
! \frac{U_{i,j+1}}{D^u_{i,j+1}} & \mbox{ else. }
! \end{array}
! \right.
! \end{equation}
!
! First advection term in (\ref{VMOM}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! \left(mn\,\partial_{\cal X}\left(\frac{UV}{Dn}\right)\right)_{i,j,k}\approx \\ \\
! \displaystyle
! \quad
! \frac{
! \frac12(U_{i,j+1}+U_{i,j})\tilde v_{i,j}\Delta y^+_{i,j}-
! \frac12(U_{i-1,j+1}+U_{i-1,j})\tilde v_{i-1,j}\Delta y^+_{i-1,j}
! }{\Delta x^v_{i,j}\Delta y^v_{i,j}}
! \end{array}
! \end{equation}
!
! For the upwind scheme used here, the interfacial
! velocities which are defined on
! X-points are here
! calculated as:
!
! \begin{equation}
! \tilde v_{i,j}=
! \left\{
! \begin{array}{ll}
! \displaystyle
! \frac{V_{i,j}}{D^v_{i,j}} & \mbox{ for } \frac12(U_{i+1,j}+U_{i,j})>0\\ \\
! \displaystyle
! \frac{V_{i+1,j}}{D^v_{i+1,j}} & \mbox{ else. }
! \end{array}
! \right.
! \end{equation}
!
! Second advection term in (\ref{VMOM}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! \left(mn\,\partial_{\cal Y}\left(\frac{V^2}{Dm}\right)\right)_{i,j,k}\approx \\ \\
! \quad
! \displaystyle
! \frac{
! \frac12(V_{i,j+1}+V_{i,j})\tilde v_{i,j+1}\Delta x^c_{i,j+1}-
! \frac12(V_{i,j}+V_{i,j-1})\tilde v_{i,j}\Delta x^c_{i,j}
! }{\Delta x^v_{i,j}\Delta y^v_{i,j}}
! \end{array}
! \end{equation}
!
! For the upwind scheme used here, the interfacial velocities which are defined
! on T-points are here
! calculated as:
!
! \begin{equation}
! \tilde v_{i,j}=
! \left\{
! \begin{array}{ll}
! \displaystyle
! \frac{V_{i,j-1}}{D^v_{i,j-1}} & \mbox{ for } \frac12(V_{i,j}+V_{i,j-1})>0\\ \\
! \displaystyle
! \frac{V_{i,j}}{D^v_{i,j}} & \mbox{ else. }
! \end{array}
! \right.
! \end{equation}
!
! When working with the option {\tt SLICE\_MODEL}, the calculation of
! all gradients in $y$-direction is suppressed.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,az,au,av,ax
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxv,dyu
#else
   use domain, only: dx,dy
#endif
   use m2d, only: dtm,vel_adv_split2d,vel_adv_scheme
   use variables_2d, only: UEx,VEx,fadv,Uadv,Vadv,Dadv,DUadv,DVadv
   use advection, only: UPSTREAM,J7,do_advection
   use halo_zones, only: update_2d_halo,wait_halo,U_TAG,V_TAG
   use getm_timers, only: tic,toc,TIM_UVADV,TIM_UVADVH
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(in)        :: U,V
   REALTYPE,dimension(E2DFIELD),target,intent(in) :: DU,DV
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer :: i,j
   REALTYPE,dimension(:,:),pointer :: pDadv
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'uv_advect() # ',Ncall
#endif
#ifdef SLICE_MODEL
   j = jmax/2 ! this MUST NOT be changed!!!
#endif
   call tic(TIM_UVADV)

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j)                                         &
!$OMP          PRIVATE(i)


!  Here begins dimensional split advection for u-velocity

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
   do j=jmin-HALO,jmax+HALO
#endif
      do i=imin-HALO,imax+HALO-1
!        the velocity to be transported
         fadv(i,j) = U(i,j)/DU(i,j)
         if (vel_adv_scheme .ne. J7) then
!           Note (KK): Uadv defined on T-points (future U-points)
!                      Vadv defined on X-points (future V-points)
            Uadv(i,j) = _HALF_*( U(i,j) + U(i+1,j) )
            Vadv(i,j) = _HALF_*( V(i,j) + V(i+1,j) )
         end if
!        Note (KK): DU only valid until imax+1
!                   therefore DUadv only valid until imax
         DUadv(i,j) = _HALF_*( DU(i,j) + DU(i+1,j) )
!        Note (KK): DV only valid until jmax+1
!                   therefore DVadv only valid until jmax+1
         DVadv(i,j) = _HALF_*( DV(i,j) + DV(i+1,j) )
      end do
#ifndef SLICE_MODEL
   end do
#endif
!$OMP END DO NOWAIT

   if (vel_adv_scheme .eq. J7) then
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin-HALO,jmax+HALO
#endif
         do i=imin-HALO,imax+HALO-1
!           Note (KK): [U|V]adv defined on T-points (future U-points)
!                      Dadv defined on X-points (future V-points)
!                      note that Dadv is shifted to j+1 !!!
            Uadv(i,j) = _HALF_*( U(i,j)*DYU + U(i+1,j)*DYUIP1 )
            if (j .ne. jmin-HALO) then
               if (az(i+1,j) .eq. 1) then
                  Vadv(i,j) = _HALF_*( V(i+1,j-1)*DXVPM + V(i+1,j)*DXVIP1 )
               else if(az(i+1,j) .eq. 2) then
!                 Note (KK): can be included into case above, when
!                            V is properly mirrored across n/s open bdys
                  if (av(i+1,j) .eq. 2) then ! southern open bdy
                     Vadv(i,j) = V(i+1,j)*_HALF_*( DXVPM + DXVIP1 )
                  else if (av(i+1,j-1) .eq. 2) then ! northern open bdy
                     Vadv(i,j) = V(i+1,j-1)*_HALF_*( DXVPM + DXVIP1 )
                  end if
               end if
               if (ax(i,j) .eq. 0) then
                  if (au(i,j-1) .ne. 0) then
                     Dadv(i,j) = DU(i,j-1)
                  else if (au(i,j) .ne. 0) then
                     Dadv(i,j) = DU(i,j)
                  end if
               else
                  Dadv(i,j) = _HALF_*( DU(i,j-1) + DU(i,j) )
               end if
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO

#ifdef SLICE_MODEL
      Dadv(imin-HALO:imax+HALO-1,j+1) = Dadv(imin-HALO:imax+HALO-1,j)
      Uadv(imin-HALO:imax+HALO-1,j+1) = Uadv(imin-HALO:imax+HALO-1,j)
      Vadv(imin-HALO:imax+HALO-1,j+1) = Vadv(imin-HALO:imax+HALO-1,j)
#endif

#ifndef SLICE_MODEL
!     OMP-NOTE (KK): j loop must not be changed and cannot be threaded!
      do j=jmin-HALO,jmax+HALO-1
#endif
!$OMP DO SCHEDULE(RUNTIME)
         do i=imin-HALO,imax+HALO-1
!           Note (KK): [U|V]adv defined on V-points (future X-points)
!                      no change of Uadv for northern closed bdy
!                      Dadv defined on U-points (future T-points)
!                      note that Dadv is not shifted anymore !
!           KK-TODO: Vadv for western/eastern open bdys ?
            if (az(i+1,j) .eq. 0) then ! southern closed bdy
               Uadv(i,j) = Uadv(i,j+1)
            else if (az(i+1,j+1) .ne. 0) then
               Uadv(i,j) = _HALF_*( Uadv(i,j) + Uadv(i,j+1) )
               Vadv(i,j) = _HALF_*( Vadv(i,j) + Vadv(i,j+1) )
            end if
            Dadv(i,j) = _HALF_*( Dadv(i,j) + Dadv(i,j+1) )
         end do
!$OMP END DO
#ifndef SLICE_MODEL
      end do
#endif

!$OMP SINGLE
      pDadv => Dadv
!$OMP END SINGLE

   else

!$OMP SINGLE
      pDadv => DU
!$OMP END SINGLE

   end if

!$OMP SINGLE
   if (vel_adv_scheme.ne.UPSTREAM .and. vel_adv_scheme.ne.J7) then
!     we need to update fadv(imax+HALO,jmin-HALO:jmax+HALO)
      call tic(TIM_UVADVH)
      call update_2d_halo(fadv,fadv,au,imin,jmin,imax,jmax,U_TAG)
      call wait_halo(U_TAG)
      call toc(TIM_UVADVH)
   end if

   call do_advection(dtm,fadv,Uadv,Vadv,DUadv,DVadv,pDadv,pDadv,  &
                     vel_adv_scheme,vel_adv_split2d,_ZERO_,U_TAG, &
                     advres=UEx)
!$OMP END SINGLE


!  Here begins dimensional split advection for v-velocity

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
   do j=jmin-HALO,jmax+HALO-1
#endif
      do i=imin-HALO,imax+HALO
!        the velocity to be transported
         fadv(i,j) = V(i,j)/DV(i,j)
         if (vel_adv_scheme .ne. J7) then
!           Note (KK): Uadv defined on X-points (future U-points)
!                      Vadv defined on T-points (future V-points)
            Uadv(i,j) = _HALF_*( U(i,j) + U(i,j+1) )
            Vadv(i,j) = _HALF_*( V(i,j) + V(i,j+1) )
         end if
!        Note (KK): DU only valid until imax+1
!                   therefore DUadv only valid until imax+1
         DUadv(i,j) = _HALF_*( DU(i,j) + DU(i,j+1) )
!        Note (KK): DV only valid until jmax+1
!                   therefore DVadv only valid until jmax
         DVadv(i,j) = _HALF_*( DV(i,j) + DV(i,j+1) )
      end do
#ifndef SLICE_MODEL
   end do
#endif
!$OMP END DO NOWAIT

   if (vel_adv_scheme .eq. J7) then
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin-HALO,jmax+HALO-1
#endif
         do i=imin-HALO,imax+HALO
!           Note (KK): [U|V]adv defined on T-points (future V-points)
!                      Dadv defined on X-points (future U-points)
!                      note that Dadv is shifted to i+1 !!!
            if (i .ne. imin-HALO) then
               if (az(i,j+1) .eq. 1) then
                  Uadv(i,j) = _HALF_*( U(i-1,j+1)*DYUMP + U(i,j+1)*DYUJP1 )
               else if(az(i,j+1) .eq. 2) then
!                 Note (KK): can be included into case above, when
!                            U is properly mirrored across w/e open bdys
                  if (au(i,j+1) .eq. 2) then ! western open bdy
                     Uadv(i,j) = U(i,j+1)*_HALF_*( DYUMP + DYUJP1 )
                  else if (au(i-1,j+1) .eq. 2) then ! eastern open bdy
                     Uadv(i,j) = U(i-1,j+1)*_HALF_*( DYUMP + DYUJP1 )
                  end if
               end if
               if (ax(i,j) .eq. 0) then
                  if (av(i-1,j) .ne. 0) then
                     Dadv(i,j) = DV(i-1,j)
                  else if (av(i,j) .ne. 0) then
                     Dadv(i,j) = DV(i,j)
                  end if
               else
                  Dadv(i,j) = _HALF_*( DV(i-1,j) + DV(i,j) )
               end if
            end if
            Vadv(i,j) = _HALF_*( V(i,j)*DXV + V(i,j+1)*DXVJP1 )
         end do
#ifdef SLICE_MODEL
!$OMP END DO
#endif

!        OMP-NOTE (KK): i loop must not be changed and cannot be threaded!
         do i=imin-HALO,imax+HALO-1
!           Note (KK): [U|V]adv defined on U-points (future X-points)
!                      no change of Vadv for eastern closed bdy
!                      Dadv defined on V-points (future T-points)
!                      note that Dadv is not shifted anymore !
!           KK-TODO: Uadv for northern/southern open bdys ?
            if (az(i,j+1) .eq. 0) then ! western closed bdy
               Vadv(i,j) = Vadv(i+1,j)
            else if (az(i+1,j+1) .ne. 0) then
               Uadv(i,j) = _HALF_*( Uadv(i,j) + Uadv(i+1,j) )
               Vadv(i,j) = _HALF_*( Vadv(i,j) + Vadv(i+1,j) )
            end if
            Dadv(i,j) = _HALF_*( Dadv(i,j) + Dadv(i+1,j) )
         end do
#ifndef SLICE_MODEL
      end do
!$OMP END DO
#endif

!$OMP SINGLE
      pDadv => Dadv
!$OMP END SINGLE

   else

!$OMP SINGLE
      pDadv => DV
!$OMP END SINGLE

   end if

!$OMP END PARALLEL

   if (vel_adv_scheme.ne.UPSTREAM .and. vel_adv_scheme.ne.J7) then
!     we need to update fadv(imin-HALO:imax+HALO,jmax+HALO)
      call tic(TIM_UVADVH)
      call update_2d_halo(fadv,fadv,av,imin,jmin,imax,jmax,V_TAG)
      call wait_halo(V_TAG)
      call toc(TIM_UVADVH)
   end if

   call do_advection(dtm,fadv,Uadv,Vadv,DUadv,DVadv,pDadv,pDadv,  &
                     vel_adv_scheme,vel_adv_split2d,_ZERO_,V_TAG, &
                     advres=VEx)

   call toc(TIM_UVADV)
#ifdef DEBUG
   write(debug,*) 'Leaving uv_advect()'
   write(debug,*)
#endif
   return
   end subroutine uv_advect
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
