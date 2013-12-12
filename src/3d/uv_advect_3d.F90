#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: uv_advect_3d - 3D momentum advection \label{sec-uv-advect-3d}
!
! !INTERFACE:
   subroutine uv_advect_3d()
!
! !DESCRIPTION:
!
! Wrapper to prepare and do calls to routine {\tt do\_advection\_3d}
! (see section \ref{sec-do-advection-3d} on page
! \pageref{sec-do-advection-3d}) to calculate the advection terms of the
! 3D velocities.
!
! If {\tt save\_numerical\_analyses} is set to {\tt .true.}, the
! numerical dissipation is calculated using the method suggested
! by \cite{BURCHARD12}.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax,az,au,av,ax
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxv,dyu
#else
   use domain, only: dx,dy
#endif
   use m3d, only: vel3d_adv_split,vel3d_adv_hor,vel3d_adv_ver
   use variables_3d, only: dt,uu,vv,ww,hn,hvel,hun,hvn,uuEx,vvEx
   use advection, only: NOADV,UPSTREAM,J7
   use advection_3d, only: do_advection_3d
   use halo_zones, only: update_3d_halo,wait_halo,U_TAG,V_TAG
   use variables_3d, only: do_numerical_analyses
   use variables_3d, only: numdis3d,numdis2d
#ifdef _MOMENTUM_TERMS_
   use domain, only: dry_u,dry_v
   use variables_3d, only: adv_u,adv_v
#endif
   use getm_timers, only: tic,toc,TIM_UVADV3D,TIM_UVADV3DH
!$ use omp_lib
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                             :: i,j,k
   REALTYPE,dimension(I3DFIELD)        :: fadv3d,uuadv,vvadv,wwadv
   REALTYPE,dimension(I3DFIELD),target :: hnadv,huadv,hvadv
   REALTYPE,dimension(:,:,:),pointer   :: phadv,phuadv,phvadv
   REALTYPE,dimension(I3DFIELD)        :: work3d,hires
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'uv_advect_3d() # ',Ncall
#endif
#ifdef SLICE_MODEL
   j = jmax/2 ! this MUST NOT be changed!!!
#endif
   call tic(TIM_UVADV3D)

   if (vel3d_adv_hor.ne.NOADV .or. vel3d_adv_ver.ne.NOADV) then

!     Note (KK): wwadv(kmax) will be overwritten by ww(kmax) anyway
      wwadv(:,:,0) = _ZERO_

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j)                                         &
!$OMP          PRIVATE(i,k)


!     Here begins dimensional split advection for u-velocity

!$OMP SINGLE
!     KK-TODO: _POINTER_REMAP_, but this requires that h[u|v]n become
!              pointers like mask_[u|v]flux in do_advection_3d and similar
!              for D[U|V] in do_advection and in uv_advect
!#ifdef _POINTER_REMAP_
!      phuadv(imin-HALO:,jmin-HALO:,0:) => hvel(1+_IRANGE_HALO_,_JRANGE_HALO_,_KRANGE_)
!#else
      phuadv => huadv
      phuadv(_IRANGE_HALO_-1,:,:) = hvel(1+_IRANGE_HALO_,:,:)
!#endif
      phvadv => hvadv
!$OMP END SINGLE

      do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
         do j=jmin-HALO,jmax+HALO
#endif
            do i=imin-HALO,imax+HALO-1
!              the velocity to be transported
               fadv3d(i,j,k) = uu(i,j,k)/hun(i,j,k)
               if (vel3d_adv_hor .ne. J7) then
                  uuadv(i,j,k) = _HALF_*( uu(i,j,k) + uu(i+1,j,k) )
                  vvadv(i,j,k) = _HALF_*( vv(i,j,k) + vv(i+1,j,k) )
               end if
               wwadv(i,j,k) = _HALF_*( ww(i,j,k) + ww(i+1,j,k) )
!              Note (KK): hvn only valid until jmax+1
!                         therefore phvadv only valid until jmax+1
               phvadv(i,j,k) = _HALF_*( hvn(i,j,k) + hvn(i+1,j,k) )
            end do
#ifndef SLICE_MODEL
         end do
#endif
!$OMP END DO NOWAIT
      end do

      if (vel3d_adv_hor .eq. J7) then
         do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
            do j=jmin-HALO,jmax+HALO
#endif
               do i=imin-HALO,imax+HALO-1
!                 Note (KK): [uu|vv]adv defined on T-points (future U-points)
!                            hnadv defined on X-points (future V-points)
!                            note that hnadv is shifted to j+1 !!!
                  uuadv(i,j,k) = _HALF_*( uu(i,j,k)*DYU + uu(i+1,j,k)*DYUIP1 )
                  if (j .ne. jmin-HALO) then
                     if (az(i+1,j) .eq. 1) then
                        vvadv(i,j,k) = _HALF_*( vv(i+1,j-1,k)*DXVPM + vv(i+1,j,k)*DXVIP1 )
                     else if(az(i+1,j) .eq. 2) then
!                       Note (KK): can be included into case above, when
!                                  vv is properly mirrored across n/s open bdys
                        if (av(i+1,j) .eq. 2) then ! southern open bdy
                           vvadv(i,j,k) = vv(i+1,j,k)*_HALF_*( DXVPM + DXVIP1 )
                        else if (av(i+1,j-1) .eq. 2) then ! northern open bdy
                           vvadv(i,j,k) = vv(i+1,j-1,k)*_HALF_*( DXVPM + DXVIP1 )
                        end if
                     end if
                     if (ax(i,j) .eq. 0) then
                        if (au(i,j-1) .ne. 0) then
                           hnadv(i,j,k) = hun(i,j-1,k)
                        else if (au(i,j) .ne. 0) then
                           hnadv(i,j,k) = hun(i,j,k)
                        end if
                     else
                        hnadv(i,j,k) = _HALF_*( hun(i,j-1,k) + hun(i,j,k) )
                     end if
                  end if
               end do
#ifndef SLICE_MODEL
            end do
#endif
!$OMP END DO
         end do

#ifdef SLICE_MODEL
         hnadv(imin-HALO:imax+HALO-1,j+1,:) = hnadv(imin-HALO:imax+HALO-1,j,:)
         uuadv(imin-HALO:imax+HALO-1,j+1,:) = uuadv(imin-HALO:imax+HALO-1,j,:)
         vvadv(imin-HALO:imax+HALO-1,j+1,:) = vvadv(imin-HALO:imax+HALO-1,j,:)
#endif

         do k=1,kmax
#ifndef SLICE_MODEL
!           OMP-NOTE (KK): j loop must not be changed and cannot be threaded!
            do j=jmin-HALO,jmax+HALO-1
#endif
!$OMP DO SCHEDULE(RUNTIME)
               do i=imin-HALO,imax+HALO-1
!                 Note (KK): [uu|vv]adv defined on V-points (future X-points)
!                            no change of uuadv for northern closed bdy
!                            hnadv defined on U-points (future T-points)
!                            note that hnadv is not shifted anymore !
!                 KK-TODO: hadv for western/eastern open bdys ?
                  if (az(i+1,j) .eq. 0) then ! southern closed bdy
                     uuadv(i,j,k) = uuadv(i,j+1,k)
                  else if (az(i+1,j+1) .ne. 0) then
                     uuadv(i,j,k) = _HALF_*( uuadv(i,j,k) + uuadv(i,j+1,k) )
                     vvadv(i,j,k) = _HALF_*( vvadv(i,j,k) + vvadv(i,j+1,k) )
                  end if
                  hnadv(i,j,k) = _HALF_*( hnadv(i,j,k) + hnadv(i,j+1,k) )
               end do
!$OMP END DO
#ifndef SLICE_MODEL
            end do
#endif
         end do

!$OMP SINGLE
         phadv => hnadv
!$OMP END SINGLE

      else

!$OMP SINGLE
         phadv => hun
!$OMP END SINGLE

      end if

!$OMP SINGLE
      if (vel3d_adv_hor.ne.NOADV .and. vel3d_adv_hor.ne.UPSTREAM .and. vel3d_adv_hor.ne.J7) then
!        we need to update fadv3d(imax+HALO,jmin-HALO:jmax+HALO)
         call tic(TIM_UVADV3DH)
         call update_3d_halo(fadv3d,fadv3d,au,imin,jmin,imax,jmax,kmax,U_TAG)
         call wait_halo(U_TAG)
         call toc(TIM_UVADV3DH)
      end if
!$OMP END SINGLE

      if (do_numerical_analyses) then
         do k=1,kmax ! calculate square of u-velocity before advection step
!$OMP DO SCHEDULE(RUNTIME)
            do j=jmin-HALO,jmax+HALO
               do i=imin-HALO,imax+HALO
                  work3d(i,j,k) = fadv3d(i,j,k)**2
               end do
            end do
!$OMP END DO NOWAIT
         end do
!$OMP BARRIER
      end if

!$OMP SINGLE
      call do_advection_3d(dt,fadv3d,uuadv,vvadv,wwadv,phuadv,phvadv,phadv,phadv,    &
                           vel3d_adv_split,vel3d_adv_hor,vel3d_adv_ver,_ZERO_,U_TAG, &
                           advres=uuEx)
!$OMP END SINGLE

#ifdef _MOMENTUM_TERMS_
      do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
         do j=jmin,jmax
            do i=imin,imax
               adv_u(i,j,k) = dry_u(i,j) * uuEx(i,j,k)
            end do
         end do
!$OMP END DO NOWAIT
      end do
#endif

      if (do_numerical_analyses) then

!$OMP SINGLE
         call do_advection_3d(dt,work3d,uuadv,vvadv,wwadv,phuadv,phvadv,phadv,phadv,    &
                              vel3d_adv_split,vel3d_adv_hor,vel3d_adv_ver,_ZERO_,U_TAG, &
                              hires=hires)

         numdis2d = _ZERO_
!$OMP END SINGLE

         do k=1,kmax ! calculate kinetic energy dissipaion rate for u-velocity
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
            do j=jmin,jmax
#endif
               do i=imin,imax
                  work3d(i,j,k) = ( work3d(i,j,k) - fadv3d(i,j,k)**2 ) / dt
               end do
#ifndef SLICE_MODEL
            end do
#endif
!$OMP END DO NOWAIT
         end do

!$OMP BARRIER
!$OMP SINGLE
         call update_3d_halo(work3d,work3d,au,imin,jmin,imax,jmax,kmax,U_TAG)
         call wait_halo(U_TAG)
!$OMP END SINGLE

         do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
            do j=jmin,jmax
#endif
               do i=imin,imax
                  if (az(i,j) .eq. 1) then
                     numdis3d(i,j,k) = _HALF_*( work3d(i-1,j,k) + work3d(i,j,k) )
                     numdis2d(i,j) = numdis2d(i,j)                             &
                                    +_HALF_*( work3d(i-1,j,k)*hires(i-1,j,k)   &
                                             +work3d(i  ,j,k)*hires(i  ,j,k) )
                  end if
               end do
#ifndef SLICE_MODEL
            end do
#endif
!$OMP END DO
         end do

      end if


!     Here begins dimensional split advection for v-velocity

!$OMP SINGLE
!     KK-TODO: _POINTER_REMAP_, but this requires that h[u|v]n become
!              pointers like mask_[u|v]flux in do_advection_3d and similar
!              for D[U|V] in do_advection and in uv_advect
      phuadv => huadv
!#ifdef _POINTER_REMAP_
!      phvadv(imin-HALO:,jmin-HALO:,0:) => hvel(_IRANGE_HALO_,1+_JRANGE_HALO_,_KRANGE_)
!#else
      phvadv => hvadv
      phvadv(:,_JRANGE_HALO_-1,:) = hvel(:,1+_JRANGE_HALO_,:)
!#endif
!$OMP END SINGLE


      do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
         do j=jmin-HALO,jmax+HALO-1
#endif
            do i=imin-HALO,imax+HALO
!              the velocity to be transported
               fadv3d(i,j,k) = vv(i,j,k)/hvn(i,j,k)
               if (vel3d_adv_hor .ne. J7) then
                  uuadv(i,j,k) = _HALF_*( uu(i,j,k) + uu(i,j+1,k) )
                  vvadv(i,j,k) = _HALF_*( vv(i,j,k) + vv(i,j+1,k) )
               end if
               wwadv(i,j,k) = _HALF_*( ww(i,j,k) + ww(i,j+1,k) )
!              Note (KK): hun only valid until imax+1
!                         therefore phuadv only valid until imax+1
               phuadv(i,j,k) = _HALF_*( hun(i,j,k) + hun(i,j+1,k) )
            end do
#ifndef SLICE_MODEL
         end do
#endif
!$OMP END DO NOWAIT
      end do

      if (vel3d_adv_hor .eq. J7) then
         do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
            do j=jmin-HALO,jmax+HALO-1
#endif
               do i=imin-HALO,imax+HALO
!                 Note (KK): [uu|vv]adv defined on T-points (future V-points)
!                            hnadv defined on X-points (future U-points)
!                            note that hnadv is shifted to i+1 !!!
                  if (i .ne. imin-HALO) then
                     if (az(i,j+1) .eq. 1) then
                        uuadv(i,j,k) = _HALF_*( uu(i-1,j+1,k)*DYUMP + uu(i,j+1,k)*DYUJP1 )
                     else if(az(i,j+1) .eq. 2) then
!                       Note (KK): can be included into case above, when
!                                  uu is properly mirrored across w/e open bdys
                        if (au(i,j+1) .eq. 2) then ! western open bdy
                           uuadv(i,j,k) = uu(i,j+1,k)*_HALF_*( DYUMP + DYUJP1 )
                        else if (au(i-1,j+1) .eq. 2) then ! eastern open bdy
                           uuadv(i,j,k) = uu(i-1,j+1,k)*_HALF_*( DYUMP + DYUJP1 )
                        end if
                     end if
                     if (ax(i,j) .eq. 0) then
                        if (av(i-1,j) .ne. 0) then
                           hnadv(i,j,k) = hvn(i-1,j,k)
                        else if (av(i,j) .ne. 0) then
                           hnadv(i,j,k) = hvn(i,j,k)
                        end if
                     else
                        hnadv(i,j,k) = _HALF_*( hvn(i-1,j,k) + hvn(i,j,k) )
                     end if
                  end if
                  vvadv(i,j,k) = _HALF_*( vv(i,j,k)*DXV + vv(i,j+1,k)*DXVJP1 )
               end do
#ifdef SLICE_MODEL
!$OMP END DO
#endif

!              OMP-NOTE (KK): i loop must not be changed and cannot be threaded!
               do i=imin-HALO,imax+HALO-1
!                 Note (KK): [uu|vv]adv defined on U-points (future X-points)
!                            no change of vvadv for eastern closed bdy
!                            hnadv defined on V-points (future T-points)
!                            note that hnadv is not shifted anymore !
!                 KK-TODO: uuadv for northern/southern open bdys ?
                  if (az(i,j+1) .eq. 0) then ! western closed bdy
                     vvadv(i,j,k) = vvadv(i+1,j,k)
                  else if (az(i+1,j+1) .ne. 0) then
                     uuadv(i,j,k) = _HALF_*( uuadv(i,j,k) + uuadv(i+1,j,k) )
                     vvadv(i,j,k) = _HALF_*( vvadv(i,j,k) + vvadv(i+1,j,k) )
                  end if
                  hnadv(i,j,k) = _HALF_*( hnadv(i,j,k) + hnadv(i+1,j,k) )
               end do
#ifndef SLICE_MODEL
            end do
!$OMP END DO NOWAIT
#endif
         end do

!$OMP SINGLE
         phadv => hnadv
!$OMP END SINGLE

      else

!$OMP SINGLE
         phadv => hvn
!$OMP END SINGLE

      end if

!$OMP SINGLE
      if (vel3d_adv_hor.ne.NOADV .and. vel3d_adv_hor.ne.UPSTREAM .and. vel3d_adv_hor.ne.J7) then
!        we need to update fadv3d(imin-HALO:imax+HALO,jmax+HALO)
         call tic(TIM_UVADV3DH)
         call update_3d_halo(fadv3d,fadv3d,av,imin,jmin,imax,jmax,kmax,V_TAG)
         call wait_halo(V_TAG)
         call toc(TIM_UVADV3DH)
      end if
!$OMP END SINGLE

      if (do_numerical_analyses) then
         do k=1,kmax ! calculate square of v-velocity before advection step
!$OMP DO SCHEDULE(RUNTIME)
            do j=jmin-HALO,jmax+HALO
               do i=imin-HALO,imax+HALO
                  work3d(i,j,k) = fadv3d(i,j,k)**2
               end do
            end do
!$OMP END DO NOWAIT
         end do
!$OMP BARRIER
      end if

!$OMP SINGLE
      call do_advection_3d(dt,fadv3d,uuadv,vvadv,wwadv,phuadv,phvadv,phadv,phadv,    &
                           vel3d_adv_split,vel3d_adv_hor,vel3d_adv_ver,_ZERO_,V_TAG, &
                           advres=vvEx)
!$OMP END SINGLE

#ifdef _MOMENTUM_TERMS_
      do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
         do j=jmin,jmax
            do i=imin,imax
               adv_v(i,j,k) = dry_v(i,j) * vvEx(i,j,k)
            end do
         end do
!$OMP END DO NOWAIT
      end do
#endif

      if (do_numerical_analyses) then

!$OMP SINGLE
         call do_advection_3d(dt,work3d,uuadv,vvadv,wwadv,phuadv,phvadv,phadv,phadv,    &
                              vel3d_adv_split,vel3d_adv_hor,vel3d_adv_ver,_ZERO_,V_TAG, &
                              hires=hires)
!$OMP END SINGLE

         do k=1,kmax ! calculate kinetic energy dissipaion rate for v-velocity
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
            do j=jmin,jmax
#endif
               do i=imin,imax
                  work3d(i,j,k) = ( work3d(i,j,k) - fadv3d(i,j,k)**2 ) / dt
               end do
#ifndef SLICE_MODEL
            end do
#endif
!$OMP END DO NOWAIT
         end do

!$OMP BARRIER
!$OMP SINGLE
#ifdef SLICE_MODEL
         work3d(imin:imax,j-1,1:kmax) = work3d(imin:imax,j,1:kmax)
#else
         call update_3d_halo(work3d,work3d,av,imin,jmin,imax,jmax,kmax,V_TAG)
         call wait_halo(V_TAG)
#endif
!$OMP END SINGLE

         do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
            do j=jmin,jmax
#endif
               do i=imin,imax
                  if (az(i,j) .eq. 1) then
                     numdis3d(i,j,k) = numdis3d(i,j,k)                             &
                                       +_HALF_*( work3d(i,j-1,k) + work3d(i,j,k) )
                     numdis2d(i,j) = numdis2d(i,j)                             &
                                    +_HALF_*( work3d(i,j-1,k)*hires(i,j-1,k)   &
                                             +work3d(i,j  ,k)*hires(i,j  ,k) )
                  end if
               end do
#ifndef SLICE_MODEL
            end do
#endif
!$OMP END DO
         end do

      end if

!$OMP END PARALLEL

   end if

   call toc(TIM_UVADV3D)
#ifdef DEBUG
   write(debug,*) 'Leaving uv_advect_3d()'
   write(debug,*)
#endif
   return
   end subroutine uv_advect_3d
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
