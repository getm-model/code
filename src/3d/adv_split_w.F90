#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:  adv_split_w - vertical advection of 3D quantities \label{sec-w-split-adv}
!
! !INTERFACE:
   subroutine adv_split_w(dt,f,hi,adv3d,ww,ho,              &
                          splitfac,scheme,action,tag_3d,az, &
                          mask_finalise)
!  Note (KK): Keep in sync with interface in advection_3d.F90
!
! !DESCRIPTION:
!
! Executes an advection step in vertical direction. The 1D advection
! equation
!
! \begin{equation}\label{adv_w_step}
! h^n_{i,j,k} c^n_{i,j,k} =
! h^o_{i,j,k} c^o_{i,j,k}
! - \Delta t
! \left(w_{i,j,k}\tilde c^w_{i,j,k}-w_{i,j,k-1}\tilde c^w_{i,j,k-1}\right),
! \end{equation}
!
! is accompanied by an fractional step for the 1D continuity equation
!
! \begin{equation}\label{adv_w_step_h}
! h^n_{i,j,k}  =
! h^o_{i,j,k}
! - \Delta t
! \left(w_{i,j,k}\tilde -w_{i,j,k-1}\right).
! \end{equation}
!
! Here, $n$ and $o$ denote values before and after this operation,
! respectively, $n$ denote intermediate values when other
! 1D advection steps come after this and $o$ denotes intermediate
! values when other 1D advection steps came before this.
!
! The interfacial fluxes $\tilde c^w_{i,j,k}$ are calculated by means of
! monotone and non-monotone schemes which are described in detail in
! section \ref{sec-u-split-adv} on page \pageref{sec-u-split-adv}.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax
   use advection, only: adv_interfacial_reconstruction
   use advection, only: NOADV,UPSTREAM
   use advection, only: NOSPLIT_FINALISE,SPLIT_UPDATE
   use advection_3d, only: adv_ver_iterations,W_TAG
   use halo_zones, only: U_TAG,V_TAG
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                             :: dt,splitfac
   REALTYPE,dimension(I3DFIELD),intent(in)         :: ww,ho
   integer,intent(in)                              :: scheme,action,tag_3d
   integer,dimension(E2DFIELD),intent(in)          :: az
   logical,dimension(E2DFIELD),intent(in),optional :: mask_finalise
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(I3DFIELD),target,intent(inout) :: f,hi,adv3d
!
! !LOCAL VARIABLES:
   logical            :: iterate,use_limiter
   integer            :: i,j,k,kshift,it,iters,iters_new,rc
   REALTYPE           :: itersm1,dti,dtik,hio,advn,fuu,fu,fd,splitfack
   REALTYPE,dimension(:),allocatable        :: wflux
   REALTYPE,dimension(:),allocatable,target :: cfl0
   REALTYPE,dimension(:),pointer            :: hiaux,advaux,faux,cfls
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'adv_split_w() # ',Ncall
#endif
#ifdef SLICE_MODEL
   j = jmax/2 ! this MUST NOT be changed!!!
#endif

   if (tag_3d .eq. W_TAG) then
      kshift = 1
   else
      kshift = 0
   end if

   use_limiter = .false.
   dti = splitfac*dt

   iterate=.false.
   if (action .eq. SPLIT_UPDATE) then
      if (adv_ver_iterations .gt. 1) iterate=.true.
   end if

   if (.not. iterate) then
#ifdef NO_BAROTROPIC
      stop 'adv_split_w: do enable iterations with compiler option NO_BAROTROPIC'
#endif
   end if

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j,use_limiter)                             &
!$OMP          PRIVATE(rc,wflux,hiaux,advaux,faux,cfl0,cfls)           &
!$OMP          PRIVATE(itersm1,dtik,splitfack)                         &
!$OMP          PRIVATE(i,k,it,iters,iters_new,hio,advn,fuu,fu,fd)


   if (scheme .ne. NOADV) then

!     Each thread allocates its own HEAP storage:
      allocate(wflux(0:kmax),stat=rc)    ! work array
      if (rc /= 0) stop 'adv_split_w: Error allocating memory (wflux)'

#ifdef _POINTER_REMAP_
      if (iterate) then
#endif
         allocate(hiaux(0:kmax),stat=rc)    ! work array
         if (rc /= 0) stop 'adv_split_w: Error allocating memory (hiaux)'

         allocate(advaux(0:kmax),stat=rc)    ! work array
         if (rc /= 0) stop 'adv_split_w: Error allocating memory (advaux)'

         allocate(faux(0:kmax),stat=rc)    ! work array
         if (rc /= 0) stop 'adv_split_w: Error allocating memory (faux)'
#ifdef _POINTER_REMAP_
      end if
#endif

      if (scheme.ne.UPSTREAM .or. iterate) then
         allocate(cfl0(0:kmax),stat=rc)    ! work array
         if (rc /= 0) stop 'adv_split_w: Error allocating memory (cfl0)'

         if (iterate) then
            allocate(cfls(0:kmax),stat=rc)    ! work array
            if (rc /= 0) stop 'adv_split_w: Error allocating memory (cfls)'
         else
            cfls => cfl0
         end if
      end if

      wflux(0) = _ZERO_
      wflux(kmax) = _ZERO_

!     Note (KK): as long as h[u|v]n([i|j]max+HALO) are trash (SMALL)
!                they have to be excluded from the loop to avoid
!                unnecessary iterations and warnings

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin-HALO,jmax+HALO-1
#endif
         do i=imin-HALO,imax+HALO-1
            if (az(i,j) .eq. 1) then
!              Note (KK): exclude vertical advection of normal velocity at open bdys

               if (scheme.ne.UPSTREAM .or. iterate) then
                  do k=1-kshift,kmax-1
                     cfl0(k) = abs(ww(i,j,k))*dti/(_HALF_*(hi(i,j,k)+hi(i,j,k+1)))
                  end do
               end if

               iters = 1
               itersm1 = _ONE_
               dtik = dti
               splitfack = splitfac

               it = 0

               do while (it .lt. iters)

                  if (it .eq. 0) then
#ifdef _POINTER_REMAP_
                     if (iterate) then
#endif
                        hiaux  = hi   (i,j,:)
                        advaux = adv3d(i,j,:)
                        faux   = f    (i,j,:)
#ifdef _POINTER_REMAP_
                     else
                        hiaux (0:) => hi   (i,j,:)
                        advaux(0:) => adv3d(i,j,:)
                        faux  (0:) => f    (i,j,:)
                     end if
#endif
                  end if
                  it = it + 1

                  if (iterate) then
                     if (it .eq. 1) then
                        cfls = cfl0
                     else if (scheme.ne.UPSTREAM .or. iters.lt.adv_ver_iterations) then
                        do k=1-kshift,kmax-1
                           cfls(k) = abs(ww(i,j,k))*dti/(_HALF_*(hiaux(k)+hiaux(k+1)))
                        end do
                     end if
                     if (iters .lt. adv_ver_iterations) then
!                       estimate number of iterations by maximum cfl number in water column
                        iters_new = max(1,maxval(ceiling(cfls(1-kshift:kmax-1))))
                        if (iters_new .gt. iters) then
                           if (iters_new .gt. adv_ver_iterations) then
!$OMP CRITICAL
                              STDERR 'adv_split_w: too many iterations needed at'
                              STDERR 'i=',i,' j=',j,':',iters_new
!$OMP END CRITICAL
                              iters = adv_ver_iterations
                           else
                              iters = iters_new
                           end if
                           itersm1 = _ONE_ / iters
                           dtik = dti * itersm1
                           splitfack = splitfac * itersm1
                           if (it .gt. 1) then
!$OMP CRITICAL
                              STDERR 'adv_split_w: restart iterations during it=',it
                              STDERR 'i=',i,' j=',j,':',iters
!$OMP END CRITICAL
                              it = 0
                              cycle
                           end if
                        end if
                     end if
                  end if

!                 Calculating w-interface fluxes !
                  do k=1-kshift,kmax-1
!                    Note (KK): overwrite zero flux at k=0 in case of W_TAG
                     if (ww(i,j,k) .gt. _ZERO_) then
                        fu = faux(k)               ! central
                        if (scheme .ne. UPSTREAM) then
!                          also fall back to upstream near boundaries
                           use_limiter = (k .gt. 1-kshift)
                        end if
                        if (use_limiter) then
                           fuu = faux(k-1)            ! upstream
                           fd  = faux(k+1)            ! downstream
                        end if
                     else
                        fu = faux(k+1)               ! central
                        if (scheme .ne. UPSTREAM) then
!                          also fall back to upstream near boundaries
                           use_limiter = (k .lt. kmax-1)
                        end if
                        if (use_limiter) then
                           fuu = faux(k+2)            ! upstream
                           fd  = faux(k  )            ! downstream
                        end if
                     end if
                     if (use_limiter) then
                        fu = adv_interfacial_reconstruction(scheme,cfls(k)*itersm1,fuu,fu,fd)
                     end if
                     wflux(k) = ww(i,j,k)*fu
                  end do
                  do k=1,kmax-kshift
!                    Note (KK): in case of W_TAG do not advect at k=kmax
                     hio = hiaux(k)
                     hiaux(k) = hio - dtik*(ww(i,j,k  )-ww(i,j,k-1))
                     advn = splitfack*(wflux(k  )-wflux(k-1))
                     advaux(k) = advaux(k) + advn
                     if (action .eq. SPLIT_UPDATE) then
                        faux(k) = ( hio*faux(k) - dt*advn ) / hiaux(k)
                     end if
                  end do
               end do
#ifdef _POINTER_REMAP_
               if (iterate) then
#endif
                  hi   (i,j,1:kmax-kshift) = hiaux (1:kmax-kshift)
                  adv3d(i,j,1:kmax-kshift) = advaux(1:kmax-kshift)
#ifndef _POINTER_REMAP_
                  if (action .eq. SPLIT_UPDATE) then
#endif
                  f    (i,j,1:kmax-kshift) = faux  (1:kmax-kshift)
#ifndef _POINTER_REMAP_
                  end if
#endif
#ifdef _POINTER_REMAP_
               end if
#endif
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO

!     Each thread must deallocate its own HEAP storage:
      deallocate(wflux,stat=rc)
      if (rc /= 0) stop 'adv_split_w: Error deallocating memory (wflux)'

#ifdef _POINTER_REMAP_
      if (iterate) then
#endif
         deallocate(hiaux,stat=rc)
         if (rc /= 0) stop 'adv_split_w: Error deallocating memory (hiaux)'
         deallocate(advaux,stat=rc)
         if (rc /= 0) stop 'adv_split_w: Error deallocating memory (advaux)'
         deallocate(faux,stat=rc)
         if (rc /= 0) stop 'adv_split_w: Error deallocating memory (faux)'
#ifdef _POINTER_REMAP_
      end if
#endif
      if (scheme.ne.UPSTREAM .or. iterate) then
         deallocate(cfl0,stat=rc)
         if (rc /= 0) stop 'adv_split_w: Error deallocating memory (cfl0)'

         if (iterate) then
            deallocate(cfls,stat=rc)
            if (rc /= 0) stop 'adv_split_w: Error deallocating memory (cfls)'
         end if
      end if


   end if

   if (action .eq. NOSPLIT_FINALISE) then
      do k=1,kmax-kshift
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
         do j=jmin-HALO,jmax+HALO
#endif
            do i=imin-HALO,imax+HALO
               if (mask_finalise(i,j)) then
                  f(i,j,k) = ( ho(i,j,k)*f(i,j,k) - dt*adv3d(i,j,k) ) / hi(i,j,k)
               end if
            end do
#ifndef SLICE_MODEL
         end do
#endif
!$OMP END DO NOWAIT
      end do
   end if

!$OMP END PARALLEL

#ifdef DEBUG
   write(debug,*) 'Leaving adv_split_w()'
   write(debug,*)
#endif
   return
   end subroutine adv_split_w
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
