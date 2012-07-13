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
   use advection, only: adv_tvd_limiter
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
   REALTYPE,dimension(I3DFIELD),intent(inout)      :: f,hi,adv3d
!
! !LOCAL VARIABLES:
   logical            :: iterate,use_limiter
   integer            :: i,j,k,kshift,it,iters,rc
   REALTYPE           :: dti,dtik,hio,advn,cfl,r,limit,fu,fc,fd,splitfack
   REALTYPE,dimension(:),allocatable :: wflux
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
   j = jmax/2 ! thus MUST NOT be changed!!!
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
      iters = 1
      dtik = dti
      splitfack = splitfac
   end if

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j,use_limiter,iters,dtik,splitfack)        &
!$OMP          PRIVATE(i,k,it,rc,wflux,hio,advn,cfl,r,limit,fu,fc,fd)

   if (scheme .ne. NOADV) then

!     Each thread allocates its own HEAP storage:
      allocate(wflux(0:kmax),stat=rc)    ! work array
      if (rc /= 0) stop 'adv_split_w: Error allocating memory (wflux)'

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
               if (iterate) then
!                 estimate number of iterations by maximum cfl number in water column
                  cfl = _ZERO_
                  do k=1-kshift,kmax-1
                     cfl = max(cfl,abs(ww(i,j,k))*dti/(_HALF_*(hi(i,j,k)+hi(i,j,k+1))))
                  end do
                  iters = max(1,ceiling(cfl))
                  if (iters .gt. adv_ver_iterations) then
!$OMP CRITICAL
                     STDERR 'adv_split_w: too many iterations needed at'
                     STDERR 'i=',i,' j=',j,':',iters
!$OMP END CRITICAL
                  end if
                  !iters =  min(200,int(cfl)+1) !original
                  iters = min(iters,adv_ver_iterations)
                  dtik = dti/iters
                  splitfack = splitfac/iters
               end if
               do it=1,iters
!                 Calculating w-interface fluxes !
                  do k=1-kshift,kmax-1
!                    Note (KK): overwrite zero flux at k=0 in case of W_TAG
                     if (ww(i,j,k) .gt. _ZERO_) then
                        fc = f(i,j,k  )               ! central
                        if (scheme .ne. UPSTREAM) then
!                          also fall back to upstream near boundaries
                           use_limiter = (k .gt. 1-kshift)
                        end if
                        if (use_limiter) then
                           cfl = ww(i,j,k)*dtik/(_HALF_*(hi(i,j,k)+hi(i,j,k+1)))
                           fu = f(i,j,k-1)            ! upstream
                           fd = f(i,j,k+1)            ! downstream
                           if (abs(fd-fc) .gt. 1.d-10) then
                              r = (fc-fu)/(fd-fc)     ! slope ratio
                           else
                              r = (fc-fu)*1.d10
                           end if
                        end if
                     else
                        fc = f(i,j,k+1)               ! central
                        if (scheme .ne. UPSTREAM) then
!                          also fall back to upstream near boundaries
                           use_limiter = (k .lt. kmax-1)
                        end if
                        if (use_limiter) then
                           cfl = -ww(i,j,k)*dtik/(_HALF_*(hi(i,j,k)+hi(i,j,k+1)))
                           fu = f(i,j,k+2)            ! upstream
                           fd = f(i,j,k  )            ! downstream
                           if (abs(fc-fd) .gt. 1.d-10) then
                              r = (fu-fc)/(fc-fd)     ! slope ratio
                           else
                              r = (fu-fc)*1.d10
                           end if
                        end if
                     end if
                     wflux(k) = ww(i,j,k)*fc
                     if (use_limiter) then
                        limit = adv_tvd_limiter(scheme,cfl,r)
                        wflux(k) = wflux(k) + ww(i,j,k)*_HALF_*limit*(_ONE_-cfl)*(fd-fc)
                     end if
                  end do
                  do k=1,kmax-kshift
!                    Note (KK): in case of W_TAG do not advect at k=kmax
                     hio = hi(i,j,k)
                     hi(i,j,k) = hio - dtik*(ww(i,j,k  )-ww(i,j,k-1))
                     advn = splitfack*(wflux(k  )-wflux(k-1))
                     adv3d(i,j,k) = adv3d(i,j,k) + advn
                     if (action .eq. SPLIT_UPDATE) then
                        f(i,j,k) = ( hio*f(i,j,k) - dt*advn ) / hi(i,j,k)
                     end if
                  end do
               end do
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO

!     Each thread must deallocate its own HEAP storage:
      deallocate(wflux,stat=rc)
      if (rc /= 0) stop 'adv_split_w: Error deallocating memory (wflux)'

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
