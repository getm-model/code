#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:  adv_w_split_3d - 1D z-advection \label{sec-w-split-adv}
!
! !INTERFACE:
   subroutine adv_w_split_3d(dt,f,hi,adv3d,ww,ho,                       &
                             az,au,av,splitfac,scheme,onestep_finalise)
!  Note (KK): Keep in sync with interface in advection.F90
!
! !DESCRIPTION:
!
! Here, the $z$-directional split 1D advection step is executed
! with a number of options for the numerical scheme. The basic
! advection equation is accompanied by an fractional step
! for the continuity equation and both equations look as follows:
!
! \begin{equation}\label{adv_w_step}
! h^n_{i,j,k} c^n_{i,j,k} =
! h^o_{i,j,k} c^o_{i,j,k}
! - \Delta t
! \left(w_{i,j,k}\tilde c^w_{i,j,k}-w_{i,j,k-1}\tilde c^w_{i,j,k-1}\right),
! \end{equation}
!
! with the 1D continuity equation
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
! {\tt u\_split\_adv}, see section \ref{sec-u-split-adv} on
! page \pageref{sec-u-split-adv}.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax
   use advection, only: UPSTREAM,P2,SUPERBEE,MUSCL,P2_PDM
   use advection_3d, only: itersmax_adv
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                        :: dt,splitfac
   REALTYPE,dimension(I3DFIELD),intent(in)    :: ww,ho
   integer,dimension(E2DFIELD),intent(in)     :: az,au,av
   integer,intent(in)                         :: scheme
   logical,intent(in),optional                :: onestep_finalise
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(I3DFIELD),intent(inout) :: f,hi,adv3d
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   logical            :: iterate
   integer            :: i,j,k,it,iters,rc
   REALTYPE           :: hio,advn,cfl,x,r,Phi,limit,fu,fc,fd,itersm1
   REALTYPE,dimension(:),allocatable :: flux1d
   REALTYPE,parameter :: one6th=_ONE_/6
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'adv_w_split_3d() # ',Ncall
#endif

   iterate=.false.
   if (.not.present(onestep_finalise)) then
      if (itersmax_adv .gt. 1) iterate=.true.
   end if

   if (.not. iterate) then
#ifdef NO_BAROTROPIC
      stop 'adv_w_split_3d: do enable iterations with compiler option NO_BAROTROPIC'
#endif
      iters = 1
      itersm1 = _ONE_
   end if

!$OMP PARALLEL DEFAULT(SHARED)                                        &
!$OMP FIRSTPRIVATE(iters,itersm1)                                     &
!$OMP PRIVATE(i,j,k,it,rc,flux1d,hio,advn,cfl,x,r,Phi,limit,fu,fc,fd)

!  Each thread allocates its own HEAP storage:
   allocate(flux1d(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'adv_w_split_3d: Error allocating memory (flux1d)'

   flux1d(0) = _ZERO_
   flux1d(kmax) = _ZERO_

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j).eq.1 .or. (az(i,j).eq.2 .and. (    au(i-1,j  ).eq.1 &
                                                    .or.au(i  ,j  ).eq.1 &
                                                    .or.av(i  ,j-1).eq.1 &
                                                    .or.av(i  ,j  ).eq.1 ))) then
!           Note (KK): exclude tracer open bdy cells but
!                      include all velocity open bdys
            if (iterate) then
!              estimate number of iterations by maximum cfl number in water column
               cfl = _ZERO_
               do k=1,kmax-1
                  cfl = max(cfl,abs(splitfac*ww(i,j,k))*dt/(_HALF_*(hi(i,j,k)+hi(i,j,k+1))))
               end do
               !iters =  min(200,int(cfl)+1) !original
               iters = min(max(1,ceiling(cfl)),itersmax_adv)
               itersm1 = _ONE_/iters
            end if
            do it=1,iters
!              Calculating w-interface fluxes !
               do k=1,kmax-1
                  if (ww(i,j,k) .gt. _ZERO_) then
                     fc = f(i,j,k  )               ! central
                     if (scheme .ne. UPSTREAM) then
                        cfl = itersm1*splitfac*ww(i,j,k)*dt/(_HALF_*(hi(i,j,k)+hi(i,j,k+1)))
                        if (k .eq. 1) then
                           fu = f(i,j,1)
                        else
                           fu = f(i,j,k-1)         ! upstream
                        end if
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
                        cfl = -itersm1*splitfac*ww(i,j,k)*dt/(_HALF_*(hi(i,j,k)+hi(i,j,k+1)))
                        if (k .eq. kmax-1) then
                           fu = f(i,j,kmax)
                        else
                           fu = f(i,j,k+2)         ! upstream
                        end if
                        fd = f(i,j,k  )            ! downstream
                        if (abs(fc-fd) .gt. 1.d-10) then
                           r = (fu-fc)/(fc-fd)     ! slope ratio
                        else
                           r = (fu-fc)*1.d10
                        end if
                     end if
                  end if
                  flux1d(k) = ww(i,j,k)*fc
                  if (scheme .ne. UPSTREAM) then
                     select case (scheme)
                        case ((P2),(P2_PDM))
                           x = one6th*(_ONE_-_TWO_*cfl)
                           Phi = (_HALF_+x) + (_HALF_-x)*r
                           if (scheme.eq.P2) then
                              limit = Phi
                           else
                              limit = max(_ZERO_,min(Phi,_TWO_/(_ONE_-cfl),_TWO_*r/(cfl+1.d-10)))
                           end if
                        case (SUPERBEE)
                           limit = max(_ZERO_,min(_ONE_, _TWO_*r),min(r,_TWO_))
                        case (MUSCL)
                           limit = max(_ZERO_,min(_TWO_,_TWO_*r,_HALF_*(_ONE_+r)))
                        case default
                           stop 'adv_w_split_3d: invalid scheme'
                     end select
                     flux1d(k) = flux1d(k) + ww(i,j,k)*_HALF_*limit*(_ONE_-cfl)*(fd-fc)
                  end if
               end do

!              Doing the w-advection step
               do k=1,kmax
                  hio = hi(i,j,k)
                  hi(i,j,k) = hio - itersm1*splitfac*dt*(ww(i,j,k  )-ww(i,j,k-1))
                  advn = itersm1*splitfac*(flux1d(k  )-flux1d(k-1))
                  adv3d(i,j,k) = adv3d(i,j,k) + advn
                  if (present(onestep_finalise)) then
!                    Note (KK): do not update f in case of onestep_finalise=.false. !!!
                     if (onestep_finalise) then
                        f(i,j,k) = ( ho(i,j,k)*f(i,j,k) - dt*adv3d(i,j,k) ) / hi(i,j,k)
                     end if
                  else
                     f(i,j,k) = ( hio*f(i,j,k) - dt*advn ) / hi(i,j,k)
                  end if
               end do
            end do
         end if
      end do
   end do
!$OMP END DO

!  Each thread must deallocate its own HEAP storage:
   deallocate(flux1d,stat=rc)
   if (rc /= 0) stop 'adv_w_split_3d: Error deallocating memory (flux1d)'

!$OMP END PARALLEL

#ifdef DEBUG
   write(debug,*) 'Leaving adv_w_split_3d()'
   write(debug,*)
#endif
   return
   end subroutine adv_w_split_3d
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
