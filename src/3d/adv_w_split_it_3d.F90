#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:  adv_w_split_it_3d -  iterated 1D z-advection
!             \label{sec-w-split-it-adv}
!
! !INTERFACE:
   subroutine adv_w_split_it_3d(dt,f,hi,ww,az,splitfac,method)
!
! !DESCRIPTION:
!
! Here, the same one-dimensional advection step as in {\tt w\_split\_adv}
! (see section \ref{sec-w-split-adv} on page \pageref{sec-w-split-adv})
! is applied, but with an
! iteration in time in case that the vertical Courant number exceeds
! unity at any interface of the water column under calculation.
!
! The number of time steps is calculated as
!
! \begin{equation}
! N^{\Delta_t}_{i,j} = \max_k\left\{\mbox{\tt int}\left(C_{i,j,k}+1\right)
! \right\},
! \end{equation}
!
! with the Courant number
!
! \begin{equation}
! C_{i,j,k}=\left|w_{i,j,k}\right|
! \frac{\Delta t}{\frac12 \left(h^n_{i,j,k}+h^n_{i,j,k+1}
! \right)}
! \end{equation}
!
! and the truncation function {\tt int}.
!
! After the number of iterations $N^{\Delta_t}_{i,j}$ is calculated,
! the vertical advection step is calculated $N^{\Delta_t}_{i,j}$
! with a time step of $\Delta t / N^{\Delta_t}_{i,j}$. By doing so,
! it is avoided that the model blows up due to violation
! of the CFL criterium, which could happen fast in case of
! high vertical resolution together with processes such as upwelling
! or fast sinking material. The good thing about this procedure is that
! only the water column is punished by higher numerical load
! in which the potential violation of the CFL criterium
! occurs.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax
   use advection_3d, only: cu
   use advection_3d, only: UPSTREAM_SPLIT,P2,SUPERBEE,MUSCL,P2_PDM
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                        :: dt,splitfac
   REALTYPE,dimension(I3DFIELD),intent(in)    :: ww
   integer,dimension(E2DFIELD),intent(in)     :: az
   integer,intent(in)                         :: method
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(I3DFIELD),intent(inout) :: f,hi
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer         :: i,ii,j,jj,k,kk,it
   REALTYPE        :: hio,c,alpha,beta,x,r,Phi,limit,fu,fc,fd,cmax
   logical         :: READY
   REALTYPE,parameter :: one6th=_ONE_/6
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'adv_w_split_it_3d() # ',Ncall
#endif

! OMP note: Initialization (mem copy) done in serial:
   !cu = _ZERO_

!$OMP PARALLEL DEFAULT(SHARED)                                          &
!$OMP       PRIVATE(i,ii,j,jj,k,kk,it,READY)                            &
!$OMP       PRIVATE(hio,c,alpha,beta,x,r,Phi,limit,fu,fc,fd,cmax)

! OMP TODO: The present loops (j-i-k) gains only a small speedup from
!  threading. Likely, the many jumps in memory limits the performance.

! Calculating w-interface fluxes !
   select case (method)
      case (UPSTREAM_SPLIT)
!$OMP MASTER
         cu(:,:,kmax) = _ZERO_
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(RUNTIME)
         do j=jmin,jmax
            do i=imin,imax
               if (az(i,j) .eq. 1) then
                  cmax= _ZERO_
                  it=1
                  ready=.false.
222               do ii=1,it
                     do k=1,kmax-1
                        !cu(i,j,k) = _ZERO_
                        if (ww(i,j,k) .gt. _ZERO_) then
                           c=ww(i,j,k)/it*dt/(_HALF_*(hi(i,j,k)+hi(i,j,k+1)))
                           if (c .gt. cmax) cmax=c
                           cu(i,j,k)=ww(i,j,k)*f(i,j,k)
                        else
                           c=-ww(i,j,k)/it*dt/(_HALF_*(hi(i,j,k)+hi(i,j,k+1)))
                           if (c .gt. cmax) cmax=c
                           cu(i,j,k)=ww(i,j,k)*f(i,j,k+1)
                        end if
                     end do
                     if (.not. READY) then
                        it=min(200,int(cmax)+1)
#ifdef DEBUG
                        if (it .gt. 1) write(95,*) i,j,it,cmax
#endif
                     end if
                     if ((it .gt. 1) .and. (.not. READY)) then
                        READY=.true.
                        goto 222
                     end if
                     do k=1,kmax   ! Doing a w-advection step
                        hio=hi(i,j,k)
                        hi(i,j,k)=hio-splitfac/float(it)               &
                                  *dt*(ww(i,j,k)-ww(i,j,k-1))
                        f(i,j,k)=(f(i,j,k)*hio-splitfac/float(it)      &
                                  *dt*(cu(i,j,k)-cu(i,j,k-1)))/hi(i,j,k)
                     end do
                  end do
               end if
            end do
         end do
!$OMP END DO
      case ((P2),(Superbee),(MUSCL),(P2_PDM))
!$OMP MASTER
         cu(:,:,0) = _ZERO_
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(RUNTIME)
         do j=jmin,jmax
            do i=imin,imax
               if (az(i,j) .eq. 1) then
                  cmax= _ZERO_
                  it=1
                  ready=.false.
111               do ii=1,it
                     do k=1,kmax
                        cu(i,j,k) = _ZERO_
                        if (ww(i,j,k) .gt. _ZERO_) then
                           if (k.lt.kmax) then
                              c=ww(i,j,k)/float(it)*dt/(_HALF_*(hi(i,j,k)+hi(i,j,k+1)))
                           else
                              c=ww(i,j,k)/float(it)*dt/hi(i,j,k)
                           end if
                           if (c .gt. cmax) cmax=c
                           if (k .gt. 1) then
                              fu=f(i,j,k-1)         ! upstream
                           else
                              fu=f(i,j,k)
                           end if
                           fc=f(i,j,k  )            ! central
                           if (k.lt.kmax) then
                              fd=f(i,j,k+1)         ! downstream
                           else
                              fd=f(i,j,k)           ! downstream
                           end if
                           if (abs(fd-fc) .gt. 1.d-10) then
                              r=(fc-fu)/(fd-fc)     ! slope ratio
                           else
                              r=(fc-fu)*1.d10
                           end if
                        else
                           if (k.lt.kmax) then
                              c=-ww(i,j,k)/float(it)*dt/(_HALF_*(hi(i,j,k)+hi(i,j,k+1)))
                           else
                              c=-ww(i,j,k)/float(it)*dt/hi(i,j,k)
                           end if
                           if (c .gt. cmax) cmax=c
                           if (k .lt. kmax-1) then
                              fu=f(i,j,k+2)         ! upstream
                           else
                              if (k.lt.kmax) then
                                 fu=f(i,j,k+1)
                              else
                                 fu=f(i,j,k)
                              end if
                           end if
                           if (k.lt.kmax) then
                              fc=f(i,j,k+1)            ! central
                           else
                              fc=f(i,j,k)              ! central
                           end if
                           fd=f(i,j,k  )            ! downstream
                           if (abs(fc-fd) .gt. 1.d-10) then
                              r=(fu-fc)/(fc-fd)     ! slope ratio
                           else
                           r=   (fu-fc)*1.d10
                           end if
                        end if
                        select case (method)
                           case ((P2),(P2_PDM))
                              x = one6th*(_ONE_-_TWO_*c)
                              Phi=(_HALF_+x)+(_HALF_-x)*r
                              if (method.eq.P2) then
                                 limit=Phi
                              else
                                 limit=max(_ZERO_,min(Phi,_TWO_/(_ONE_-c),_TWO_*r/(c+1.d-10)))
                              end if
                           case (Superbee)
                              limit=max(_ZERO_, min(_ONE_, _TWO_*r), min(r,_TWO_) )
                           case (MUSCL)
                              limit=max(_ZERO_,min(_TWO_,_TWO_*r,_HALF_*(_ONE_+r)))
                           case default
                              FATAL 'This is not so good - do_advection_3d()'
                              stop 'adv_w_split_it_3d'
                        end select
                        cu(i,j,k)=ww(i,j,k)*(fc+_HALF_*limit*(_ONE_-c)*(fd-fc))
                     end do
                     if (.not. READY) then
                        it=min(200,int(cmax)+1)
#ifdef DEBUG
                        if (it .gt. 1) write(95,*) i,j,it,cmax
#endif
                     end if
                     if ((it .gt. 1) .and. (.not. READY)) then
                        READY=.true.
                        goto 111
                     end if

                     do k=1,kmax   ! Doing a w-advection step
                        hio=hi(i,j,k)
                        hi(i,j,k)=hio                  &
                                 -splitfac/float(it)*dt*(ww(i,j,k)-ww(i,j,k-1))
                        f(i,j,k)=(f(i,j,k)*hio         &
                                -splitfac/float(it)*dt*(cu(i,j,k)-cu(i,j,k-1)))/hi(i,j,k)
                     end do
                  end do
               end if
            end do
         end do
!$OMP END DO
   end select
!$OMP END PARALLEL

#ifdef DEBUG
   write(debug,*) 'Leaving adv_w_split_it_3d()'
   write(debug,*)
#endif
   return
   end subroutine adv_w_split_it_3d

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
