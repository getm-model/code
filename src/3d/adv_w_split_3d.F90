#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:  adv_w_split_3d - 1D z-advection \label{sec-w-split-adv}
!
! !INTERFACE:
   subroutine adv_w_split_3d(dt,f,ww,az,splitfac,method)
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
   use advection_3d, only: hi,hio,cu
   use advection_3d, only: UPSTREAM_SPLIT,P2,SUPERBEE,MUSCL,P2_PDM
   use advection_3d, only: ONE6TH
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE , intent(in)               :: ww(I3DFIELD),dt
   integer , intent(in)                :: az(E2DFIELD)
   REALTYPE, intent(in)                :: splitfac
   integer, intent(in)                 :: method
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout)             :: f(I3DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer         :: i,j,k
   REALTYPE        :: c,x,r,Phi,limit,fu,fc,fd
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'adv_w_split_3d() # ',Ncall
#endif
#ifdef NO_BAROTROPIC
   STDERR 'Do not use adv_w_split_3d() with the compiler option NO_BAROTROPIC'
   STDERR 'Use adv_w_split_it_3d() instead, choose option ITERATE_VERT_ADV'
   stop 'adv_w_split_3d()'
#endif

   cu(:,:,kmax)=_ZERO_

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,c,x,r,Phi,limit,fu,fc,fd)

! Calculating w-interface fluxes !

   select case (method)
      case (UPSTREAM_SPLIT)
         do k=1,kmax-1
!$OMP DO SCHEDULE(RUNTIME)
            do j=jmin,jmax
               do i=imin,imax
                  if (az(i,j) .eq. 1) then
                     if (ww(i,j,k) .gt. _ZERO_) then
                        cu(i,j,k)=ww(i,j,k)*f(i,j,k)
                     else
                        cu(i,j,k)=ww(i,j,k)*f(i,j,k+1)
                     end if
                  else
                     cu(i,j,k) = _ZERO_
                  end if
               end do
            end do
!$OMP END DO NOWAIT
         end do
      case ((P2),(Superbee),(MUSCL),(P2_PDM))
         do k=1,kmax-1
!$OMP DO SCHEDULE(RUNTIME)
            do j=jmin,jmax
               do i=imin,imax
                  cu(i,j,k) = _ZERO_
                  if (az(i,j) .eq. 1) then
                     if (ww(i,j,k) .gt. _ZERO_) then
                        c=ww(i,j,k)*dt/(_HALF_*(hi(i,j,k)+hi(i,j,k+1)))
                        if (k .gt. 1) then
                           fu=f(i,j,k-1)         ! upstream
                        else
                           fu=f(i,j,k)
                        end if
                        fc=f(i,j,k  )            ! central
                        fd=f(i,j,k+1)            ! downstream
                        if (abs(fd-fc) .gt. 1.d-10) then
                           r=(fc-fu)/(fd-fc)     ! slope ratio
                        else
                           r=(fc-fu)*1.d10
                        end if
                     else
                        c=-ww(i,j,k)*dt/(_HALF_*(hi(i,j,k)+hi(i,j,k+1)))
                        if (k .lt. kmax-1) then
                           fu=f(i,j,k+2)         ! upstream
                        else
                           fu=f(i,j,k+1)
                        end if
                        fc=f(i,j,k+1)            ! central
                        fd=f(i,j,k  )            ! downstream
                        if (abs(fc-fd) .gt. 1.d-10) then
                           r=(fu-fc)/(fc-fd)     ! slope ratio
                        else
                           r=(fu-fc)*1.d10
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
                           stop 'adv_w_split_3d'
                     end select
                     cu(i,j,k)=ww(i,j,k)*(fc+_HALF_*limit*(_ONE_-c)*(fd-fc))
                  else
                     cu(i,j,k) = _ZERO_
                  end if
               end do
            end do
!$OMP END DO NOWAIT
         end do
   end select

!$OMP BARRIER
   do k=1,kmax   ! Doing a w-advection step
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin,imax
            if (az(i,j) .eq. 1) then
               hio(i,j,k)=hi(i,j,k)
               hi(i,j,k)=hio(i,j,k)-splitfac*dt*(ww(i,j,k)-ww(i,j,k-1))
               f(i,j,k)=(f(i,j,k)*hio(i,j,k)-        &
                         splitfac*dt*(cu(i,j,k)-cu(i,j,k-1)))/hi(i,j,k)
            end if
         end do
      end do
!$OMP END DO
   end do
!$OMP END PARALLEL

#ifdef DEBUG
   write(debug,*) 'Leaving adv_w_split_3d()'
   write(debug,*)
#endif
   return
   end subroutine adv_w_split_3d

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
