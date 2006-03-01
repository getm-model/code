!$Id: w_split_it_adv.F90,v 1.3 2006-03-01 15:54:08 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:  w_split_it_adv -  iterated 1D z-advection 
!             \label{sec-w-split-it-adv}
!
! !INTERFACE:
   subroutine w_split_it_adv(dt,f,ww,az,splitfac,method)
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
! 
! !USES:
   use domain, only: imin,imax,jmin,jmax
   use domain, only: iimin,iimax,jjmin,jjmax,kmax
   use advection_3d, only: hi,hio,cu
   use advection_3d, only: UPSTREAM_SPLIT,P2,SUPERBEE,MUSCL,P2_PDM
   use advection_3d, only: ONE,TWO,ONE6TH
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
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer         :: i,ii,j,jj,k,kk,it
   REALTYPE        :: c,alpha,beta,x,r,Phi,limit,fu,fc,fd,cmax
   logical         :: READY
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'w_split_it_adv() # ',Ncall
#endif

   cu = _ZERO_

! Calculating w-interface fluxes !

   select case (method)
      case (UPSTREAM_SPLIT)
         do j=jjmin,jjmax
            do i=iimin,iimax
               if (az(i,j) .eq. 1) then
                  do k=1,kmax-1
                     cu(i,j,k) = _ZERO_
                     if (ww(i,j,k) .gt. _ZERO_) then
                        cu(i,j,k)=ww(i,j,k)*f(i,j,k)
                     else
                        cu(i,j,k)=ww(i,j,k)*f(i,j,k+1)
                     end if
                  end do
                  do k=1,kmax   ! Doing a w-advection step
                     hio(i,j,k)=hi(i,j,k)
                     hi(i,j,k)=hio(i,j,k)-splitfac*dt*(ww(i,j,k)-ww(i,j,k-1))
                     f(i,j,k)=(f(i,j,k)*hio(i,j,k)-        &
                               splitfac*dt*(cu(i,j,k)-cu(i,j,k-1)))/hi(i,j,k)
                  end do
               end if
            end do
         end do
      case ((P2),(Superbee),(MUSCL),(P2_PDM))
         do j=jjmin,jjmax
            do i=iimin,iimax
               if (az(i,j) .eq. 1) then
                  cmax= _ZERO_
                  it=1
                  ready=.false.
111               do ii=1,it
                     do k=1,kmax-1
                        cu(i,j,k) = _ZERO_
                        if (ww(i,j,k) .gt. _ZERO_) then
                           c=ww(i,j,k)/float(it)*dt/(0.5*(hi(i,j,k)+hi(i,j,k+1)))
                           if (c .gt. cmax) cmax=c
                           if (k .gt. 1) then
                              fu=f(i,j,k-1)         ! upstream
                           else
                              fu=f(i,j,k)
                           end if
                           fc=f(i,j,k  )            ! central
                           fd=f(i,j,k+1)            ! downstream
                           if (abs(fd-fc) .gt. 1e-10) then
                              r=(fc-fu)/(fd-fc)     ! slope ratio
                           else
                              r=(fc-fu)*1.e10
                           end if
                        else
                           c=-ww(i,j,k)/float(it)*dt/(0.5*(hi(i,j,k)+hi(i,j,k+1)))
                           if (c .gt. cmax) cmax=c
                           if (k .lt. kmax-1) then
                              fu=f(i,j,k+2)         ! upstream
                           else
                              fu=f(i,j,k+1)
                           end if
                           fc=f(i,j,k+1)            ! central
                           fd=f(i,j,k  )            ! downstream
                           if (abs(fc-fd) .gt. 1e-10) then
                              r=(fu-fc)/(fc-fd)     ! slope ratio
                           else
                           r=   (fu-fc)*1.e10
                           end if
                        end if
                        x = one6th*(ONE-TWO*c)
                        Phi=(0.5+x)+(0.5-x)*r
                        select case (method)
                           case ((P2),(P2_PDM))
                              x = one6th*(ONE-TWO*c)
                              Phi=(0.5+x)+(0.5-x)*r
                              if (method.eq.P2) then
                                 limit=Phi
                              else
                                 limit=max(_ZERO_,min(Phi,2./(1.-c),2.*r/(c+1.e-10)))
                              end if
                           case (Superbee)
                              limit=max(_ZERO_, min(ONE, TWO*r), min(r,TWO) )
                           case (MUSCL)
                              limit=max(_ZERO_,min(TWO,TWO*r,0.5*(ONE+r)))
                           case default
                              FATAL 'This is not so good - do_advection_3d()'
                              stop 'w_split_it_adv'
                        end select
                        cu(i,j,k)=ww(i,j,k)*(fc+0.5*limit*(1-c)*(fd-fc))
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
                        hio(i,j,k)=hi(i,j,k)
                        hi(i,j,k)=hio(i,j,k)                  &
                                 -splitfac/float(it)*dt*(ww(i,j,k)-ww(i,j,k-1))
                        f(i,j,k)=(f(i,j,k)*hio(i,j,k)         &
                                -splitfac/float(it)*dt*(cu(i,j,k)-cu(i,j,k-1)))/hi(i,j,k)
                     end do
                  end do
               end if
            end do
         end do
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving w_split_it_adv()'
   write(debug,*)
#endif
   return
   end subroutine w_split_it_adv

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
