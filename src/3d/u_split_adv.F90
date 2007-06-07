!$Id: u_split_adv.F90,v 1.5 2007-06-07 10:25:19 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:  u_split_adv - 1D x-advection \label{sec-u-split-adv} 
!
! !INTERFACE:
   subroutine u_split_adv(dt,f,uu,hun, &
                          delxu,delyu,area_inv,au,splitfac,method,az,AH)
! !DESCRIPTION:
!
! Here, the $x$-directional split 1D advection step is executed 
! with a number of options for the numerical scheme. The basic 
! advection equation is accompanied by an fractional step
! for the continuity equation and both equations look as follows:
!
! \begin{equation}\label{adv_u_step}
! h^n_{i,j,k} c^n_{i,j,k} =
! h^o_{i,j,k} c^o_{i,j,k}
! \displaystyle
! - \Delta t 
! \frac{
! p_{i,j,k}\tilde c^u_{i,j,k}\Delta y^u_{i,j}-
! p_{i-1,j,k}\tilde c^u_{i-1,j,k}\Delta y^u_{i-1,j}
! }{\Delta x^c_{i,j}\Delta y^c_{i,j}},
! \end{equation}
!
! with the 1D continuity equation
!
! \begin{equation}\label{adv_u_step_h}
! h^n_{i,j,k} =
! h^o_{i,j,k} 
! \displaystyle
! - \Delta t 
! \frac{
! p_{i,j,k}\Delta y^u_{i,j}-
! p_{i-1,j,k}\Delta y^u_{i-1,j}
! }{\Delta x^c_{i,j}\Delta y^c_{i,j}}
! \end{equation}
!
! Here, $n$ and $o$ denote values before and after this operation,
! respectively, $n$ denote intermediate values when other
! 1D advection steps come after this and $o$ denotes intermediate
! values when other 1D advection steps came before this.
! Furthermore, when this $x$-directional split step is repeated
! during the total time step (Strang splitting), the time step $\Delta t$
! denotes a fraction of the full time step.
!
! The interfacial fluxes $\tilde c^u_{i,j,k}$ are calculated by means of
! monotone Total Variation Diminishing (TVD), the first-order monotone
! upstream and the (non-monotone)
! unlimited third-order polynomial scheme
! according to:
!
! \begin{equation}\label{LaxWendroffForm}
! \tilde c_{i,j,k}=
! \left\{
! \begin{array}{ll}
! \left(c_{i,j,k}+\frac12 \tilde c_{i,j,k}^+ (1-|C_{i,j,k}|)
! (c_{i+1,j,k}-c_{i,j,k})\right) & \mbox{ for } p_{i,j,k} \geq 0, \\ \\
! \left(c_{i+1,j,k}+\frac12 \tilde c_{i,j,k}^- (1-|C_{i,j,k}|)
! (c_{i,j,k}-c_{i+1,j,k})\right) & \mbox{ else, }
! \end{array}
! \right.
! \end{equation}
!
! with the Courant number $C_{i,j,k}=u_{i,j,k}\Delta t/\Delta x$ and
!
! \begin{equation}\label{alphabeta}
! \tilde c_{i,j,k}^+=\alpha_{i,j,k}+\beta_{i,j,k}r^+_{i,j,k}, \quad
! \tilde c_{i,j,k}^-=\alpha_{i,j,k}+\beta_{i,j,k}r^-_{i,j,k},
! \end{equation}
! 
! where
! 
! \begin{equation}
! \alpha_{i,j,k}=\frac12 +\frac16(1-2|C_{i,j,k}|),\quad
! \beta _{i,j,k}=\frac12 -\frac16(1-2|C_{i,j,k}|),
! \end{equation}
! 
! and
! 
! \begin{equation}
! r^+_{i,j,k}=\frac{c_{i,j,k}-c_{i-1,j,k}}{c_{i+1,j,k}-c_{i,j,k}},
! \quad
! r^+_{i,j,k}=\frac{c_{i+2,j,k}-c_{i+1,j,k}}{c_{i+1,j,k}-c_{i,j,k}}.
! \end{equation}
!
! It should be noted that by
! formulation (\ref{LaxWendroffForm}) this so-called P$_2$ scheme 
! is cast into the
! so-called Lax-Wendroff form, which would be recovered for
! $\tilde c_{i,j,k}^+=\tilde c_{i,j,k}^-=1$.
! 
! In order to obtain a monotonic and positive scheme, the factors
! $\tilde c_{i,j,k}^+$ are limited in the following way:
! 
! \begin{equation}\label{PDM}
! \tilde c_{i,j,k}^+ \rightarrow \max \left[
! 0,\min\left(\tilde c_{i,j,k}^+,\frac{2}{1-|C_{i,j,k}|},
! \frac{2r^+_{i,j,k}}{|C_{i,j,k}|}\right)
! \right],
! \end{equation}
! 
! and, equivalently, for $\tilde c_{i,j,k}^-$.
! This so-called PDM-limiter has been described in detail
! by \cite{LEONARD91}, who named the PDM-limited P$_2$ scheme
! also ULTIMATE QUICKEST (quadratic upstream interpolation
! for convective kinematics with estimated stream terms).
! 
! Some simpler limiters which do not exploit the third-order
! polynomial properties of the discretisation (\ref{LaxWendroffForm}) have been
! listed by \cite{ZALEZAK87}. Among those are the MUSCL scheme by
! \cite{VANLEER79},
! 
! \begin{equation}\label{MUSCL}
! \tilde c_{i,j,k}^+ \rightarrow \max \left[
! 0,\min\left(
! 2,2r^+_{i,j,k},\frac{1+r^+_{i,j,k}}{2}
! \right)
! \right],
! \end{equation}
! 
! and the Superbee scheme by \cite{ROE85},
! 
! \begin{equation}\label{Superbee}
! \tilde c_{i,j,k}^+ \rightarrow \max \left[
! 0,\min(1,2r^+_{i,j,k}),\min(r^+_{i,j,k},2)
! \right].
! \end{equation}
! 
! The selector for the schemes is {\tt method}:
!
! \vspace{0.5cm}
!
! \begin{tabular}{ll}
! {\tt method = UPSTREAM\_SPLIT}: & first-order upstream (monotone) \\
! {\tt method = P2}: & third-order polynomial (non-monotone) \\
! {\tt method = P2\_PDM}: & third-order ULTIMATE-QUICKEST (monotone) \\
! {\tt method = MUSCL}: & second-order TVD (monotone) \\
! {\tt method = Superbee}: & second-order TVD (monotone) \\
! \end{tabular}
!
! \vspace{0.5cm}
!
! Furthermore, the horizontal diffusion in $y$-direction
! with the constant diffusion
! coefficient {\tt AH} is carried out here by means of a central difference
! second-order scheme.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax
   use advection_3d, only: hi,hio,cu
   use advection_3d, only: UPSTREAM_SPLIT,P2,SUPERBEE,MUSCL,P2_PDM
   use advection_3d, only: ONE,TWO,ONE6TH
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)      :: uu(I3DFIELD),hun(I3DFIELD)
   REALTYPE, intent(in)      :: delxu(I2DFIELD),delyu(I2DFIELD)
   REALTYPE, intent(in)      :: area_inv(I2DFIELD),dt
   integer, intent(in)       :: au(E2DFIELD), az(E2DFIELD)
   REALTYPE, intent(in)      :: splitfac
   integer, intent(in)       :: method
   REALTYPE, intent(in)      :: AH
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout)   :: f(I3DFIELD)
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer         :: i,ii,j,jj,k,kk
   REALTYPE        :: c,alpha,beta,x,r,Phi,limit,fu,fc,fd
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'u_split_adv() # ',Ncall
#endif

   cu = _ZERO_

! Calculating u-interface fluxes !
   select case (method)
      case (UPSTREAM_SPLIT)
         do k=1,kmax
            do j=jmin,jmax
               do i=imin-1,imax
                  cu(i,j,k) = _ZERO_
                  if (au(i,j) .gt. 0) then
                     if (uu(i,j,k) .gt. _ZERO_) then
                        cu(i,j,k)=uu(i,j,k)*f(i,j,k)
                     else
                        cu(i,j,k)=uu(i,j,k)*f(i+1,j,k)
                     end if
                  end if
               end do
            end do
         end do
      case ((P2),(Superbee),(MUSCL),(P2_PDM))
         do k=1,kmax
            do j=jmin,jmax
               do i=imin-1,imax
                  cu(i,j,k) = _ZERO_
                  if (au(i,j) .gt. 0) then
                     if (uu(i,j,k) .gt. _ZERO_) then
                        c=uu(i,j,k)/hun(i,j,k)*dt/delxu(i,j)
                        if (au(i-1,j) .gt. 0) then
                           fu=f(i-1,j,k)         ! upstream
                        else
                           fu=f(i  ,j,k)
                        end if
                        fc=f(i  ,j,k)            ! central
                        fd=f(i+1,j,k)            ! downstream
                        if (abs(fd-fc) .gt. 1e-10) then
                           r=(fc-fu)/(fd-fc)     ! slope ratio
                        else
                           r=(fc-fu)*1.e10
                        end if
                     else
                        c=-uu(i,j,k)/hun(i,j,k)*dt/delxu(i,j)
                        if (au(i+1,j) .gt. 0) then
                           fu=f(i+2,j,k)         ! upstream
                        else
                           fu=f(i+1,j,k)
                        end if
                        fc=f(i+1,j,k)            ! central
                        fd=f(i  ,j,k)            ! downstream
                        if (abs(fc-fd) .gt. 1e-10) then
                           r=(fu-fc)/(fc-fd)     ! slope ratio
                        else
                           r=(fu-fc)*1.e10
                        end if
                     end if
                     select case (method)
                        case ((P2),(P2_PDM))
                           x = one6th*(ONE-TWO*c)
                           Phi=(0.5+x)+(0.5-x)*r
                           if (method.eq.P2) then
                              limit=Phi
                           else
                              limit=max(_ZERO_,min(Phi,2./(1.-c), &
                                        2.*r/(c+1.e-10)))
                           end if
                        case (Superbee)
                           limit=max(_ZERO_,min(ONE,TWO*r),min(r,TWO))
                        case (MUSCL)
                           limit=max(_ZERO_,min(TWO,TWO*r,0.5*(ONE+r)))
                        case default
                           FATAL 'Not so good - do_advection_3d()'
                           stop 'u_split_adv'
                     end select
                     cu(i,j,k)=uu(i,j,k)*(fc+0.5*limit*(1-c)*(fd-fc))
!Horizontal diffusion
                     if ( AH.gt.0. .and. az(i,j).gt.0 .and. az(i+1,j).gt.0 ) then
                        cu(i,j,k) = cu(i,j,k)-AH*hun(i,j,k) &
                                     *(f(i+1,j,k)-f(i,j,k))/delxu(i,j)
                     end if
                  end if
               end do
            end do
         end do
   end select

   do k=1,kmax   ! Doing the u-advection step
      do j=jmin,jmax
         do i=imin,imax
            if (az(i,j) .eq. 1) then
               hio(i,j,k)=hi(i,j,k)
               hi(i,j,k)=hio(i,j,k)                           &
                        -splitfac*dt*(uu(i,j,k)*delyu(i,j)    &
                                     -uu(i-1,j,k)*delyu(i-1,j))*area_inv(i,j)
               f(i,j,k)=(f(i,j,k)*hio(i,j,k)           &
                        -splitfac*dt*(cu(i,j,k)*delyu(i,j)    &
                        -cu(i-1,j,k)*delyu(i-1,j))*area_inv(i,j) &
                        )/hi(i,j,k)
            end if
         end do
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving u_split_adv()'
   write(debug,*)
#endif
   return
   end subroutine u_split_adv

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
