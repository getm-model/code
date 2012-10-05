#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: gotm - a wrapper to call GOTM \label{sec-gotm}
!
! !INTERFACE:
   subroutine gotm()
!
! !DESCRIPTION:
!
! Here, the turbulence module of the General Ocean Turbulence Model (GOTM,
! see {\tt www.gotm.net} and \cite{UMLAUFea05}) is called. First, all
! necessary parameters are transformed to suit with a 1D water column model,
! i.e., 3D fields are transformed to a vertical vector, 2D horizontal
! fields are converted to a scalar. The transformed 3D fields are
! the layer heights {\tt hn $\rightarrow$ h}, the shear squared
! {\tt SS $\rightarrow$ SS1d},
! the buoyancy frequency squared {\tt NN $\rightarrow$ NN1d},
! the turbulent kinetic energy {\tt tke $\rightarrow$ tke1d},
! the dissipation rate {\tt eps $\rightarrow$ eps1d}
! (from which the integral length scale {\tt L1d} is calculated), the
! eddy viscosity {\tt num $\rightarrow$ num1d}, and the eddy diffusivity
! {\tt nuh $\rightarrow$ nuh1d}. The scalars are the surface and bottom friction
! velocities, {\tt u\_taus} and {\tt u\_taub}, respectively, the
! surface roughness parameter {\tt z0s} (which is currently hard-coded),
! and the bottom roughess parameter {\tt z0b}.
! Then, the GOTM turbulence module {\tt do\_turbulence} is called with
! all the transformed parameters discussed above. Finally, the
! vertical vectors {\tt tke1d}, {\tt eps1d}, {\tt num1d} and {\tt nuh1d}
! are transformed back to 3D fields.
!
! In case that the compiler option {\tt STRUCTURE\_FRICTION} is switched on,
! the additional turbulence production by structures in the water column is calculated
! by calculating the total production as
! \begin{equation}
! P_{tot} = P +C \left(u^2+v^2\right)^{3/2},
! \end{equation}
! with the shear production $P$, and the structure friction coefficient $C$. The
! latter is calculated in the routine {\tt structure\_friction\_3d.F90}.
!
! There are furthermore a number of compiler options provided, e.g.\
! for an older GOTM version, for barotropic calcuations,
! and for simple parabolic viscosity profiles circumventing the GOTM
! turbulence module.
!
! !USES:
   use halo_zones, only: update_3d_halo,wait_halo,H_TAG
   use domain, only: imin,imax,jmin,jmax,kmax,az,min_depth,crit_depth,z0
   use variables_2d, only: D,z
   use variables_3d, only: dt,kmin,ho,hn,tke,eps,SS,num,taus,taub,zub,zvb
#ifndef NO_BAROCLINIC
   use variables_3d, only: NN,nuh
#endif
   use variables_3d, only: avmback,avhback
#ifdef STRUCTURE_FRICTION
   use variables_3d, only: uu,vv,hun,hvn,sf
#endif
   use turbulence, only: do_turbulence,cde
   use turbulence, only: tke1d => tke, eps1d => eps, L1d => L
   use turbulence, only: num1d => num, nuh1d => nuh
   use getm_timers, only: tic, toc, TIM_GOTM, TIM_GOTMTURB, TIM_GOTMH
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   REALTYPE                  :: u_taus,u_taub,z0s,z0b
   REALTYPE                  :: h(0:kmax),dry,zz
   REALTYPE                  :: NN1d(0:kmax),SS1d(0:kmax)
   REALTYPE                  :: xP(0:kmax)
!EOP
!-----------------------------------------------------------------------
!BOC
! Note: For ifort we need to explicitly state that this routine is
! single-thread only. Presently I don't know why that is necessary,
! but if I use ifort -omp without any OMP-statements in this file,
! then the result is garbage.
! The OMP SINGLE or OMP MASTER statements helps, but sometimes it *still*
! messes up, in the sense that NaN "suddenly" appears on output.
! Apparently, writing out array-copy explicitly helps.
!    BJB 2009-09-17.
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'gotm() # ',Ncall
#endif
   call tic(TIM_GOTM)

   xP = _ZERO_
#ifdef NO_BAROCLINIC
   NN1d = _ZERO_
#endif
   do j=jmin,jmax
      do i=imin,imax

         if (az(i,j) .ge. 1 ) then

#ifdef STRUCTURE_FRICTION
! BJB-TODO: Change all constants to double
            do k=1,kmax
               xP(k)= _QUART_*(                                                &
               (uu(i  ,j  ,k)/hun(i  ,j  ,k))**2*(sf(i  ,j  ,k)+sf(i+1,j  ,k)) &
              +(uu(i-1,j  ,k)/hun(i-1,j  ,k))**2*(sf(i-1,j  ,k)+sf(i  ,j  ,k)) &
              +(vv(i  ,j  ,k)/hvn(i  ,j  ,k))**2*(sf(i  ,j  ,k)+sf(i  ,j+1,k)) &
              +(vv(i  ,j-1,k)/hvn(i  ,j-1,k))**2*(sf(i  ,j-1,k)+sf(i  ,j  ,k)))
            end do
#endif
            u_taus = sqrt(taus(i,j))
            u_taub = sqrt(taub(i,j))
            do k=0,kmax
               h(k)    = hn(i,j,k)
               SS1d(k) = SS(i,j,k)
#ifndef NO_BAROCLINIC
               NN1d(k) = NN(i,j,k)
#endif
               tke1d(k)=tke(i,j,k)
               eps1d(k)=eps(i,j,k)
               L1d(k)  =cde*tke1d(k)**1.5/eps1d(k)
               num1d(k)=num(i,j,k)
#ifndef NO_BAROCLINIC
               nuh1d(k)=nuh(i,j,k)
#endif
            end do
            z0s = _TENTH_
            z0b = _HALF_*( max( z0(i,j) , zub(i-1,j  ) , zub(i,j) ) &
                          +max( z0(i,j) , zvb(i  ,j-1) , zvb(i,j) ) )
            if (z0s .gt. D(i,j)/10.) z0s= D(i,j)/10.

#ifdef PARABOLIC_VISCOSITY
            zz = _ZERO_
            do k=1,kmax-1
               zz=zz+hn(i,j,k)
! BJB-TODO: Get rid of **1.5 and **2
               tke1d(k)=max(1.e-10,3.333333*taub(i,j)*(_ONE_-zz/D(i,j)))
               L1d(k)=0.4*(zz+z0b)*sqrt(_ONE_-zz/D(i,j))
               eps1d(k)=0.16431677*sqrt(tke1d(k)*tke1d(k)*tke1d(k))/L1d(k)
               num1d(k)=0.09*tke1d(k)*tke1d(k)/eps1d(k)
#ifndef NO_BAROCLINIC
               nuh1d(k)=num1d(k)
#endif
            end do
#else
            ! If we do tic/toc for do_turbulence, then we can
            ! easily get into the millions of system_clock calls,
            ! as the call is deeply in loops
            !call tic(TIM_GOTMTURB)
            call do_turbulence(kmax,dt,D(i,j),u_taus,u_taub,z0s,z0b,h, &
                               NN1d,SS1d,xP)
            !call toc(TIM_GOTMTURB)
#endif
            do k=0,kmax
               tke(i,j,k) = tke1d(k)
               eps(i,j,k) = eps1d(k)
               num(i,j,k) = num1d(k) + avmback
#ifndef NO_BAROCLINIC
               nuh(i,j,k) = nuh1d(k) + avhback
#endif
            end do
         end if
      end do
   end do

   call tic(TIM_GOTMH)
   call update_3d_halo(num,num,az,imin,jmin,imax,jmax,kmax,H_TAG)
   call wait_halo(H_TAG)
#ifndef NO_BAROCLINIC
   call update_3d_halo(nuh,nuh,az,imin,jmin,imax,jmax,kmax,H_TAG)
   call wait_halo(H_TAG)
#endif
   call toc(TIM_GOTMH)

   call toc(TIM_GOTM)
#ifdef DEBUG
   write(debug,*) 'Leaving gotm()'
   write(debug,*)
#endif
   return
   end subroutine gotm
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
