!$Id: gotm.F90,v 1.13 2006-03-01 15:54:08 kbk Exp $
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
! There are furthermore a number of compiler options provided, e.g.\
! for an older GOTM version, for barotropic calcuations,
! and for simple parabolic viscosity profiles circumventing the GOTM
! turbulence module. 
!
! !USES:
   use halo_zones, only: update_3d_halo,wait_halo,H_TAG
   use domain, only: iimin,iimax,jjmin,jjmax,kmax,az,min_depth,crit_depth
   use variables_2d, only: D,zub,zvb,z
   use variables_3d, only: dt,kmin,ho,hn,tke,eps,SS,num,taus,taub
#ifndef NO_BAROCLINIC
   use variables_3d, only: NN,nuh
#endif
   use variables_3d, only: avmback,avhback
   use turbulence, only: do_turbulence,cde
   use turbulence, only: tke1d => tke, eps1d => eps, L1d => L
   use turbulence, only: num1d => num, nuh1d => nuh
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   REALTYPE                  :: u_taus,u_taub,z0s,z0b
   REALTYPE                  :: h(0:kmax),dry,zz
#ifdef GOTM_3_0
   REALTYPE                  :: P(0:kmax),B(0:kmax)
#endif
   REALTYPE                  :: NN1d(0:kmax),SS1d(0:kmax)
#ifndef GOTM_3_0
   REALTYPE                  :: xP(0:kmax)
#endif
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'gotm() # ',Ncall
#endif

#ifndef GOTM_3_0
   xP = _ZERO_
#endif
   do j=jjmin,jjmax
      do i=iimin,iimax

         if (az(i,j) .eq. 1 ) then

            u_taus = sqrt(taus(i,j))
            u_taub = sqrt(taub(i,j))

            h = hn(i,j,:)
            SS1d = SS(i,j,:)
#ifndef NO_BAROCLINIC
            NN1d = NN(i,j,:)
#endif

            tke1d=tke(i,j,:)
            eps1d=eps(i,j,:)
            L1d=cde*tke1d**1.5/eps1d
            num1d=num(i,j,:)
#ifndef NO_BAROCLINIC
            nuh1d=nuh(i,j,:)
#endif

            z0s = 0.1
            z0b=0.5*(max(zub(i-1,j),zub(i,j))+max(zvb(i,j-1),zvb(i,j)))
            if (z0s .gt. D(i,j)/10.) z0s=D(i,j)/10.

#ifdef PARABOLIC_VISCOSITY
            zz = _ZERO_
            do k=1,kmax-1
               zz=zz+hn(i,j,k)
               tke1d(k)=max(1.e-10,3.333333*u_taub**2*(1.-zz/D(i,j)))
               L1d(k)=0.4*(zz+z0b)*sqrt(1.-zz/D(i,j))
               eps1d(k)=0.16431677*tke1d(k)**1.5/L1d(k)
               num1d(k)=0.09*tke1d(k)**2/eps1d(k)
#ifndef NO_BAROCLINIC
               nuh1d(k)=num1d(k)
#endif
            end do
#else
#ifdef GOTM_3_0

            P  = num(i,j,:)*SS1d
#ifndef NO_BAROCLINIC
            B  = -nuh(i,j,:)*NN1d
#endif

            call do_turbulence(kmax,dt,D(i,j),u_taus,u_taub,z0s,z0b,h, &
                               NN1d,SS1d,P,B)
#else
            call do_turbulence(kmax,dt,D(i,j),u_taus,u_taub,z0s,z0b,h, &
                               NN1d,SS1d,xP)
#endif
#endif

            tke(i,j,:) = tke1d
            eps(i,j,:) = eps1d
            num(i,j,:) = num1d + avmback
#ifndef NO_BAROCLINIC
            nuh(i,j,:) = nuh1d + avhback
#endif
         end if
      end do
   end do
   
   call update_3d_halo(num,num,az,iimin,jjmin,iimax,jjmax,kmax,H_TAG)
   call wait_halo(H_TAG)
#ifndef NO_BAROCLINIC
   call update_3d_halo(nuh,nuh,az,iimin,jjmin,iimax,jjmax,kmax,H_TAG)
   call wait_halo(H_TAG)
#endif


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
