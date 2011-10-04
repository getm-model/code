#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ip_shchepetkin_mcwilliams
!
! !INTERFACE:
   subroutine ip_shchepetkin_mcwilliams()
!
! !DESCRIPTION:
!
! Here, the pressure gradient is calculated according to the
! method and the algorithm suggested by Shchepetkin and McWilliams,
! 2003. This method uses a nonconservative  Density-Jacobian scheme,
! based on  cubic polynomial fits for the bouyancy "buoy" and "zz",
! the vertical position of rho-points, as functions of its respective
! array indices. The cubic polynomials are monotonized by using harmonic
! mean instead of linear averages to interpolate slopes.
! Exact anti-symmetry of the density Jacobian
! \begin{equation}
!        J(rho,zz)=-J(zz,rho)
! \end{equation}
! is retained for the density/bouyancy Jacobian in the pressure
! gradient formulation in x-direction for a non aligned vertical
! coordinate $\sigma$, the atmospheric pressure $p_0$ and the sea
! surface elevation $\eta$:
! \begin{equation}
! \label{eq: shchepetkin pgf}
! - {1 \over \rho_0} \partial_x p = \underbrace{-{1 \over \rho_0} \partial_x p_0 - g \partial_x\eta}_{uu\_momentum}
!                          + \underbrace{buoy(\eta) \partial_x\eta + \int_z^\eta J(buoy,zz)\mbox{d}\sigma}_{idpdx}
! \end{equation}
! Details about the calculation of the integral over the Jacobian in
! (\ref{eq: shchepetkin pgf}) can be found in Shchepetkin and McWilliams,
! 2003.
!
! If parameter OneFifth (below) is set to zero, the scheme should
! become identical to standard Jacobian.
!
! !USES:
   use internal_pressure
   use variables_3d, only: hn,buoy,sseo
   use domain, only: H,az,au,av
!$ use omp_lib
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   REALTYPE                  :: dR(I3DFIELD)
   REALTYPE                  :: dZ(I3DFIELD)
   REALTYPE                  :: P(I3DFIELD)
   REALTYPE                  :: dxm1,dym1,cff,cff1,cff2
   REALTYPE                  :: AJ
   REALTYPE                  :: eps=1.e-10
   REALTYPE                  :: OneFifth = 0.2
   REALTYPE                  :: FC(I2DFIELD)
   REALTYPE                  :: dZx(I2DFIELD)
   REALTYPE                  :: dRx(I2DFIELD)
!EOP
!-----------------------------------------------------------------------
!BOC

! OMP-NOTE: Initialization to _ZERO_ of arrays are done by master
!   threads little at a time. Typically while the remaining threads
!   execute the loop prior to the one where the arrays are needed.
!   The local work arrays do presently not need initialization.
!    BJB 2009-09-24.

   zz(:,:,0) = _ZERO_

#if ! ( defined(SPHERICAL) || defined(CURVILINEAR) )
   dxm1 = _ONE_/DXU
   dym1 = _ONE_/DYV
#endif

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP    PRIVATE(i,j,k, cff,cff1,cff2,AJ)

! BJB-NOTE: We do not need to initialize these:
!!$OMP MASTER
!   dR=_ZERO_
!   dZ=_ZERO_
!   P=_ZERO_
!!$OMP END MASTER

!$OMP DO SCHEDULE(RUNTIME)
!  First, the rho-point heights are calculated
   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO
         if (az(i,j) .ge. 1) then
            zz(i,j,1)=-H(i,j)+0.5*hn(i,j,1)
            do k=2,kmax
               zz(i,j,k)=zz(i,j,k-1)+0.5*(hn(i,j,k-1)+hn(i,j,k))
            end do
         end if
      end do
   end do
!$OMP END DO

! BJB-NOTE: dZx and dRx do not need to be initialized
!!$OMP MASTER
!   dZx=_ZERO_
!   dRx=_ZERO_
!!$OMP END MASTER

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax+1
      do i=imin,imax+1
         if (az(i,j) .ge. 1) then
            do k=1,kmax-1
               dR(i,j,k)=buoy(i,j,k+1)-buoy(i,j,k)
               dZ(i,j,k)=zz(i,j,k+1)-zz(i,j,k)
            end do
            dR(i,j,kmax)=dR(i,j,kmax-1)
            dZ(i,j,kmax)=dZ(i,j,kmax-1)
            dR(i,j,0)=dR(i,j,1)
            dZ(i,j,0)=dZ(i,j,1)

            do k=kmax,1,-1
               cff=2.0*dR(i,j,k)*dR(i,j,k-1)
               if (cff .gt. eps) then
                  dR(i,j,k)=cff/(dR(i,j,k)+dR(i,j,k-1))
               else
                  dR(i,j,k)=_ZERO_
               end if
               dZ(i,j,k)=2.0*dZ(i,j,k)*dZ(i,j,k-1)/(dZ(i,j,k)+dZ(i,j,k-1))
            end do

            cff1=1.0/(zz(i,j,kmax)-zz(i,j,kmax-1))
            cff2=0.5*(buoy(i,j,kmax)-buoy(i,j,kmax-1))*0.5*hn(i,j,kmax)*cff1

            P(i,j,kmax)=(buoy(i,j,kmax)+cff2)*0.5*hn(i,j,kmax)

            do k=kmax-1,1,-1
               P(i,j,k)=P(i,j,k+1)+0.5*((buoy(i,j,k+1)+buoy(i,j,k))* &
                    (zz(i,j,k+1)-zz(i,j,k))-OneFifth*((dR(i,j,k+1)-dR(i,j,k))* &
                    (zz(i,j,k+1)-zz(i,j,k)-1./12.*(dZ(i,j,k+1)+dZ(i,j,k)))- &
                    (dZ(i,j,k+1)-dZ(i,j,k))*(buoy(i,j,k+1)-buoy(i,j,k)- &
                    1./12.*(dR(i,j,k+1)+dR(i,j,k)))))
            end do
         end if
      end do
   end do
!$OMP END DO

! Compute pressure gradient term as it
! appears on the right hand side of u equation
! OMP-NOTE: If we want to thread the k-loop, then each thread
!    needs to allocate it's own (i,j)-size work arrays,
!    BJB 2009-09-24.
   do k=kmax,1,-1
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin-HALO+1,imax+HALO
            dZx(i,j)=zz(i,j,k)-zz(i-1,j,k)
            dRx(i,j)=buoy(i,j,k)-buoy(i-1,j,k)
         end do
      end do
!$OMP END DO

! BJB-NOTE: I am not sure that this intialization is necessary.
!    BJB 2009-09-25.
!$OMP MASTER
   idpdx(:,:,0) = _ZERO_
   idpdy(:,:,0) = _ZERO_
!$OMP END MASTER


!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin,imax+HALO-1
            cff=2.0*dZx(i,j)*dZx(i+1,j)
            if (cff .gt. eps) then
               cff1=1.0/(dZx(i,j)+dZx(i+1,j))
               dZx(i,j)=cff*cff1
            else
               dZx(i,j)=_ZERO_
            end if
            cff1=2.0*dRx(i,j)*dRx(i+1,j)
            if (cff1 .gt. eps) then
               cff2=1.0/(dRx(i,j)+dRx(i+1,j))
               dRx(i,j)=cff1*cff2
            else
               dRx(i,j)=_ZERO_
            end if
         end do
      end do
!$OMP END DO


!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin,imax+HALO-2
            if (au(i,j) .ge. 1) then
#if defined(SPHERICAL) || defined(CURVILINEAR)
               dxm1=_ONE_/DXU
#endif
               AJ = P(i+1,j,k) - P(i,j,k)
               FC(i,j) = 0.5*((buoy(i+1,j,k)+buoy(i,j,k))* &
                  (zz(i+1,j,k)-zz(i,j,k))-OneFifth*((dRx(i+1,j)-dRx(i,j))* &
                  (zz(i+1,j,k)-zz(i,j,k)-1.0/12.0*(dZx(i+1,j)+dZx(i,j)))- &
                  (dZx(i+1,j)-dZx(i,j))*(buoy(i+1,j,k)-buoy(i,j,k)-1.0/12.0*(dRx(i+1,j)+dRx(i,j)))))
               AJ = AJ + FC(i,j)
               idpdx(i,j,k)=AJ*dxm1*hun(i,j,k) - hun(i,j,k)*dxm1*(sseo(i+1,j)-sseo(i,j))*0.5*(buoy(i+1,j,kmax)+buoy(i,j,kmax))
            end if
         end do
      end do
!$OMP END DO
   end do

! Compute pressure gradient term as it
! appears on the right hand side of v equation
! OMP-NOTE: The following loops cannot be threaded along the j-dimension.
!    Threading along k seems to be non-trivial too. Thus, we reverse the
!    loop order on OMP and thread over i.
   do k=kmax,1,-1
#ifdef GETM_OMP
!$OMP DO SCHEDULE(RUNTIME)
      do i=imin,imax
         do j=jmin-HALO+1,jmax+HALO
#else
      do j=jmin-HALO+1,jmax+HALO
         do i=imin,imax
#endif

            dZx(i,j)=zz(i,j,k)-zz(i,j-1,k)
            dRx(i,j)=buoy(i,j,k)-buoy(i,j-1,k)
         end do
      end do
#ifdef GETM_OMP
!$OMP END DO
#endif

#ifdef GETM_OMP
!$OMP DO SCHEDULE(RUNTIME)
      do i=imin,imax
         do j=jmin,jmax+HALO-1
#else
      do j=jmin,jmax+HALO-1
         do i=imin,imax
#endif
            cff=2.0*dZx(i,j)*dZx(i,j+1)
            if (cff .gt. eps) then
               cff1=1.0/(dZx(i,j)+dZx(i,j+1))
               dZx(i,j)=cff*cff1
            else
               dZx(i,j)=_ZERO_
            end if
            cff1=2.0*dRx(i,j)*dRx(i,j+1)
            if (cff1 .gt. eps) then
               cff2=1.0/(dRx(i,j)+dRx(i,j+1))
               dRx(i,j)=cff1*cff2
            else
               dRx(i,j)=_ZERO_
            end if
         end do
      end do
#ifdef GETM_OMP
!$OMP END DO
#endif

#ifdef GETM_OMP
!$OMP DO SCHEDULE(RUNTIME)
      do i=imin,imax
         do j=jmin,jmax+HALO-2
#else
      do j=jmin,jmax+HALO-2
         do i=imin,imax
#endif
            if (av(i,j) .ge. 1) then
#if defined(SPHERICAL) || defined(CURVILINEAR)
               dym1 = _ONE_/DYV
#endif
               AJ = P(i,j+1,k) - P(i,j,k)
               FC(i,j) = 0.5*((buoy(i,j+1,k)+buoy(i,j,k))* &
                      (zz(i,j+1,k)-zz(i,j,k))-OneFifth*((dRx(i,j+1)-dRx(i,j))* &
                      (zz(i,j+1,k)-zz(i,j,k)- &
                      1.0/12.0*(dZx(i,j+1)+dZx(i,j)))-(dZx(i,j+1)-dZx(i,j))* &
                      (buoy(i,j+1,k)-buoy(i,j,k)-1.0/12.0*(dRx(i,j+1)+dRx(i,j)))))
               AJ = AJ + FC(i,j)
               idpdy(i,j,k)=AJ*dym1*hvn(i,j,k) - hvn(i,j,k)*dym1*(sseo(i,j+1)-sseo(i,j))*0.5*(buoy(i,j+1,kmax)+buoy(i,j,kmax))
            end if
         end do
      end do
#ifdef GETM_OMP
!$OMP END DO
#endif
   end do

!$OMP END PARALLEL

   return
   end subroutine ip_shchepetkin_mcwilliams
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Richard Hofmeister                              !
!-----------------------------------------------------------------------
