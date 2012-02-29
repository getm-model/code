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
!$ use omp_lib
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
! !LOCAL VARIABLES:
   integer                           :: i,j,k,rc
   REALTYPE                          :: zi,cff,FC
   REALTYPE,dimension(:),allocatable :: hw,delrho,dZ,dR
   REALTYPE,dimension(I2DFIELD)      :: delzx,delrhox,dZx,dRx
   REALTYPE,dimension(I3DFIELD)      :: P
   REALTYPE,parameter                :: OneFifth=0.2d0,OneTwelfth=_ONE_/12.0d0
!EOP
!-----------------------------------------------------------------------
!BOC

!  KK-TODO: put this in a central place (if needed at all)
   idpdx(:,:,0) = _ZERO_
#ifndef SLICE_MODEL
   idpdy(:,:,0) = _ZERO_
#endif

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          PRIVATE(i,j,k,rc)                                       &
!$OMP          PRIVATE(zi,cff,FC)                                      &
!$OMP          PRIVATE(hw,delrho,dZ,dR)

! OMP-NOTE: Each thread allocates its own HEAP storage for the
!    vertical work storage:
   allocate(hw(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'ip_shchepetkin_mcwilliams(): Error allocating memory (hw)'
   allocate(delrho(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'ip_shchepetkin_mcwilliams(): Error allocating memory (delrho)'
   allocate(dZ(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'ip_shchepetkin_mcwilliams(): Error allocating memory (dZ)'
   allocate(dR(0:kmax),stat=rc)    ! work array
   if (rc /= 0) stop 'ip_shchepetkin_mcwilliams(): Error allocating memory (dR)'

!$OMP DO SCHEDULE(RUNTIME)
!  First, the rho-point heights are calculated
   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO
         if (az(i,j) .ge. 1) then
            zi = -H(i,j)
            do k=1,kmax
               zz(i,j,k) = zi + _HALF_*hn(i,j,k)
               zi = zi + hn(i,j,k)
               hw(k-1) = _HALF_ * ( hn(i,j,k-1) + hn(i,j,k) )
               delrho(k-1) = buoy(i,j,k) - buoy(i,j,k-1)
            end do
            hw(kmax) = _HALF_ * hn(i,j,kmax) ! only used for convenience

            do k=2,kmax-1

               dZ(k) = _TWO_*hw(k-1)*hw(k) / ( hw(k-1) + hw(k) )

               cff = _TWO_*delrho(k-1)*delrho(k)
               if (cff .gt. _ZERO_) then
                  dR(k) = cff / ( delrho(k-1) + delrho(k) )
               else
                  dR(k) = _ZERO_
               end if

            end do

!           Note (KK): these are the BCs implemented in the ROMS code
            dZ(   1) =     hw(     1)
            dZ(kmax) =     hw(kmax-1)
            dR(   1) = delrho(     1)
            dR(kmax) = delrho(kmax-1)

!           Note (KK): these are the BCs suggested in the original paper
!            dZ(   1) = _HALF_ * ( _THREE_*    hw(     1) - dZ(     2) )
!            dZ(kmax) = _HALF_ * ( _THREE_*    hw(kmax-1) - dZ(kmax-1) )
!            dR(   1) = _HALF_ * ( _THREE_*delrho(     1) - dR(     2) )
!            dR(kmax) = _HALF_ * ( _THREE_*delrho(kmax-1) - dR(kmax-1) )

            P(i,j,kmax) = (buoy(i,j,kmax)+_HALF_*hw(kmax)*delrho(kmax-1)/hw(kmax-1))*hw(kmax)
            do k=kmax-1,1,-1
               P(i,j,k) =   P(i,j,k+1)                                                                 &
                          + _HALF_*(                                                                   &
                                     (buoy(i,j,k)+buoy(i,j,k+1))*hw(k)                                 &
                                    -OneFifth*(                                                        &
                                                (dR(k+1)-dR(k))*(    hw(k)-OneTwelfth*(dZ(k)+dZ(k+1))) &
                                               -(dZ(k+1)-dZ(k))*(delrho(k)-OneTwelfth*(dR(k)+dR(k+1))) &
                                              )                                                        &
                                   )
            end do
         end if
      end do
   end do
!$OMP END DO


!  Compute pressure gradient term as it
!  appears on the right hand side of u equation

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO
         if (au(i,j) .eq. 0) then
            delzx(i,j) = _ZERO_
            delrhox(i,j) = _ZERO_
         end if
      end do
   end do
!$OMP END DO

   do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO-1
            if (au(i,j) .ge. 1) then
               delzx(i,j) = zz(i+1,j,k) - zz(i,j,k)
               delrhox(i,j) = buoy(i+1,j,k) - buoy(i,j,k)
            end if
         end do
         do i=imin-HALO+1,imax+HALO-1
            if (az(i,j) .ge. 1) then
               cff = _TWO_*delzx(i-1,j)*delzx(i,j)
               if (cff .gt. _ZERO_) then
                  dZx(i,j) = cff / ( delzx(i-1,j) + delzx(i,j) )
               else
                  dZx(i,j) = _ZERO_
               end if

               cff = _TWO_*delrhox(i-1,j)*delrhox(i,j)
               if (cff .gt. _ZERO_) then
                  dRx(i,j) = cff / ( delrhox(i-1,j) + delrhox(i,j) )
               else
                  dRx(i,j) = _ZERO_
               end if
            end if
         end do
         do i=imin-HALO+1,imax+HALO-2
            if (au(i,j) .ge. 1) then
               FC = _HALF_*(                                                                                  &
                             (buoy(i,j,k)+buoy(i+1,j,k))*delzx(i,j)                                           &
                            -OneFifth*(                                                                       &
                                        (dRx(i+1,j)-dRx(i,j))*(  delzx(i,j)-OneTwelfth*(dZx(i,j)+dZx(i+1,j))) &
                                       -(dZx(i+1,j)-dZx(i,j))*(delrhox(i,j)-OneTwelfth*(dRx(i,j)+dRx(i+1,j))) &
                                      )                                                                       &
                           )
               idpdx(i,j,k)=_HALF_*(hn(i,j,k)+hn(i+1,j,k))*(P(i+1,j,k)-P(i,j,k)+FC)/DXU
            end if
         end do
      end do
!$OMP END DO
   end do


#ifndef SLICE_MODEL

!  Compute pressure gradient term as it
!  appears on the right hand side of v equation

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO-1
      do i=imin-HALO,imax+HALO
         if (av(i,j) .eq. 0) then
            delzx(i,j) = _ZERO_
            delrhox(i,j) = _ZERO_
         end if
      end do
   end do
!$OMP END DO

   do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO
            if (av(i,j) .ge. 1) then
               delzx(i,j) = zz(i,j+1,k) - zz(i,j,k)
               delrhox(i,j) = buoy(i,j+1,k) - buoy(i,j,k)
            end if
         end do
      end do
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO+1,jmax+HALO-1
         do i=imin-HALO,imax+HALO
            if (az(i,j) .ge. 1) then
               cff = _TWO_*delzx(i,j-1)*delzx(i,j)
               if (cff .gt. _ZERO_) then
                  dZx(i,j) = cff / ( delzx(i,j-1) + delzx(i,j) )
               else
                  dZx(i,j) = _ZERO_
               end if

               cff = _TWO_*delrhox(i,j-1)*delrhox(i,j)
               if (cff .gt. _ZERO_) then
                  dRx(i,j) = cff / ( delrhox(i,j-1) + delrhox(i,j) )
               else
                  dRx(i,j) = _ZERO_
               end if
            end if
         end do
      end do
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO+1,jmax+HALO-2
         do i=imin-HALO,imax+HALO
            if (av(i,j) .ge. 1) then
               FC = _HALF_*(                                                                                  &
                             (buoy(i,j,k)+buoy(i,j+1,k))*delzx(i,j)                                           &
                            -OneFifth*(                                                                       &
                                        (dRx(i,j+1)-dRx(i,j))*(  delzx(i,j)-OneTwelfth*(dZx(i,j)+dZx(i,j+1))) &
                                       -(dZx(i,j+1)-dZx(i,j))*(delrhox(i,j)-OneTwelfth*(dRx(i,j)+dRx(i,j+1))) &
                                      )                                                                       &
                           )
               idpdy(i,j,k)=_HALF_*(hn(i,j,k)+hn(i,j+1,k))*(P(i,j+1,k)-P(i,j,k)+FC)/DYV
            end if
         end do
      end do
!$OMP END DO
   end do

#endif


!$OMP END PARALLEL

   return
   end subroutine ip_shchepetkin_mcwilliams
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2007 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
