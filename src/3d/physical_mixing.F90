#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: physical_mixing()
!
! !INTERFACE:
   subroutine physical_mixing(F,diffusivity,pm3d,pm2d,AH_method)
!
! !DESCRIPTION:
!
! Here, the physical tracer variance decay for the tracer $F$,
! ${ D}^{\mbox{\scriptsize phys}}\left(\langle F \rangle^2 \right)$, 
! due to horizontal and vertical
! mixing is calculated as proposed in \cite{BURCHARDea08b}:
! \begin{equation}
! { D}^{\mbox{\scriptsize phys}}\left(F^2 \right)=
!  2 K_h  \left(\partial_x  F \right)^2
! +2 K_h  \left(\partial_y  F \right)^2
! +2 K_v  \left(\partial_z  F \right)^2.
! \end{equation}
!
!
! !USES:
   use domain,       only: imin,imax,jmin,jmax,kmax,az
   use variables_3d, only: dt,nuh,hn

   IMPLICIT NONE

! !INPUT PARAMETERS:
   REALTYPE, intent(in)  :: F(I3DFIELD)
   REALTYPE, intent(in)  :: diffusivity
   integer, intent(in)   :: AH_method
!
! !OUTPUT PARAMETERS
   REALTYPE, intent(out) :: pm3d(I3DFIELD)
   REALTYPE, intent(out) :: pm2d(I2DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Hannes Rennau
!
! !LOCAL VARIABLES:
   REALTYPE                  :: dupper,dlower
   integer                   :: i,j,k

!EOP
!-----------------------------------------------------------------------
!BOC

   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1) then
            if (AH_method .eq. 0) then
               pm3d(i,j,:) = _ZERO_
            end if
            pm2d(i,j)=_ZERO_
            dlower=_ZERO_
            do k=1,kmax
               if (k .eq. kmax) then
                  dupper=_ZERO_
               else
                  dupper=_TWO_*(diffusivity+nuh(i,j,k))*(F(i,j,k+1)-F(i,j,k))**2 &
                         /(_HALF_*(hn(i,j,k+1)+hn(i,j,k)))**2
               end if
               pm3d(i,j,k) = pm3d(i,j,k) + _HALF_*(dupper+dlower)
               dlower=dupper
               pm2d(i,j)=pm2d(i,j)+pm3d(i,j,k)*hn(i,j,k)
            end do
         end if
      end do
   end do

   return
   end subroutine physical_mixing
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2010 - Hannes Rennau, Richard Hofmeister               !
!-----------------------------------------------------------------------
