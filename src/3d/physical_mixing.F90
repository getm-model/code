#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: physical_mixing()
!
! !INTERFACE:
   subroutine physical_mixing(F,diffusivity,phymix_3d,phymix_int,AH_method)
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
   REALTYPE, intent(out) :: phymix_3d(I3DFIELD)
   REALTYPE, intent(out) :: phymix_int(I2DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Hannes Rennau
!
! !LOCAL VARIABLES:
   REALTYPE                  :: dupper,dlower,phymix_sum,phymix
   integer                   :: i,j,k

!EOP
!-----------------------------------------------------------------------
!BOC

   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1) then
            phymix_sum = _ZERO_
            dlower=_ZERO_
            do k=1,kmax
               if (k .eq. kmax) then
                  dupper=_ZERO_
               else
                  dupper=_TWO_*(diffusivity+nuh(i,j,k))* &
                         ((F(i,j,k+1)-F(i,j,k))/(_HALF_*(hn(i,j,k)+hn(i,j,k+1))))**2
               end if
               phymix = _HALF_ * ( dlower + dupper )
               if (AH_method .eq. 0) then
                  phymix_3d(i,j,k) = phymix ! avoid cumulative sum in time
               else
                  phymix_3d(i,j,k) = phymix_3d(i,j,k) + phymix
               end if
               phymix_sum = phymix_sum + phymix_3d(i,j,k)*hn(i,j,k)
               dlower=dupper
            end do
            phymix_int(i,j) = phymix_sum
         end if
      end do
   end do

   return
   end subroutine physical_mixing
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2010 - Hannes Rennau, Richard Hofmeister               !
!-----------------------------------------------------------------------
