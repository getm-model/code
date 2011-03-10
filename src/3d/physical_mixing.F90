#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: physical_mixing() 
!
! !INTERFACE:
   subroutine physical_mixing(F,diffusivity,pm3d,pm2d)
!
! !DESCRIPTION:
!

! !USES:
   use domain,       only: imin,imax,jmin,jmax,kmax,H
   use variables_3d, only: dt,nuh,hn,ssen
   IMPLICIT NONE

! !INPUT PARAMETERS:
   REALTYPE, intent(in)  :: F(I3DFIELD)
   REALTYPE, intent(in)  :: diffusivity
!
! !INPUT PARAMETERS
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
         pm2d(i,j)=_ZERO_
         dlower=_ZERO_
         do k=1,kmax-1
           if (k.eq.kmax) then
             dupper=_ZERO_
           else
             dupper=_TWO_*(diffusivity+nuh(i,j,k))*(F(i,j,k+1)-F(i,j,k))**2 &
                      /(_HALF_*(hn(i,j,k+1)+hn(i,j,k)))**2
           end if
           pm3d(i,j,k)=_HALF_*(dupper+dlower)
           dlower=dupper
           pm2d(i,j)=pm2d(i,j)+pm3d(i,j,k)*hn(i,j,k)
         end do
      end do
   end do 
   return
   end subroutine physical_mixing
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2010 - Hannes Rennau, Richard Hofmeister               !
!-----------------------------------------------------------------------

