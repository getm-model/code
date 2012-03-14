#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: physical_mixing()
!
! !INTERFACE:
   subroutine physical_mixing(F,AH,diffusivity,pm3d,pm2d)
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
   use domain,       only: imin,imax,jmin,jmax,kmax,H,au,av
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxu,dyv
#else
   use domain, only: dx,dy
#endif
   use variables_3d, only: dt,nuh,hn,ssen

   IMPLICIT NONE

! !INPUT PARAMETERS:
   REALTYPE, intent(in)  :: F(I3DFIELD)
   REALTYPE, intent(in)  :: AH
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
   REALTYPE                  :: aux(I3DFIELD)

!EOP
!-----------------------------------------------------------------------
!BOC
   pm3d=_ZERO_
   if (AH .gt. _ZERO_) then
      ! Physical dissipation on U-points
      do k=1,kmax
         do j=jmin,jmax
            do i=imin-1,imax
               if (au(i,j) .gt. 0) then
                  aux(i,j,k)=_TWO_*AH*((F(i+1,j,k)-F(i,j,k))/DXU)**2  
               else
                  aux(i,j,k)=_ZERO_
               end if
            end do
         end do
      end do
      do k=1,kmax
         do j=jmin,jmax
            do i=imin,imax
               pm3d(i,j,k)=pm3d(i,j,k)+_HALF_*(aux(i,j,k)+aux(i-1,j,k))  
            end do
         end do
      end do
      ! Physical dissipation on V-points
      do k=1,kmax
         do j=jmin-1,jmax
            do i=imin,imax
               if (av(i,j) .gt. 0) then
                  aux(i,j,k)=_TWO_*AH*((F(i,j+1,k)-F(i,j,k))/DYV)**2  
               else
                  aux(i,j,k)=_ZERO_
               end if
            end do
         end do
      end do
      do k=1,kmax
         do j=jmin,jmax
            do i=imin,imax
               pm3d(i,j,k)=pm3d(i,j,k)+_HALF_*(aux(i,j,k)+aux(i,j-1,k))  
            end do
         end do
      end do
   end if
   do j=jmin,jmax
      do i=imin,imax
         pm2d(i,j)=_ZERO_
         dlower=_ZERO_
         do k=1,kmax
           if (k .eq. kmax) then
             dupper=_ZERO_
           else
             dupper=_TWO_*(diffusivity+nuh(i,j,k))*(F(i,j,k+1)-F(i,j,k))**2 &
                      /(_HALF_*(hn(i,j,k+1)+hn(i,j,k)))**2
           end if
           pm3d(i,j,k)=pm3d(i,j,k)+_HALF_*(dupper+dlower)
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
