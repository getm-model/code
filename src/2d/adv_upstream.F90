#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:  adv_upstream - 2D upstream advection \label{sec-upstream-adv}
!
! !INTERFACE:
   subroutine adv_upstream(dt,f,Di,adv,U,V,Do,Dn, &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                           dxv,dyu,dxu,dyv,arcd1, &
#endif
                           az,AH)
!
! !DESCRIPTION:
!
! Here the advection terms are calculated by a upstream scheme and
! the advection is done in one single time step.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
#if !( defined(SPHERICAL) || defined(CURVILINEAR) )
   use domain, only: dx,dy,ard1
#endif
   use advection, only: flux
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                        :: dt,AH
   REALTYPE,dimension(E2DFIELD),intent(in)    :: U,V,Do,Dn
#if defined(SPHERICAL) || defined(CURVILINEAR)
   REALTYPE,dimension(E2DFIELD),intent(in)    :: dxv,dyu,dxu,dyv,arcd1
#endif
   integer,dimension(E2DFIELD),intent(in)     :: az
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(inout) :: f,Di,adv
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer :: i,j
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'adv_upstream() # ',Ncall
#endif

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)

!  Calculating u-interface fluxes !
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin-1,imax
         if (U(i,j) .gt. _ZERO_) then
            flux(i,j) = U(i,j)*f(i,j)
         else
            flux(i,j) = U(i,j)*f(i+1,j)
         end if
         if ((AH.gt._ZERO_).and.(az(i,j).gt.0).and.(az(i+1,j).gt.0)) then
            flux(i,j) = flux(i,j) - AH*( f(i+1,j)                                 &
                                        -f(i  ,j))/DXU*_HALF_*(Dn(i+1,j)+Dn(i,j))
         end if
      end do
   end do
!$OMP END DO

!  Updating the advection term for u-advection !
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         adv(i,j) = (flux(i  ,j)*DYU - flux(i-1,j)*DYUIM1)*ARCD1
      end do
   end do
!$OMP END DO

#ifndef SLICE_MODEL
!  Calculating v-interface fluxes !
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-1,jmax
      do i=imin,imax
         if (V(i,j) .gt. _ZERO_) then
            flux(i,j) = V(i,j)*f(i,j)
         else
            flux(i,j) = V(i,j)*f(i,j+1)
         end if
         if ((AH.gt._ZERO_).and.(az(i,j).gt.0).and.(az(i,j+1).gt.0)) then
            flux(i,j) = flux(i,j) - AH*( f(i,j+1)                                 &
                                        -f(i,j  ))/DYV*_HALF_*(Dn(i,j+1)+Dn(i,j))
         end if
      end do
   end do
!$OMP END DO

!  Updating the advection term for v-advection !
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         adv(i,j) = adv(i,j) + (flux(i,j  )*DXV - flux(i,j-1)*DXVJM1)*ARCD1
      end do
   end do
!$OMP END DO
#endif

!  Doing the full advection in one step
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1) then
            f(i,j)=(f(i,j)*Do(i,j)-dt*adv(i,j))/Dn(i,j)
         end if
      end do
   end do
!$OMP END DO

!$OMP END PARALLEL

#ifdef DEBUG
   write(debug,*) 'Leaving adv_upstream()'
   write(debug,*)
#endif
   return
   end subroutine adv_upstream
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
