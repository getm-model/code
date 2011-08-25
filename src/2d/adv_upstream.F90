#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:  adv_upstream - 2D upstream advection \label{sec-upstream-adv}
!
! !INTERFACE:
   subroutine adv_upstream(dt,f,U,V,Do,Dn, &
                           delxv,delyu,delxu,delyv,area_inv,az,AH)
! !DESCRIPTION:
!
! Here the advection terms are calculated by a upstream scheme and
! the advection is done in one single time step.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
   use advection, only: cu
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                        :: dt,AH
   REALTYPE,dimension(E2DFIELD),intent(in)    :: U,V,Do,Dn
   REALTYPE,dimension(E2DFIELD),intent(in)    :: delxv,delyu,delxu,delyv
   REALTYPE,dimension(E2DFIELD),intent(in)    :: area_inv
   integer,dimension(E2DFIELD),intent(in)     :: az
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(inout) :: f
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer         :: rc,i,ii,j,jj
#ifdef USE_ALLOCATED_ARRAYS
   REALTYPE, dimension(:,:), allocatable       :: adv
#else
   REALTYPE        :: adv(I2DFIELD)
#endif
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'adv_upstream() # ',Ncall
#endif

#ifdef USE_ALLOCATED_ARRAYS
   allocate(adv(I2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'adv_upstream: Error allocating memory (adv)'
#endif

! Note: We do not need to initialize adv.
!   Tested BJB 2009-09-25.


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)

!  Calculating u-interface fluxes !
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin-1,imax
         if (U(i,j) .gt. _ZERO_) then
            cu(i,j)=U(i,j)*f(i,j)
         else
            cu(i,j)=U(i,j)*f(i+1,j)
         end if
         if ((AH.gt._ZERO_).and.(az(i,j).gt.0).and.(az(i+1,j).gt.0))    &
            cu(i,j)=cu(i,j)-AH*(f(i+1,j)-f(i,j))/delxu(i,j) &
                      *_HALF_*(Dn(i+1,j)+Dn(i,j))
      end do
   end do
!$OMP END DO NOWAIT

!$OMP BARRIER

!  Updating the advection term for u-advection !
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         adv(i,j)=( cu(i  ,j)*delyu(i  ,j)    &
                   -cu(i-1,j)*delyu(i-1,j))*area_inv(i,j)
      end do
   end do
!$OMP END DO NOWAIT

!$OMP BARRIER

#ifndef SLICE_MODEL
!  Calculating v-interface fluxes !
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-1,jmax
      do i=imin,imax
         if (V(i,j) .gt. _ZERO_) then
            cu(i,j)=V(i,j)*f(i,j)
         else
            cu(i,j)=V(i,j)*f(i,j+1)
         end if
         if ((AH.gt._ZERO_).and.(az(i,j).gt.0).and.(az(i,j+1).gt.0))   &
            cu(i,j)=cu(i,j)-AH*(f(i,j+1)-f(i,j))/delyv(i,j)   &
                      *_HALF_*(Dn(i,j+1)+Dn(i,j))
      end do
   end do
!$OMP END DO NOWAIT

!$OMP BARRIER

!  Updating the advection term for v-advection !
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         adv(i,j)=adv(i,j)+( cu(i,j  )*delxv(i,j  )   &
                            -cu(i,j-1)*delxv(i,j-1))*area_inv(i,j)
      end do
   end do
!$OMP END DO NOWAIT
#endif

!$OMP BARRIER

!  Doing the full advection in one step
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1)                                        &
            f(i,j)=(f(i,j)*Do(i,j)-dt*adv(i,j))/Dn(i,j)
      end do
   end do
!$OMP END DO

!$OMP END PARALLEL

#ifdef USE_ALLOCATED_ARRAYS
#ifdef FORTRAN90
   deallocate(adv,stat=rc)    ! work array
   if (rc /= 0) stop 'adv_upstream: Error de-allocating memory (adv)'
#endif
#endif

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
