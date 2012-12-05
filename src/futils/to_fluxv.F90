#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: to_fluxv() - calculates volume flux in local y-direction
!
! !INTERFACE:
   subroutine to_fluxv(imin,jmin,imax,jmax,av, &
#if defined(CURVILINEAR) || defined(SPHERICAL)
                       dxv,                    &
#else
                       dx,                     &
#endif
                       V,missing,fluxv)
!
! !DESCRIPTION:
!
! !USES:
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)                       :: imin,jmin,imax,jmax
   integer,dimension(E2DFIELD),intent(in)   :: av
#if defined(CURVILINEAR) || defined(SPHERICAL)
   REALTYPE,dimension(E2DFIELD),intent(in)  :: dxv
#else
   REALTYPE,intent(in)                      :: dx
#endif
   REALTYPE,dimension(E2DFIELD),intent(in)  :: V
   REALTYPE,intent(in)                      :: missing
!
! !OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(out) :: fluxv
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
   integer :: i,j
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'to_fluxv() # ',Ncall
#endif
#ifdef SLICE_MODEL
!  Note (KK): this value MUST NOT be changed !!!
   j = jmax/2
#endif

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j)                                         &
!$OMP          PRIVATE(i)

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
   do j=jmin-HALO,jmax+HALO
#endif
      do i=imin-HALO,imax+HALO
         if (av(i,j) .gt. 0) then
            fluxv(i,j) = DXV*V(i,j)
         else
            fluxv(i,j) = missing
         end if
      end do
#ifndef SLICE_MODEL
   end do
#endif
!$OMP END DO

!$OMP END PARALLEL

#ifdef SLICE_MODEL
   fluxv(:,j-1) = fluxv(:,j)
   fluxv(:,j+1) = fluxv(:,j)
#endif

   return
   end subroutine to_fluxv
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2012 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
