#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: to_fluxu() - calculates volume flux in local x-direction
!
! !INTERFACE:
   subroutine to_fluxu(imin,jmin,imax,jmax,au, &
#if defined(CURVILINEAR) || defined(SPHERICAL)
                       dyu,                    &
#else
                       dy,                     &
#endif
                       U,missing,fluxu)
!
! !DESCRIPTION:
!
! !USES:
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)                       :: imin,jmin,imax,jmax
   integer,dimension(E2DFIELD),intent(in)   :: au
#if defined(CURVILINEAR) || defined(SPHERICAL)
   REALTYPE,dimension(E2DFIELD),intent(in)  :: dyu
#else
   REALTYPE,intent(in)                      :: dy
#endif
   REALTYPE,dimension(E2DFIELD),intent(in)  :: U
   REALTYPE,intent(in)                      :: missing
!
! !OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(out) :: fluxu
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
   write(debug,*) 'to_fluxu() # ',Ncall
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
         if (au(i,j) .gt. 0) then
            fluxu(i,j) = DYU*U(i,j)
         else
            fluxu(i,j) = missing
         end if
      end do
#ifndef SLICE_MODEL
   end do
#endif
!$OMP END DO

!$OMP END PARALLEL

#ifdef SLICE_MODEL
   fluxu(:,j+1) = fluxu(:,j)
#endif

   return
   end subroutine to_fluxu
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2012 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
