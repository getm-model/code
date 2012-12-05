#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: to_fluxw() - calculates vertical volume flux
!
! !INTERFACE:
   subroutine to_fluxw(imin,jmin,imax,jmax,kmin,kmax,az, &
#if defined(CURVILINEAR) || defined(SPHERICAL)
                       arcd1,                            &
#else
                       ard1,                             &
#endif
                       ww,missing,fluxw)
!
! !DESCRIPTION:
!
! !USES:
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)                       :: imin,jmin,imax,jmax,kmax
   integer,dimension(I2DFIELD),intent(in)   :: kmin
   integer,dimension(E2DFIELD),intent(in)   :: az
#if defined(CURVILINEAR) || defined(SPHERICAL)
   REALTYPE,dimension(E2DFIELD),intent(in)  :: arcd1
#else
   REALTYPE,intent(in)                      :: ard1
#endif
   REALTYPE,dimension(I3DFIELD),intent(in)  :: ww
   REALTYPE,intent(in)                      :: missing
!
! !OUTPUT PARAMETERS:
   REALTYPE,dimension(I3DFIELD),intent(out) :: fluxw
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
   integer :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'to_fluxw() # ',Ncall
#endif
#ifdef SLICE_MODEL
!  Note (KK): this value MUST NOT be changed !!!
   j = jmax/2
#endif

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j)                                         &
!$OMP          PRIVATE(i,k)

   do k=1,kmax-1

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin-HALO,jmax+HALO
#endif
         do i=imin-HALO,imax+HALO
            if (az(i,j).gt.0 .and. kmin(i,j).le.k) then
               fluxw(i,j,k) = ww(i,j,k)/ARCD1
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO

   end do

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
   do j=jmin-HALO,jmax+HALO
#endif
      do i=imin-HALO,imax+HALO
         if (az(i,j) .eq. 0) then
!           Note (KK): can be skipped if fluxw is initialised with missing
            fluxw(i,j,:) = missing
         else
            do k=0,kmin(i,j)-2
               fluxw(i,j,k) = missing
            end do
            fluxw(i,j,kmin(i,j)-1) = _ZERO_
            fluxw(i,j,kmax       ) = _ZERO_
         end if
      end do
#ifndef SLICE_MODEL
   end do
#endif
!$OMP END DO

!$OMP END PARALLEL

#ifdef SLICE_MODEL
   fluxw(:,j+1,:) = fluxw(:,j,:)
#endif

   return
   end subroutine to_fluxw
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2012 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
