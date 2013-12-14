#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: flux_center2interface -
!
! !INTERFACE:
   subroutine flux_center2interface(tagc,fluxc,ttag,fluxi)
!
! !USES:
   use halo_zones, only: U_TAG,V_TAG
   use domain    , only: imin,imax,jmin,jmax,au,av
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain         , only: dxc,dyc,dxu,dyu,dxv,dyv
#else
   use domain         , only: dx,dy
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)                       :: tagc,ttag
   REALTYPE,dimension(E2DFIELD),intent(in)  :: fluxc
!
! !OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(out) :: fluxi
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES
   integer :: i,j
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'flux_center2interface() # ',Ncall
#endif

   select case(tagc)

      case (U_TAG)

         select case(ttag)

            case (U_TAG)
               do j=jmin-HALO,jmax+HALO
                  do i=imin-HALO,imax+HALO-1
                     if (au(i,j).eq.1 .or. au(i,j).eq.2) then
                        fluxi(i,j) = _HALF_*( DYC*fluxc(i,j) + DYCIP1*fluxc(i+1,j) ) / DYU
                     end if
                  end do
               end do
            case (V_TAG)
               do j=jmin-HALO,jmax+HALO-1
                  do i=imin-HALO,imax+HALO
                     if (av(i,j).eq.1 .or. av(i,j).eq.2) then
                        fluxi(i,j) = _HALF_*( DYC*fluxc(i,j) + DYCJP1*fluxc(i,j+1) ) / DYV
                     end if
                  end do
               end do

         end select

      case (V_TAG)

         select case(ttag)

            case (U_TAG)
               do j=jmin-HALO,jmax+HALO
                  do i=imin-HALO,imax+HALO-1
                     if (au(i,j).eq.1 .or. au(i,j).eq.2) then
                        fluxi(i,j) = _HALF_*( DXC*fluxc(i,j) + DXCIP1*fluxc(i+1,j) ) / DXU
                     end if
                  end do
               end do
            case (V_TAG)
               do j=jmin-HALO,jmax+HALO-1
                  do i=imin-HALO,imax+HALO
                     if (av(i,j).eq.1 .or. av(i,j).eq.2) then
                        fluxi(i,j) = _HALF_*( DXC*fluxc(i,j) + DXCJP1*fluxc(i,j+1) ) / DXV
                     end if
                  end do
               end do

         end select

   end select

#ifdef DEBUG
   write(debug,*) 'Leaving flux_center2interface()'
   write(debug,*)
#endif
   return
   end subroutine flux_center2interface
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2013 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
