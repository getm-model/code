!$Id: structure_friction_3d.F90,v 1.1 2008-03-26 13:25:51 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: structure_friction_3d - 
! \label{sec-structure_friction-3d}
!
! !INTERFACE:
   subroutine structure_friction_3d()
!
! !DESCRIPTION:
! Here, the quadratic friction term resulting from a structure in the water column is
! calculated. This term will be added as additional forcing to the
! three-dimensional momentum equations, where it is treated numerically
! implicitly. Therefore here, only the following terms is calculated:
! \begin{equation}
! \mbox{\tt sf} = C(z) \sqrt{u(z)^2+v(z)^2},
! \end{equation}
! with the friction coefficient $C$ bearing the physical unit [1/m]. 
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax
   use variables_3d, only: uu,vv,sf,huo,hvo
!  #define CALC_HALO_WW
#ifndef CALC_HALO_WW
   use domain, only: az
   use halo_zones, only: update_3d_halo,wait_halo,z_TAG
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                  :: dtm1
   integer                   :: i,j,k
#ifdef STRUCTURE_FRICTION
   REALTYPE                  :: cds(I2DFIELD)
#endif

!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'structure_friction_3d() # ',Ncall
#endif

!  cds is the friction coefficient which a structure imposes to a
!  flow field. It has the unit [1/m] and should be a 3d field. For the
!  first, a dependency on the vertical coordinate is not assumed.
!  cds should be later read from a netcdf input file.
   cds=0
   do j=jmin,jmax
      do i=imin,imax
         if ((i.eq.90).and.(j.eq.6)) then
            cds(i,j)=0.01
         end if   
      end do
   end do
   sf=_ZERO_
   do k=1,kmax
#ifdef CALC_HALO_WW
      do j=jmin-1,jmax+1
         do i=imin-1,imax+1
#else
      do j=jmin,jmax
         do i=imin,imax
#endif
               sf(i,j,k)=cds(i,j)*sqrt(0.5*(                                &
                        (uu(i  ,j  ,k)/(huo(i  ,j  ,k)))**2                 &
                       +(uu(i-1,j  ,k)/(huo(i-1,j  ,k)))**2                 &
                       +(vv(i  ,j  ,k)/(hvo(i  ,j  ,k)))**2                 &
                       +(vv(i  ,j-1,k)/(hvo(i  ,j-1,k)))**2                 &
                       ))
         end do
      end do
   end do

#ifndef CALC_HALO_WW
   call update_3d_halo(sf,sf,az,imin,jmin,imax,jmax,kmax,z_TAG)
   call wait_halo(z_TAG)
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving structure_friction_3d()'
   write(debug,*)
#endif
   return
   end subroutine structure_friction_3d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
