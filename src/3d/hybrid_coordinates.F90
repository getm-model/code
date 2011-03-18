#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:  hybrid vertical coordinates - (z-sigma)
! \label{sec-hybrid-coordinates}
!
! !INTERFACE:
   subroutine hybrid_coordinates(first)
!
! !DESCRIPTION:
! For Hans to do
!
! !USES:
#if 0
   use domain, only: imin,imax,jmin,jmax,kmax,H,HU,HV,az,au,av,min_depth
   use domain, only: ga,ddu,ddl,d_gamma,gamma_surf
   use variables_3d, only: dt,kmin,kumin,kvmin,ho,hn,huo,hun,hvo,hvn
   use variables_3d, only: sseo,ssen,ssuo,ssun,ssvo,ssvn
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical, intent(in)                 :: first
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
#if 0
   integer         :: i,j,k,rc,kk
   REALTYPE, save, dimension(:),     allocatable  :: dga,be,sig,zlev
   REALTYPE, save, dimension(:,:,:), allocatable  :: gga
#endif
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'coordinates() # ',Ncall
#endif

STDERR 'hybrid_coordinates()'

   if (first) then
!  do the initialisation 
   end if ! first

! and here do the updates

#ifdef DEBUG
   write(debug,*) 'Leaving hybrid_coordinates()'
   write(debug,*)
#endif
   return
   end subroutine hybrid_coordinates
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
