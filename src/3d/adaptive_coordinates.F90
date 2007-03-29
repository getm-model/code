!$Id: adaptive_coordinates.F90,v 1.1 2007-03-29 12:28:22 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:  adaptive vertical coordinates
! \label{sec-adaptive-coordinates}
!
! !INTERFACE:
   subroutine adaptive_coordinates(first)
!
! !DESCRIPTION:
!  For Richard to do
!
! !USES:
#if 0
   use domain, only: iimin,iimax,jjmin,jjmax,kmax,H,HU,HV,az,au,av,min_depth
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
!  Original author(s): Richard Hofmeister and Hans Burchard
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

STDERR 'adaptive_coordinates()'

   if (first) then
!  do the initialisation 
   end if ! first

! and here do the updates

#ifdef DEBUG
   write(debug,*) 'Leaving adaptive_coordinates()'
   write(debug,*)
#endif
   return
   end subroutine adaptive_coordinates
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
