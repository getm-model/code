!$Id: variables_2d.F90,v 1.4 2003-04-23 12:09:44 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: variables_2d - global variables for 2D model
!
! !INTERFACE:
   module variables_2d
!
! !DESCRIPTION:
!  This modules contains declarations for all variables related to 2D
!  hydrodynamical calculations. Information about the calculation domain
!  is included from the \emph{domain.F90} module.
!  The module contains public subroutines to initialise and cleanup.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,H,HU,HV,min_depth
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
#ifdef STATIC
#include "static_2d.h"
#else
#include "dynamic_declarations_2d.h"
#endif
   integer                             :: size2d_field
   integer                             :: mem2d
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: variables_2d.F90,v $
!  Revision 1.4  2003-04-23 12:09:44  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.3  2003/04/07 15:50:20  kbk
!  initialise variables
!
!  Revision 1.1.1.1  2002/05/02 14:00:47  gotm
!  recovering after CVS crash
!
!  Revision 1.1  2001/05/03 19:30:41  bbh
!  2D variables seperated from m2d
!
! !LOCAL VARIABLES:
   integer                   :: rc
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_variables_2d - initialise 2D relatedstuff.
!
! !INTERFACE:
   subroutine init_variables_2d(runtype)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Allocates memory for 2D related fields.
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: runtype
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!
!  See log for module.
!
! !LOCAL VARIABLES:
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_variables_2d() # ',Ncall
#endif

   LEVEL2 'init_variables_2d'
   size2d_field=((imax+HALO)-(imin+HALO)+1)*((jmax+HALO)-(jmin+HALO)+1)
   mem2d=n2d_fields*size2d_field*REAL_SIZE

!  Allocates memory for the public data members - if not static
#ifndef STATIC
#include "dynamic_allocations_2d.h"
#endif

   z = _ZERO_; zu = _ZERO_; zv = _ZERO_
   U = _ZERO_; DU = _ZERO_; fU = _ZERO_; SlUx = _ZERO_; Slru = _ZERO_
   V = _ZERO_; DV = _ZERO_; fV = _ZERO_; SlVx = _ZERO_; Slrv = _ZERO_

   Uint = _ZERO_; Vint = _ZERO_
   UEx = _ZERO_; VEx = _ZERO_
   ru = _ZERO_; rv = _ZERO_
   res_du = _ZERO_; res_u = _ZERO_; res_dv = _ZERO_; res_v =  _ZERO_
   surfdiv = _ZERO_

#ifdef DEBUG
   write(debug,*) 'Leaving init_variables_2d()'
   write(debug,*)
#endif
   return
   end subroutine init_variables_2d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clean_variables_2d - cleanup after 2D run.
!
! !INTERFACE:
   subroutine clean_variables_2d()
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This routine cleans up after a 2D integration. Close open files etc.
!
! !REVISION HISTORY:
!  See log for module.
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'clean_variables_2d() # ',Ncall
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving clean_variables_2d()'
   write(debug,*)
#endif
   return
   end subroutine clean_variables_2d
!EOC

!-----------------------------------------------------------------------

   end module variables_2d

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
