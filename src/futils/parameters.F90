!$Id: parameters.F90,v 1.1.1.1 2002-05-02 14:01:19 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  Parameters - setup the entire model
!
! !INTERFACE:
   module parameters
!
! !DESCRIPTION:
!  Intended as a place holder for various global constants shared between
!  different parts of the model. In the future this module might be phased
!  out and the defined parameters will be moved to other - maybe - more
!  appropriate places.
!
! !USES:
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
!  Physical Constants
   REALTYPE, parameter	:: g = 9.81
   REALTYPE, parameter	:: rho_0 = 1025.
   REALTYPE, parameter	:: cp = 3985.
   REALTYPE, parameter	:: kappa = 0.4
   REALTYPE, parameter	:: avmmol = 1.8e-6
!
!  Turbulence related constants - see www.gotm.net
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: parameters.F90,v $
!  Revision 1.1.1.1  2002-05-02 14:01:19  gotm
!  recovering after CVS crash
!
!  Revision 1.6  2001/10/22 08:48:30  bbh
!  Am moved from paramters.F90 to m2d.F90
!
!  Revision 1.5  2001/10/12 09:00:23  bbh
!  AM in here - should be moved sometime
!
!  Revision 1.4  2001/07/26 13:52:02  bbh
!  Cleaning
!
!  Revision 1.3  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.2  2001/05/18 08:41:49  bbh
!  Added rho_0
!
!  Revision 1.1.1.1  2001/04/17 08:43:09  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------

   end module parameters

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
