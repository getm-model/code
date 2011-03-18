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
   REALTYPE, parameter                 :: g = 9.81
   REALTYPE, parameter                 :: rho_0 = 1025.
   REALTYPE, parameter                 :: cp = 3985.
   REALTYPE, parameter                 :: kappa = 0.4
   REALTYPE, parameter                 :: avmmol = 1.8e-6
   REALTYPE, parameter                 :: avmolt=  1.4e-7
   REALTYPE, parameter                 :: avmols = 1.1e-9
!
!  Turbulence related constants - see www.gotm.net
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-----------------------------------------------------------------------

   end module parameters

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
