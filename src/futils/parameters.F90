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
   REALTYPE :: g = 9.81d0
   REALTYPE :: rho_0 = 1025.0d0
   REALTYPE :: cp = 3985.0d0
   REALTYPE :: kappa = 0.4d0
!
!  Turbulence related constants - see www.gotm.net
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_parameters - reads in parameters from namelist
!
! !INTERFACE:
   subroutine init_parameters()
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Reads the namelist.
!
! !LOCAL VARIABLES:
   namelist /parameters/ &
             g,rho_0,cp,kappa,avmmol,avmolt,avmols
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_parameters() # ',Ncall
#endif

   LEVEL1 'init_parameters'

   read(NAMLST,parameters)

   if (g .lt. _ZERO_) then
      LEVEL2 'changed sign of negative g'
      g = -g
   end if
   LEVEL2 'g = ',real(g)

   if (rho_0 .le. _ZERO_) then
      stop 'init_parameters(): invalid rho_0'
   else
      LEVEL2 'rho_0 = ',real(rho_0)
   end if

   LEVEL2 'cp = ',real(cp)
   LEVEL2 'kappa = ',real(kappa)

#ifdef DEBUG
   write(debug,*) 'Leaving init_parameters()'
   write(debug,*)
#endif
   return
   end subroutine init_parameters
!EOC
!-----------------------------------------------------------------------

   end module parameters

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
