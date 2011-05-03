!$Id: variables_les.F90,v 1.10 2009-09-30 11:28:44 bjb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: variables_les - global variables for les
!
! !INTERFACE:
   module variables_les
!
! !DESCRIPTION:
!  The module contains public subroutines to initialise and cleanup.
!  Depending whether the compiler option {STATIC} is set or not,
!  memory for les variables is statically or dynamically allocated, see
!  {\tt PUBLIC DATA MEMBERS}.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax

   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
#ifdef STATIC
#include "static_les.h"
#else
#include "dynamic_declarations_les.h"
#endif
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
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
! !IROUTINE: init_variables_les - initialise les related stuff.
!
! !INTERFACE:
   subroutine init_variables_les(runtype,Am_method,Am_const)

   IMPLICIT NONE
!
! !DESCRIPTION:
!  Allocates memory (unless {\tt STATIC} is set) for les related fields,
!  by an include statement.
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: runtype
   integer, intent(in)                 :: Am_method
   REALTYPE, intent(in)                :: Am_const
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !LOCAL VARIABLES:
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_variables_les() # ',Ncall
#endif

   if (Am_method .eq. 2) LEVEL2 'init_variables_les()'

!  Allocates memory for the public data members - if not static
#ifndef STATIC
#include "dynamic_allocations_les.h"
#endif

   if (Am_method .eq. 1) then
      Am2d  = Am_const
      AmX2d = Am_const
#ifndef NO_3D
      if (runtype .ne. 1) then
         Am3d  = Am_const
         AmX3d = Am_const
         AmU3d = Am_const
         AmV3d = Am_const
      end if
#endif
   else
      Am2d  = _ZERO_
      AmX2d = _ZERO_
#ifndef NO_3D
      if (runtype .ne. 1) then
         Am3d  = _ZERO_
         AmX3d = _ZERO_
         AmU3d = _ZERO_
         AmV3d = _ZERO_
      end if
#endif
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving init_variables_les()'
   write(debug,*)
#endif
   return
   end subroutine init_variables_les
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clean_variables_les - cleanup after run.
!
! !INTERFACE:
   subroutine clean_variables_les()
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This routine is currently empty.
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'clean_variables_les() # ',Ncall
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving clean_variables_les()'
   write(debug,*)
#endif
   return
   end subroutine clean_variables_les
!EOC

!-----------------------------------------------------------------------

   end module variables_les

!-----------------------------------------------------------------------
! Copyright (C) 2011 - Knut Klingbeil                                  !
!-----------------------------------------------------------------------
