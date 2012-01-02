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
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax

   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   REALTYPE,dimension(:,:),allocatable    :: AmC_2d,AmX_2d
#ifndef NO_3D
   REALTYPE,dimension(:,:,:),allocatable  :: AmC_3d,AmX_3d,AmU_3d,AmV_3d
#endif
   integer,parameter :: NO_LES=0,LES_MOMENTUM=1,LES_TRACER=2,LES_BOTH=3
   integer           :: les_mode=NO_LES
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
   integer :: rc
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
   subroutine init_variables_les(runtype)
!
! !DESCRIPTION:
!  Allocates memory (unless {\tt STATIC} is set) for les related fields,
!  by an include statement.

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in) :: runtype
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_variables_les() # ',Ncall
#endif

   LEVEL2 'init_variables_les'

   if (les_mode.eq.LES_MOMENTUM .or. les_mode.eq.LES_BOTH) then
      allocate(AmC_2d(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (AmC_2d)'
      AmC_2d  = _ZERO_

      allocate(AmX_2d(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (AmX_2d)'
      AmX_2d = _ZERO_
   end if
#ifndef NO_3D
   if (runtype .ne. 1) then
      if (les_mode.eq.LES_MOMENTUM .or. les_mode.eq.LES_BOTH) then
         allocate(AmC_3d(I3DFIELD),stat=rc)
         if (rc /= 0) stop 'init_2d: Error allocating memory (AmC_3d)'
         AmC_3d  = _ZERO_

         allocate(AmX_3d(I3DFIELD),stat=rc)
         if (rc /= 0) stop 'init_2d: Error allocating memory (AmX_3d)'
         AmX_3d = _ZERO_
      end if
      if (les_mode.eq.LES_TRACER .or. les_mode.eq.LES_BOTH) then
         allocate(AmU_3d(I3DFIELD),stat=rc)
         if (rc /= 0) stop 'init_2d: Error allocating memory (AmU_3d)'
         AmU_3d = _ZERO_
#ifndef SLICE_MODEL
         allocate(AmV_3d(I3DFIELD),stat=rc)
         if (rc /= 0) stop 'init_2d: Error allocating memory (AmV_3d)'
         AmV_3d = _ZERO_
#endif
      end if
   end if
#endif

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
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
