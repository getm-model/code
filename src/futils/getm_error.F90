#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: getm_error() - global error reporting routine
!
! !INTERFACE:
   subroutine getm_error(sub,msg)
!
! !DESCRIPTION:
!
! !USES:
#ifdef GETM_PARALLEL
    use mpi
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*)                    :: sub,msg
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: ierr
!EOP
!-----------------------------------------------------------------------
!BOC
   FATAL "Called from: ",trim(sub)
   FATAL "Message:     ",trim(msg)
#ifdef GETM_PARALLEL
   call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#else
   stop "getm_error()"
#endif

   return
   end subroutine getm_error
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2003 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
