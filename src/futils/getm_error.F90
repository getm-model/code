!$Id: getm_error.F90,v 1.2 2003-11-03 09:10:41 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: getm_error() - global error reporting routine
!
! !INTERFACE:
   subroutine getm_error(sub,msg)
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   character(len=*)                    :: sub,msg
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: getm_error.F90,v $
!  Revision 1.2  2003-11-03 09:10:41  kbk
!  now works with both serial and parallel compilation
!
!  Revision 1.1  2003/10/30 16:29:37  kbk
!  added global error handler
!
! !LOCAL VARIABLES:
   integer                   :: ierr
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef PARALLEL
   include "mpif.h"
#endif
   FATAL "Called from: ",trim(sub)
   FATAL "Message:     ",trim(msg)
#ifdef PARALLEL
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
