!$Id: getm_error.F90,v 1.4 2009-08-21 08:56:34 bjb Exp $
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
!  $Log: getm_error.F90,v $
!  Revision 1.4  2009-08-21 08:56:34  bjb
!  Fix name clash on PARALLEL with OpenMP key word
!
!  Revision 1.3  2004-04-06 16:54:33  kbk
!  cleaned a little
!
!  Revision 1.2  2003/11/03 09:10:41  kbk
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
