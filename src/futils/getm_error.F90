!$Id: getm_error.F90,v 1.1 2003-10-30 16:29:37 kbk Exp $
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
!  Revision 1.1  2003-10-30 16:29:37  kbk
!  added global error handler
!
!
! !LOCAL VARIABLES:
   integer                   :: ierr
!
!EOP
!-----------------------------------------------------------------------
!BOC
   include "mpif.h"

   FATAL "Called from: ",trim(sub)
   FATAL "Message:     ",trim(msg)

   call MPI_Abort(MPI_COMM_WORLD,1,ierr)

   return
   end subroutine getm_error
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2003 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
