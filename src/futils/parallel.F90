#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: parallel - mpi interface to 'getm'
!
! !INTERFACE:
   module kurt_parallel
!
! !DESCRIPTION:
!
! !USES:
#ifdef GETM_PARALLEL
   use halo_mpi, only: postinit_mpi,print_MPI_info,barrier,myid
#endif
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
#ifndef GETM_PARALLEL
   integer, parameter                  :: myid=-1
#endif
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
! !IROUTINE: init_parallel - initialize MPI environment
!
! !INTERFACE:
   subroutine init_parallel(runid,input_dir)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*)                    :: runid,input_dir
!
! !DESCRIPTION:
!  Initialize Parallel environment
!
! !REVISION HISTORY:
!
!  22Apr99   Karsten Bolding & Hans Burchard  Initial code.
!
! !LOCAL VARIABLES:
#ifdef GETM_PARALLEL
   logical                   :: TO_FILE=.true.
   character(len=3)          :: buf
   character(len=16)         :: pid,ext
   character(len=PATH_MAX)   :: fname
#endif
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_parallel'
#endif

#ifdef GETM_PARALLEL
   call postinit_mpi(input_dir)
   call print_MPI_info()
   if (TO_FILE) then
      if (myid .ge. 0) then
         write(buf,'(I3.3)') myid
         pid = '.' // TRIM(buf)
      else
         pid = ''
      end if
      ext   = 'stderr'
      fname = TRIM(runid) // TRIM(pid) // '.' // ext
      open(stderr,file=Fname)

      ext   = 'stdout'
      fname = TRIM(runid) // TRIM(pid) // '.' // ext
      open(stdout,file=Fname)
   end if
   call print_MPI_info()
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving init_parallel()'
   write(debug,*)
#endif
   return
   end subroutine init_parallel
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clean_parallel - close down the parallel environment
!
! !INTERFACE:
   subroutine clean_parallel()
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Initialize Parallel environment
!
! !REVISION HISTORY:
!
!  22Apr99   Karsten Bolding & Hans Burchard  Initial code.
!
! !LOCAL VARIABLES:
#ifdef GETM_PARALLEL
   integer              :: ierr
#endif
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'clean_parallel'
#endif

#ifdef GETM_PARALLEL
   LEVEL2 'At final MPI barrier'
   call barrier()
   LEVEL2 'About to finish parallel part of GETM - calling MPI_Finalize()'  
   if(myid .ge. 0) then
      call MPI_Finalize(ierr)
   end if
   LEVEL2 'MPI finalized'
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving clean_parallel()'
   write(debug,*)
#endif
   return
   end subroutine clean_parallel
!EOC

!-----------------------------------------------------------------------

   end module kurt_parallel

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
