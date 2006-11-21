!$Id: parallel.F90,v 1.3 2006-11-21 15:10:47 frv-bjb Exp $
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
#ifdef PARALLEL
   use halo_mpi, only: postinit_mpi,print_MPI_info,myid
#endif
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
#ifndef PARALLEL
   integer, parameter                  :: myid=-1
#endif
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: parallel.F90,v $
!  Revision 1.3  2006-11-21 15:10:47  frv-bjb
!  Parallel independence of INPUT_DIR for getm.inp read. Unset INPUT_DIR to use MPI working dir.
!
!  Revision 1.2  2003-04-23 12:02:43  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.1  2003/04/07 12:05:42  kbk
!  new parallel related files
!
!  Revision 1.1.1.1  2002/05/02 14:01:29  gotm
!  recovering after CVS crash
!
!  Revision 1.4  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.3  2001/05/18 12:53:08  bbh
!  Prepared for mask in update_2d_halo - but not used yet
!
!  Revision 1.2  2001/05/18 10:03:44  bbh
!  Added mask in parameter list to update_3d_halo()
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
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
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Initialize Parallel environment
!
! !REVISION HISTORY:
!
!  22Apr99   Karsten Bolding & Hans Burchard  Initial code.
!
! !LOCAL VARIABLES:
#ifdef PARALLEL
   logical                   :: TO_FILE=.true.
   character(len=3)          :: buf
   character(len=16)         :: pid,ext
   character(len=PATH_MAX)   :: fname
#endif
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_parallel'
#endif

#ifdef PARALLEL
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
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Initialize Parallel environment
!
! !REVISION HISTORY:
!
!  22Apr99   Karsten Bolding & Hans Burchard  Initial code.
!
! !LOCAL VARIABLES:
#ifdef PARALLEL
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

#ifdef PARALLEL
   if(myid .ge. 0) then
      call MPI_Finalize(ierr)
   end if
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
