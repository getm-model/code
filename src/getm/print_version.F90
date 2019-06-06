#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: print_version() - prints GETM, GOTM and optional FABM versions
!
! !INTERFACE:
   subroutine print_version()
!
! !DESCRIPTION:
!  Use git to obtain latest git hashes for GETM, GOTM and FABM.
!  Also print compiler information from GOTM.
!
! !USES:
   use getm_version, only: getm_commit_id=>git_commit_id, &
                           getm_branch_name=>git_branch_name
   use gotm_version, only: gotm_commit_id=>git_commit_id, &
                           gotm_branch_name=>git_branch_name
   use gotm_compilation
#ifdef _FABM_
   use fabm_version, only: fabm_commit_id=>git_commit_id, &
                           fabm_branch_name=>git_branch_name
#endif
   use netcdf
   IMPLICIT NONE
   integer :: v, sv
   character(len=8) :: ver
!
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL0 LINE
   LEVEL0 'GETM:   ',getm_commit_id,' (',getm_branch_name,' branch)'
   LEVEL0 'GOTM:   ',gotm_commit_id,' (',gotm_branch_name,' branch)'
#ifdef _FABM_
   LEVEL0 'FABM:   ',fabm_commit_id,' (',fabm_branch_name,' branch)'
#endif
   LEVEL0 'NetCDF: ',trim(NF90_INQ_LIBVERS())
#ifdef GETM_PARALLEL
   call MPI_get_version(v,sv)
   write(ver,'(i1,a1,i1)') v,'.',sv
   LEVEL0 'MPI:    ',trim(ver)
#endif
   LEVEL0 LINE
   LEVEL0 'Compiler: ',compiler_id,' ',compiler_version
   return
   end subroutine print_version
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2016 - Karsten Bolding (BB)                            !
!-----------------------------------------------------------------------

