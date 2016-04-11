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
   IMPLICIT NONE
!
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL0 LINE
   LEVEL0 'GETM version: ',getm_commit_id,' (',getm_branch_name,' branch)'
   LEVEL0 'GOTM version: ',gotm_commit_id,' (',gotm_branch_name,' branch)'
#ifdef _FABM_
   LEVEL0 'FABM version: ',fabm_commit_id,' (',fabm_branch_name,' branch)'
#endif
   LEVEL0 LINE
   LEVEL0 'Compiler: ',compiler_id,' ',compiler_version
   return
   end subroutine print_version
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2016 - Karsten Bolding (BB)                            !
!-----------------------------------------------------------------------

