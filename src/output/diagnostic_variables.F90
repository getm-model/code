#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: diagnostic_variables
!
! !INTERFACE:
   module diagnostic_variables
!
! !DESCRIPTION:
!  This modules serves as a container for diagnostic variables. It is the
!  responsibillity of the subroutine(s) using these variables to properly
!  allocate memory. Have a look at {\tt .../src/output/calc\_mean\_fields.F90}.
!
! !USES:
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   REALTYPE,dimension(:,:),   allocatable :: swrmean
   REALTYPE,dimension(:,:),   allocatable :: ustarmean
   REALTYPE,dimension(:,:),   allocatable :: ustar2mean
   REALTYPE,dimension(:,:,:), allocatable :: uumean
   REALTYPE,dimension(:,:,:), allocatable :: vvmean
   REALTYPE,dimension(:,:,:), allocatable :: wmean
   REALTYPE,dimension(:,:,:), allocatable :: humean
   REALTYPE,dimension(:,:,:), allocatable :: hvmean
   REALTYPE,dimension(:,:,:), allocatable :: hmean
   REALTYPE,dimension(:,:,:), allocatable :: Tmean
   REALTYPE,dimension(:,:,:), allocatable :: Smean
   REALTYPE,dimension(:,:,:), allocatable :: rhomean
   REALTYPE,dimension(:,:,:), allocatable :: bnhmean
   REALTYPE,dimension(:,:,:), allocatable :: nummix_S_mean
   REALTYPE,dimension(:,:,:), allocatable :: nummix_S_old_mean
   REALTYPE,dimension(:,:), allocatable :: nummix_S_int_mean
   REALTYPE,dimension(:,:,:), allocatable :: nummix_T_mean
   REALTYPE,dimension(:,:,:), allocatable :: nummix_T_old_mean
   REALTYPE,dimension(:,:), allocatable :: nummix_T_int_mean
   REALTYPE,dimension(:,:,:), allocatable :: numdis_3d_mean
   REALTYPE,dimension(:,:,:), allocatable :: numdis_3d_old_mean
   REALTYPE,dimension(:,:), allocatable :: numdis_int_mean
   REALTYPE,dimension(:,:,:), allocatable :: phydis_3d_mean
   REALTYPE,dimension(:,:), allocatable :: phydis_int_mean
   REALTYPE,dimension(:,:,:), allocatable :: phymix_S_mean
   REALTYPE,dimension(:,:), allocatable :: phymix_S_int_mean
   REALTYPE,dimension(:,:,:), allocatable :: phymix_T_mean
   REALTYPE,dimension(:,:), allocatable :: phymix_T_int_mean

#ifdef GETM_BIO
   REALTYPE,dimension(:,:,:,:), allocatable :: cc3dmean
#endif
#ifdef _FABM_
   REALTYPE,dimension(:,:,:,:), allocatable :: fabmmean_pel
   REALTYPE,dimension(:,:,:), allocatable :: fabmmean_ben
   REALTYPE,dimension(:,:,:,:), allocatable :: fabmmean_diag
   REALTYPE,dimension(:,:,:), allocatable :: fabmmean_diag_hz
#endif
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-----------------------------------------------------------------------

   end module diagnostic_variables

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
