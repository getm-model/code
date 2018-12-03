#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: processed_variables
!
! !INTERFACE:
   module processed_variables
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
   REALTYPE,dimension(:,:),   allocatable :: u_2d, v_2d
   REALTYPE,dimension(:,:),   allocatable :: u_2d_destag, v_2d_destag
   REALTYPE,dimension(:,:,:), allocatable :: u_3d, v_3d
   REALTYPE,dimension(:,:,:), allocatable :: u_3d_destag, v_3d_destag
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------

   end module processed_variables

!-----------------------------------------------------------------------
! Copyright (C) 2018 - Karsten Bolding & Jorn Bruggeman (BB)           !
!-----------------------------------------------------------------------
