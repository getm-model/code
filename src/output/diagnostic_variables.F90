!$Id: diagnostic_variables.F90,v 1.2 2006-01-29 20:32:34 hb Exp $
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

   REALTYPE,dimension(:,:,:), allocatable :: nummix3d_S_mean
   REALTYPE,dimension(:,:), allocatable :: nummix2d_S_mean
   REALTYPE,dimension(:,:,:), allocatable :: nummix3d_T_mean
   REALTYPE,dimension(:,:), allocatable :: nummix2d_T_mean
   REALTYPE,dimension(:,:,:), allocatable :: phymix3d_S_mean
   REALTYPE,dimension(:,:), allocatable :: phymix2d_S_mean
   REALTYPE,dimension(:,:,:), allocatable :: phymix3d_T_mean
   REALTYPE,dimension(:,:), allocatable :: phymix2d_T_mean

!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: diagnostic_variables.F90,v $
!  Revision 1.2  2006-01-29 20:32:34  hb
!  Small LaTeX corrections to source code documentation
!
!  Revision 1.1  2004-03-29 15:35:52  kbk
!  possible to store calculated mean fields
!
!
!EOP
!-----------------------------------------------------------------------

   end module diagnostic_variables

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
