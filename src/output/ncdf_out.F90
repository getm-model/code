#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdf_out - interfaces for NetCDF IO subroutines
!
! !INTERFACE:
   module ncdf_out
!
! !DESCRIPTION:
!  This file contains the interface for the files in the ../ncdf/ directory
!  that does the actual saving in the NetCDF format. To avoid circular
!  dependencies it has to be here instead of - more naturally - in the
!  ncdf-directory.
!
! !USE:
   IMPLICIT NONE
!
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_2d_ncdf,save_2d_ncdf
   public init_3d_ncdf,save_3d_ncdf
   public init_mean_ncdf,save_mean_ncdf
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-----------------------------------------------------------------------

   interface
      subroutine init_2d_ncdf(fn,title,starttime)
         character(len=*), intent(in)  :: fn,title,starttime
      end subroutine init_2d_ncdf
   end interface

   interface
      subroutine save_2d_ncdf(secs)
         REALTYPE, intent(in)          :: secs
      end subroutine save_2d_ncdf
   end interface

   interface
      subroutine init_3d_ncdf(fn,title,starttime)
         character(len=*), intent(in)  :: fn,title,starttime
      end subroutine init_3d_ncdf
   end interface

   interface
      subroutine save_3d_ncdf(secs)
         REALTYPE, intent(in)          :: secs
      end subroutine save_3d_ncdf
   end interface

   interface
      subroutine init_mean_ncdf(fn,title,starttime)
         character(len=*), intent(in)  :: fn,title,starttime
      end subroutine init_mean_ncdf
   end interface

   interface
      subroutine save_mean_ncdf(secs)
         REALTYPE, intent(in)          :: secs
      end subroutine save_mean_ncdf
   end interface

!-----------------------------------------------------------------------

   end module ncdf_out

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Karsten Bolding and Hans Burchard (BBH)         !
!-----------------------------------------------------------------------
