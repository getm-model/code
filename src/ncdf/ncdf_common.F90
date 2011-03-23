#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdf_common - interfaces for NetCDF IO subroutines
!
! !INTERFACE:
   module ncdf_common
!
! !DESCRIPTION:
!
! !USE:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-----------------------------------------------------------------------
   interface
      subroutine set_attributes(ncid,id,                               &
                                units,long_name,                       &
                                netcdf_real,                           &
                                valid_min,valid_max,valid_range,       &
                                scale_factor,add_offset,               &
                                FillValue,missing_value,               &
                                C_format,FORTRAN_format)
         integer, intent(in)           :: ncid,id
         character(len=*), optional    :: units,long_name
         integer, optional             :: netcdf_real
         REALTYPE, optional            :: valid_min,valid_max,valid_range(2)
         REALTYPE, optional            :: scale_factor,add_offset
         REALTYPE, optional            :: FillValue,missing_value
         character(len=*), optional    :: C_format,FORTRAN_format
      end subroutine set_attributes
   end interface




  interface
     subroutine init_grid_ncdf(ncid,init3d,x_dim,y_dim,z_dim)
       integer, intent(in)            :: ncid
       logical, intent(in)            :: init3d
       integer, intent(out)           :: x_dim
       integer, intent(out)           :: y_dim
       integer, intent(out), optional :: z_dim
     end subroutine init_grid_ncdf
  end interface
!
!-----------------------------------------------------------------------

   end module ncdf_common

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Karsten Bolding and Hans Burchard (BBH)         !
!-----------------------------------------------------------------------
