!$Id: ncdf_common.F90,v 1.4 2009-03-13 14:44:14 kb Exp $
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
!  $Log: ncdf_common.F90,v $
!  Revision 1.4  2009-03-13 14:44:14  kb
!  grid information in NF_DOUBLE
!
!  Revision 1.3  2005-04-25 09:32:34  kbk
!  added NetCDF IO rewrite + de-stag of velocities - Umlauf
!
!  Revision 1.2  2003/04/23 11:54:03  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.1.1.1  2002/05/02 14:01:49  gotm
!  recovering after CVS crash
!
!  Revision 1.1  2001/09/13 14:50:02  bbh
!  Cleaner and smaller NetCDF implementation + better axis support
!
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
