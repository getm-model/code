!$Id: ncdf_common.F90,v 1.1.1.1 2002-05-02 14:01:49 gotm Exp $
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
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public set_attributes
!
! !PUBLIC DATA MEMBERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ncdf_common.F90,v $
!  Revision 1.1.1.1  2002-05-02 14:01:49  gotm
!  recovering after CVS crash
!
!  Revision 1.1  2001/09/13 14:50:02  bbh
!  Cleaner and smaller NetCDF implementation + better axis support
!
!
!EOP
!-----------------------------------------------------------------------
   interface
      subroutine set_attributes(ncid,id,				&
                                units,long_name,			&
				valid_min,valid_max,valid_range,	&
				scale_factor,add_offset,		&
				FillValue,missing_value,		&
				C_format,FORTRAN_format)
         integer, intent(in)          :: ncid,id
         character(len=*), optional   :: units,long_name
         REALTYPE, optional           :: valid_min,valid_max,valid_range(2)
         REALTYPE, optional           :: scale_factor,add_offset
         REALTYPE, optional           :: FillValue,missing_value
         character(len=*), optional   :: C_format,FORTRAN_format
      end subroutine set_attributes
   end interface

!-----------------------------------------------------------------------

   end module ncdf_common

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Karsten Bolding and Hans Burchard (BBH)         !
!-----------------------------------------------------------------------
