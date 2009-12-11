!$Id: set_attributes.F90,v 1.4 2009-12-11 08:41:15 kb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Sets various attributes for a NetCDF variable.
!
! !INTERFACE:
   subroutine set_attributes(ncid,id,                            &
                             units,long_name,                    &
                             netcdf_real,                        &
                             valid_min,valid_max,valid_range,    &
                             scale_factor,add_offset,            &
                             FillValue,missing_value,            &
                             C_format,FORTRAN_format)
!
! !DESCRIPTION:
!  This routine is used to set a number of attributes for the various
!  variables. The routine make heavy use of the {\em optional} keyword.
!  The list of recognized keywords is very easy expandable. We have
!  included a sub-set of the COARDS conventions.
!
! !USES:
!  IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: ncid,id
   integer, optional                   :: netcdf_real
   character(len=*), optional          :: units,long_name
   REALTYPE, optional                  :: valid_min,valid_max,valid_range(2)
   REALTYPE, optional                  :: scale_factor,add_offset
   REALTYPE, optional                  :: FillValue,missing_value
   character(len=*), optional          :: C_format,FORTRAN_format
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See ncdfout module
!
! !LOCAL VARIABLES:
   integer                   :: len,iret
   integer                   :: ft
   REAL_4B                   :: vals(2)
!
!EOP
!-----------------------------------------------------------------------
!BOC
   include "netcdf.inc"

   if(present(netcdf_real)) then
      ft=netcdf_real
   else
      ft=NF_FLOAT
   end if

   if(present(units)) then
      len = len_trim(units)
      iret = nf_put_att_text(ncid,id,'units',len,units)
   end if

   if(present(long_name)) then
      len = len_trim(long_name)
      iret = nf_put_att_text(ncid,id,'long_name',len,long_name)
   end if

   if(present(C_format)) then
      len = len_trim(C_format)
      iret = nf_put_att_text(ncid,id,'C_format',len,C_format)
   end if

   if(present(FORTRAN_format)) then
      len = len_trim(FORTRAN_format)
      iret = nf_put_att_text(ncid,id,'FORTRAN_format',len,FORTRAN_format)
   end if

   if(present(valid_min)) then
      vals(1) = valid_min
      iret = nf_put_att_real(ncid,id,'valid_min',ft,1,vals)
   end if

   if(present(valid_max)) then
      vals(1) = valid_max
      iret = nf_put_att_real(ncid,id,'valid_max',ft,1,vals)
   end if

   if(present(valid_range)) then
      vals(1) = valid_range(1)
      vals(2) = valid_range(2)
      iret = nf_put_att_real(ncid,id,'valid_range',ft,2,vals)
   end if

   if(present(scale_factor)) then
      vals(1) = scale_factor
      iret = nf_put_att_real(ncid,id,'scale_factor',ft,1,vals)
   end if

   if(present(add_offset)) then
      vals(1) = add_offset
      iret = nf_put_att_real(ncid,id,'add_offset',ft,1,vals)
   end if

   if(present(FillValue)) then
      vals(1) = FillValue
      iret = nf_put_att_real(ncid,id,'_FillValue',ft,1,vals)
   end if

   if(present(missing_value)) then
      vals(1) = missing_value
      iret = nf_put_att_real(ncid,id,'missing_value',ft,1,vals)
   end if

   return
   end subroutine set_attributes
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Karsten Bolding and Hans Burchard (BBH)         !
!-----------------------------------------------------------------------
