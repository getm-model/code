!$Id: open_restart_ncdf.F90,v 1.4 2009-09-23 09:54:53 kb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise restart netCDf variables
!
! !INTERFACE:
   subroutine open_restart_ncdf(fname,runtype)
!
! !DESCRIPTION:
!  Opens a NetCDF formatted GETM hotstart file. All NetCDF variable
!  id's necessary for making a correct GETM hotstart are read. The id's 
!  are shared with the reading routine using the ncdf\_restart module.
!
! !USES:
   use netcdf
   use ncdf_restart
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fname
   integer, intent(in)                 :: runtype
!
! !DEFINED PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  $Log: open_restart_ncdf.F90,v $
!  Revision 1.4  2009-09-23 09:54:53  kb
!  fixed typos in DESCRIPTION
!
!  Revision 1.3  2009-04-27 08:03:02  kb
!  getm/initialise.F90
!
!  Revision 1.2  2007-10-19 07:52:36  kbk
!  zub and zvb not in hotstart files anymore
!
!  Revision 1.1  2007-09-21 13:03:42  kbk
!  added drop-in NetCDF replacement for binary hotstart file (default is binary)
!
!
! !LOCAL VARIABLES:
   integer         :: dimids(3)
!EOP
!-------------------------------------------------------------------------
!BOC
!  open netCDF file
   status = nf90_open(fname, NF90_NOWRITE, ncid)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_inq_varid(ncid, "loop", loop_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_inq_varid(ncid, "julianday", julianday_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_inq_varid(ncid, "secondsofday", secondsofday_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_inq_varid(ncid, "timestep", timestep_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_inq_varid(ncid, "z", z_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_inquire_variable(ncid,z_id,dimids = dimids)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_inquire_dimension(ncid,dimids(1),len = xlen)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_inquire_dimension(ncid,dimids(2),len = ylen)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_inq_varid(ncid, "zo", zo_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_inq_varid(ncid, "U", U_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_inq_varid(ncid, "zu", zu_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_inq_varid(ncid, "SlUx", SlUx_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_inq_varid(ncid, "Slru", Slru_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_inq_varid(ncid, "V", V_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_inq_varid(ncid, "zv", zv_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_inq_varid(ncid, "SlVx", SlVx_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_inq_varid(ncid, "Slrv", Slrv_id)
   if (status .NE. NF90_NOERR) go to 10

#ifndef NO_3D
   if (runtype .ge. 2)  then
      status = nf90_inq_varid(ncid, "ssen", ssen_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_inq_varid(ncid, "ssun", ssun_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_inq_varid(ncid, "ssvn", ssvn_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_inq_varid(ncid, "sseo", sseo_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_inq_varid(ncid, "ssuo", ssuo_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_inq_varid(ncid, "ssvo", ssvo_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_inq_varid(ncid, "Uinto", Uinto_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_inq_varid(ncid, "Vinto", Vinto_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_inq_varid(ncid, "uu", uu_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_inq_varid(ncid, "vv", vv_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_inq_varid(ncid, "ww", ww_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_inq_varid(ncid, "uuEx", uuEx_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_inq_varid(ncid, "vvEx", vvEx_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_inq_varid(ncid, "tke", tke_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_inq_varid(ncid, "eps", eps_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_inq_varid(ncid, "num", num_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_inq_varid(ncid, "nuh", nuh_id)
      if (status .NE. NF90_NOERR) go to 10

#ifndef NO_BAROCLINIC
      status = nf90_inq_varid(ncid, "T", T_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_inq_varid(ncid, "S", S_id)
      if (status .NE. NF90_NOERR) go to 10
#endif
#ifdef SPM
      status = nf90_inq_varid(ncid, "spm", spm_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_inq_varid(ncid, "spmpool", spmpool_id)
      if (status .NE. NF90_NOERR) go to 10
#endif
#ifdef GETM_BIO
      status = nf90_inq_varid(ncid, "bio", bio_id)
      if (status .NE. NF90_NOERR) go to 10
#endif
   end if
#endif

   return

   10 FATAL 'open_restart_ncdf: ',nf90_strerror(status)
   stop 'open_restart_ncdf'

   return
   end subroutine open_restart_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Karsten Bolding (BBH)                           !
!-----------------------------------------------------------------------
