!$Id: create_restart_ncdf.F90,v 1.7 2009-09-25 12:14:56 kb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Create a GETM NetCDFNetCDF  hotstart file
!
! !INTERFACE:
   subroutine create_restart_ncdf(fname,loop,runtype)
!
! !DESCRIPTION:
!  Creates a new NetCDF formatted file for storing variables necessary
!  to make a correct GETM hotstart. The created file contains dimensions
!  (xax, yax, zax) as well as the (empty) variables. Variables are named
!  corresponding to the names used in the Fortran files. Only the actual
!  domain is stored (i.e. not the halo-zones). This allows easy use of 
!  'ncmerge' to stitch a number of hotstart files together to cover the
!  entire computational domain. See read\_restart\_ncdf() for use.
!
! !USES:
   use netcdf
   use ncdf_restart
   use domain, only: ioff,joff
   use domain, only: imin,imax,jmin,jmax,kmax
   use domain, only: vert_cord
#ifdef GETM_BIO
   use bio_var, only: numc
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fname
   integer, intent(in)                 :: loop
   integer, intent(in)                 :: runtype
!
! !DEFINED PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  $Log: create_restart_ncdf.F90,v $
!  Revision 1.7  2009-09-25 12:14:56  kb
!  INCLUDE_HALOS --> SAVE_HALOS
!
!  Revision 1.6  2009-09-23 09:54:52  kb
!  fixed typos in DESCRIPTION
!
!  Revision 1.5  2009-08-21 10:39:00  kb
!  -DINCLUDE_HALOS will include halo-zones when writing/reading NetCDF hotstart files
!
!  Revision 1.4  2009-07-18 12:36:01  kb
!  fixed SPM hot-start bug - Hofmeister
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
   character(len=80)         :: history,tts
   character(len=80)         :: title
   character(len=80)         :: str_error 
!EOP
!-------------------------------------------------------------------------
!BOC
!  create netCDF file
   status = nf90_create(fname, NF90_CLOBBER, ncid)
   if (status .NE. NF90_NOERR) go to 10

!  length of netCDF dimensions
#ifdef SAVE_HALOS
   xlen = (imax+HALO)-(imin-HALO)+1
   ylen = (jmax+HALO)-(jmin-HALO)+1
#else
   xlen = imax-imin+1
   ylen = jmax-jmin+1
#endif
   zlen = kmax+1

   status = nf90_def_dim(ncid, "xax", xlen, xdim_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_dim(ncid, "yax", ylen, ydim_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_dim(ncid, "zax", zlen, zdim_id)
   if (status .NE. NF90_NOERR) go to 10

#ifdef GETM_BIO
   status = nf90_def_dim(ncid, "biodim", numc, biodim_id)
   if (status .NE. NF90_NOERR) go to 10
#endif

   status = nf90_def_var(ncid, "loop", nf90_int, loop_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "julianday", nf90_int, julianday_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "secondsofday", nf90_int, secondsofday_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "timestep", nf90_double, timestep_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "xax", nf90_double, (/ xdim_id /), xax_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "yax", nf90_double, (/ ydim_id /), yax_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "zax", nf90_double, (/ zdim_id /), zax_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "z", nf90_double, &
                            (/ xdim_id, ydim_id /), z_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "zo", nf90_double, &
                            (/ xdim_id, ydim_id /), zo_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "U", nf90_double, &
                            (/ xdim_id, ydim_id /), U_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "zu", nf90_double, &
                            (/ xdim_id, ydim_id /), zu_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "SlUx", nf90_double, &
                            (/ xdim_id, ydim_id /), SlUx_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "Slru", nf90_double, &
                            (/ xdim_id, ydim_id /), Slru_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "V", nf90_double, &
                            (/ xdim_id, ydim_id /), V_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "zv", nf90_double, &
                            (/ xdim_id, ydim_id /), zv_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "SlVx", nf90_double, &
                            (/ xdim_id, ydim_id /), SlVx_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "Slrv", nf90_double, &
                            (/ xdim_id, ydim_id /), Slrv_id)
   if (status .NE. NF90_NOERR) go to 10

#ifndef NO_3D
   if (runtype .ge. 2)  then
      status = nf90_def_var(ncid, "ssen", nf90_double, &
                               (/ xdim_id, ydim_id /), ssen_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "ssun", nf90_double, &
                               (/ xdim_id, ydim_id /), ssun_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "ssvn", nf90_double, &
                               (/ xdim_id, ydim_id /), ssvn_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "sseo", nf90_double, &
                               (/ xdim_id, ydim_id /), sseo_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "ssuo", nf90_double, &
                               (/ xdim_id, ydim_id /), ssuo_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "ssvo", nf90_double, &
                               (/ xdim_id, ydim_id /), ssvo_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "Uinto", nf90_double, &
                               (/ xdim_id, ydim_id /), Uinto_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "Vinto", nf90_double, &
                               (/ xdim_id, ydim_id /), Vinto_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "uu", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), uu_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "vv", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), vv_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "ww", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), ww_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "uuEx", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), uuEx_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "vvEx", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), vvEx_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "tke", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), tke_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "eps", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), eps_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "num", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), num_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "nuh", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), nuh_id)
      if (status .NE. NF90_NOERR) go to 10

#ifndef NO_BAROCLINIC
      status = nf90_def_var(ncid, "T", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), T_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "S", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), S_id)
      if (status .NE. NF90_NOERR) go to 10
#endif
#ifdef SPM
      status = nf90_def_var(ncid, "spm", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), spm_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "spmpool", nf90_double, &
                               (/ xdim_id, ydim_id /), spmpool_id)
      if (status .NE. NF90_NOERR) go to 10
#endif
#ifdef GETM_BIO
      status = nf90_def_var(ncid, "bio", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id, biodim_id /), &
                                bio_id)
      if (status .NE. NF90_NOERR) go to 10
#endif
   end if
#endif

!  globals
   title="GETM NetCDF hotstart file"
   status = nf90_put_att(ncid,NF90_GLOBAL,'title',trim(title))
   if (status .NE. NF90_NOERR) go to 10

   history = 'Generated by GETM, ver. '//RELEASE
   status = nf90_put_att(ncid,NF90_GLOBAL,'history',trim(history))
   if (status .NE. NF90_NOERR) go to 10

   ! leave define mode
   status = nf90_enddef(ncid)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_sync(ncid)
   if (status .NE. NF90_NOERR) go to 10

   return

   10 FATAL 'create_restart_ncdf: ',nf90_strerror(status)
   stop 'create_restart_ncdf'

   return
   end subroutine create_restart_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Karsten Bolding (BBH)                           !
!-----------------------------------------------------------------------
