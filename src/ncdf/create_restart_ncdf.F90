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
#ifdef GETM_BIO
   use bio, only: bio_calc
   use bio_var, only: numc
#endif
#ifdef _FABM_
   use getm_fabm, only: fabm_calc,model
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fname
   integer, intent(in)                 :: loop
   integer, intent(in)                 :: runtype
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
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
#ifdef _WRITE_HOT_HALOS_
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

#ifndef NO_3D
#ifdef GETM_BIO
   if (bio_calc) then
      status = nf90_def_dim(ncid, "biodim", numc, biodim_id)
      if (status .NE. NF90_NOERR) go to 10
   end if
#endif
#ifdef _FABM_
   if (fabm_calc) then
      status = nf90_def_dim(ncid, "fabm_pel_dim", size(model%info%state_variables), fabmpeldim_id)
      if (status .NE. NF90_NOERR) go to 10

      if (size(model%info%state_variables_ben).gt.0) then
         status = nf90_def_dim(ncid, "fabm_ben_dim", size(model%info%state_variables_ben), fabmbendim_id)
         if (status .NE. NF90_NOERR) go to 10
      else
         fabmbendim_id = 0
      end if
   end if
#endif
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

   status = nf90_def_var(ncid, "SlUx", nf90_double, &
                            (/ xdim_id, ydim_id /), SlUx_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "Slru", nf90_double, &
                            (/ xdim_id, ydim_id /), Slru_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "V", nf90_double, &
                            (/ xdim_id, ydim_id /), V_id)
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

      status = nf90_def_var(ncid, "ho", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), ho_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "hn", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), hn_id)
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
      if (bio_calc) then
         status = nf90_def_var(ncid, "bio", nf90_double, &
                               (/ biodim_id, xdim_id, ydim_id, zdim_id /), &
                                bio_id)
         if (status .NE. NF90_NOERR) go to 10
      end if
#endif
#ifdef _FABM_
      if (fabm_calc) then
         status = nf90_def_var(ncid, "fabm_pel", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id, fabmpeldim_id /), &
                                fabm_pel_id)
         if (status .NE. NF90_NOERR) go to 10

         if (fabmbendim_id.gt.0) then
            status = nf90_def_var(ncid, "fabm_ben", nf90_double, &
                               (/ xdim_id, ydim_id, fabmbendim_id /), &
                                fabm_ben_id)
            if (status .NE. NF90_NOERR) go to 10
         endif
      end if
#endif
   end if
#endif

!  globals
   title="GETM NetCDF hotstart file"
   status = nf90_put_att(ncid,NF90_GLOBAL,'title',trim(title))
   if (status .NE. NF90_NOERR) go to 10

   history = 'GETM, ver. '//RELEASE
   status = nf90_put_att(ncid,NF90_GLOBAL,'history',trim(history))
   if (status .NE. NF90_NOERR) go to 10

   history = GIT_REVISION
   status = nf90_put_att(ncid,NF90_GLOBAL,'git',trim(history))
   if (status .NE. NF90_NOERR) go to 10

   history = FORTRAN_VERSION
   status = nf90_put_att(ncid,NF90_GLOBAL,'compiler',trim(history))
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
