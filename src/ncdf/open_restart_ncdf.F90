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
#ifndef NO_3D
   use vertical_coordinates,only: restart_with_ho,restart_with_hn
#ifdef GETM_BIO
   use bio, only: bio_calc
   use getm_bio, only: bio_init_method
#endif
#ifdef _FABM_
   use gotm_fabm, only: fabm_calc
   use getm_fabm, only: fabm_init_method
#endif
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fname
   integer, intent(in)                 :: runtype
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL VARIABLES:
   integer         :: dimids(3)
   character(len=20) :: varnam
!EOP
!-------------------------------------------------------------------------
!BOC
!  open netCDF file
   varnam="FILE"
   status = nf90_open(fname, NF90_NOWRITE, ncid)
   if (status .NE. NF90_NOERR) go to 10

   varnam="loop"
   status = nf90_inq_varid(ncid, "loop", loop_id)
   if (status .NE. NF90_NOERR) go to 10

   varnam="julianday"
   status = nf90_inq_varid(ncid, "julianday", julianday_id)
   if (status .NE. NF90_NOERR) go to 10

   varnam="secondsofday"
   status = nf90_inq_varid(ncid, "secondsofday", secondsofday_id)
   if (status .NE. NF90_NOERR) go to 10

   varnam="timestep"
   status = nf90_inq_varid(ncid, "timestep", timestep_id)
   if (status .NE. NF90_NOERR) go to 10

!  z is required
   varnam="z"
   status = nf90_inq_varid(ncid, "z", z_id)
   if (status .NE. NF90_NOERR) go to 10

   varnam="z dimensions"
   status = nf90_inquire_variable(ncid,z_id,dimids = dimids)
   if (status .NE. NF90_NOERR) go to 10

   varnam="z dim1"
   status = nf90_inquire_dimension(ncid,dimids(1),len = xlen)
   if (status .NE. NF90_NOERR) go to 10

   varnam="z dim2"
   status = nf90_inquire_dimension(ncid,dimids(2),len = ylen)
   if (status .NE. NF90_NOERR) go to 10

!  Some of the variables are optionally in restart file
   varnam="zo"
   status = nf90_inq_varid(ncid, "zo", zo_id)
   if (status .NE. NF90_NOERR) then
      LEVEL3 'variable missing in restart file. Skipping ',varnam
      zo_id=-1
   end if

   varnam="U"
   status = nf90_inq_varid(ncid, "U", U_id)
   if (status .NE. NF90_NOERR) then
      LEVEL3 'variable missing in restart file. Skipping ',varnam
      U_id=-1
   end if

   varnam="SlUx"
   status = nf90_inq_varid(ncid, "SlUx", SlUx_id)
   if (status .NE. NF90_NOERR) then
      LEVEL3 'variable missing in restart file. Skipping ',varnam
      SlUx_id=-1
   end if

   varnam="Slru"
   status = nf90_inq_varid(ncid, "Slru", Slru_id)
   if (status .NE. NF90_NOERR) then
      LEVEL3 'variable missing in restart file. Skipping ',varnam
      Slru_id=-1
   end if

   varnam="V"
   status = nf90_inq_varid(ncid, "V", V_id)
   if (status .NE. NF90_NOERR) then
      LEVEL3 'variable missing in restart file. Skipping ',varnam
      V_id=-1
   end if

   varnam="SlVx"
   status = nf90_inq_varid(ncid, "SlVx", SlVx_id)
   if (status .NE. NF90_NOERR) then
      LEVEL3 'variable missing in restart file. Skipping ',varnam
      SlVx_id=-1
   end if

   varnam="Slrv"
   status = nf90_inq_varid(ncid, "Slrv", Slrv_id)
   if (status .NE. NF90_NOERR) then
      LEVEL3 'variable missing in restart file. Skipping ',varnam
      Slrv_id=-1
   end if

#ifndef NO_3D
   if (runtype .ge. 2)  then
      varnam="ssen"
      status = nf90_inq_varid(ncid, "ssen", ssen_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         ssen_id=-1
      end if

      varnam="ssun"
      status = nf90_inq_varid(ncid, "ssun", ssun_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         ssun_id=-1
      end if

      varnam="ssvn"
      status = nf90_inq_varid(ncid, "ssvn", ssvn_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         ssvn_id=-1
      end if

      varnam="sseo"
      status = nf90_inq_varid(ncid, "sseo", sseo_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         sseo_id=-1
      end if

      varnam="ssuo"
      status = nf90_inq_varid(ncid, "ssuo", ssuo_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         ssuo_id=-1
      end if

      varnam="ssvo"
      status = nf90_inq_varid(ncid, "ssvo", ssvo_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         ssvo_id=-1
      end if

      varnam="Uint"
      status = nf90_inq_varid(ncid, "Uint", Uint_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         Uint_id=-1
      end if

      varnam="Vint"
      status = nf90_inq_varid(ncid, "Vint", Vint_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         Vint_id=-1
      end if

      varnam="Uinto"
      status = nf90_inq_varid(ncid, "Uinto", Uinto_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         Uinto_id=-1
      end if

      varnam="Vinto"
      status = nf90_inq_varid(ncid, "Vinto", Vinto_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         Vinto_id=-1
      end if

      varnam="Uadv"
      status = nf90_inq_varid(ncid, "Uadv", Uadv_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         Uadv_id=-1
      end if

      varnam="Vadv"
      status = nf90_inq_varid(ncid, "Vadv", Vadv_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         Vadv_id=-1
      end if

      varnam="uu"
      status = nf90_inq_varid(ncid, "uu", uu_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         uu_id=-1
      end if

      varnam="vv"
      status = nf90_inq_varid(ncid, "vv", vv_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         vv_id=-1
      end if

      varnam="ww"
      status = nf90_inq_varid(ncid, "ww", ww_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         ww_id=-1
      end if

      varnam="uuEx"
      status = nf90_inq_varid(ncid, "uuEx", uuEx_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         uuEx_id=-1
      end if

      varnam="vvEx"
      status = nf90_inq_varid(ncid, "vvEx", vvEx_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         vvEx_id=-1
      end if

!  tke is required
      varnam="tke"
      status = nf90_inq_varid(ncid, "tke", tke_id)
      if (status .NE. NF90_NOERR) go to 10

!  eps is required
      varnam="eps"
      status = nf90_inq_varid(ncid, "eps", eps_id)
      if (status .NE. NF90_NOERR) go to 10

!  num is required
      varnam="num"
      status = nf90_inq_varid(ncid, "num", num_id)
      if (status .NE. NF90_NOERR) go to 10

!  nuh is required
      varnam="nuh"
      status = nf90_inq_varid(ncid, "nuh", nuh_id)
      if (status .NE. NF90_NOERR) go to 10

      varnam="ho"
      status = nf90_inq_varid(ncid, "ho", ho_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         ho_id=-1
      else
         restart_with_ho=.true.
      endif

      varnam="hn"
      status = nf90_inq_varid(ncid, "hn", hn_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         hn_id=-1
      else
         restart_with_hn=.true.
      endif

#ifndef NO_BAROCLINIC
!  T is required
      varnam="T"
      status = nf90_inq_varid(ncid, "T", T_id)
      if (status .NE. NF90_NOERR) go to 10

!  S is required
      varnam="S"
      status = nf90_inq_varid(ncid, "S", S_id)
      if (status .NE. NF90_NOERR) go to 10
#endif
#ifdef SPM
      varnam="spm"
      status = nf90_inq_varid(ncid, "spm", spm_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         spm_id=-1
      end if

      varnam="spmpool"
      status = nf90_inq_varid(ncid, "spmpool", spmpool_id)
      if (status .NE. NF90_NOERR) then
         LEVEL3 'variable missing in restart file. Skipping ',varnam
         spmpool_id=-1
      end if

#endif
#ifdef GETM_BIO
      if (bio_calc .and. bio_init_method .eq. 0) then
         varnam="bio"
         status = nf90_inq_varid(ncid, "bio", bio_id)
         if (status .NE. NF90_NOERR) go to 10
      end if
#endif

#ifdef _FABM_
      if (fabm_calc .and. fabm_init_method .eq. 0) then
         varnam="fabm_pel"
         status = nf90_inq_varid(ncid,varnam,fabm_pel_id)
         if (status .NE. NF90_NOERR) go to 10

         varnam="fabm_ben"
         status = nf90_inq_varid(ncid,varnam,fabm_ben_id)
         if (status .NE. NF90_NOERR) fabm_ben_id=0 !go to 10
      end if
#endif
   end if
#endif

   return

   10 FATAL 'open_restart_ncdf: ',varnam,' ',nf90_strerror(status)
   stop 'open_restart_ncdf'

   return
   end subroutine open_restart_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Karsten Bolding (BBH)                           !
!-----------------------------------------------------------------------
