!$Id: ncdf_2d_bdy.F90,v 1.1 2002-05-02 14:01:46 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  ncdf_2d_bdy - input in NetCDF format
!
! !INTERFACE:
   module ncdf_2d_bdy
!
! !DESCRIPTION:
!
! !USES:
   use m2d, only: dtm,bdy_times,bdy_old,bdy_new,bdy_data
   use time, only: string_to_julsecs,TimeDiff,julianday,secondsofday
   IMPLICIT NONE
!
   private
!
   public :: init_2d_bdy_ncdf,do_2d_bdy_ncdf
!
! !PRIVATE DATA MEMBERS:
   integer 	:: ncid
   integer	:: time_id,elev_id,nsets,bdy_len
   integer	:: start(2),edges(2)
   REALTYPE	:: offset
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ncdf_2d_bdy.F90,v $
!  Revision 1.1  2002-05-02 14:01:46  gotm
!  Initial revision
!
!  Revision 1.8  2001/10/22 11:43:12  bbh
!  Proper check of offset time
!
!  Revision 1.7  2001/10/17 14:57:38  bbh
!  Force offset to _ZERO_ - needs fix
!
!  Revision 1.6  2001/09/19 14:21:13  bbh
!  Cleaning
!
!  Revision 1.5  2001/07/27 06:41:35  bbh
!  Added ncdf_lon_lat.F90
!
!  Revision 1.4  2001/07/26 14:07:18  bbh
!  Typos
!
!  Revision 1.3  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.2  2001/05/18 13:04:39  bbh
!  Cosmetics
!
!  Revision 1.1  2001/05/14 12:45:56  bbh
!  Introduced module ncdf_2d_bdy
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_2d_bdy_ncdf -
!
! !INTERFACE:
   subroutine init_2d_bdy_ncdf(fname)
!
! !DESCRIPTION:
!  kurt,kurt
!
! !USES:
   use variables_3d
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)	:: fname
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See log for module
!
! !LOCAL VARIABLES:
   integer	:: err,rec_id,bdy_id
   character(len=256)	:: units
   integer	:: j1,s1
!
!EOP
!-------------------------------------------------------------------------
!BOC
   include "netcdf.inc"
#ifdef DEBUG
   write(debug,*) 'init_2d_bdy_ncdf'
   write(debug,*) 'Reading from: ',trim(fname)
#endif
   LEVEL3 'init_2d_bdy_ncdf'

   err = nf_open(fname,NCNOWRIT,ncid)
   if (err .NE. NF_NOERR) go to 10

   err = nf_inq_unlimdim(ncid, rec_id)
   if (err .ne. NF_NOERR) go to 10

   err = nf_inq_dimlen(ncid, rec_id, nsets)
   if (err .ne. NF_NOERR) go to 10

   err = nf_inq_dimid(ncid, 'nbdyp', bdy_id)
   if (err .ne. NF_NOERR)  go to 10

   err = nf_inq_dimlen(ncid, bdy_id, bdy_len)
   if (err .ne. NF_NOERR) go to 10

   err = nf_inq_varid(ncid,'time',time_id)
   if (err .NE. NF_NOERR) go to 10

   err =  nf_get_att_text(ncid,time_id,'units',units)
   if (err .NE. NF_NOERR) go to 10

#if 0
   err = nf_inq_varid(ncid,'julday',jul_id)
   if (err .NE. NF_NOERR) go to 10

   err = nf_inq_varid(ncid,'secs',secs_id)
   if (err .NE. NF_NOERR) go to 10
#endif

   err = nf_inq_varid(ncid,'elev',elev_id)
   if (err .NE. NF_NOERR) go to 10

   allocate(bdy_times(nsets),stat=err)
   if (err /= 0) stop 'init_2d_bdy_ncdf: Error allocating memory (bdy_times)'

   err = nf_get_var_real(ncid,time_id,bdy_times)
   if (err .NE. NF_NOERR) go to 10

   call string_to_julsecs(units,j1,s1)
   offset = TimeDiff(julianday,secondsofday,j1,s1)
   if( offset .lt. _ZERO_ ) then
      FATAL 'Model simulation starts before available boundary data'
      stop 'init_2d_bdy_ncdf'
   else
      LEVEL3 'Boundary offset time ',offset
   end if


!  Should have check on length of bdy_file > integration length

#ifdef DEBUG
   write(debug,*) 'Leaving init_2d_bdy_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'init_2d_bdy_ncdf: ',nf_strerror(err)
   stop
   end subroutine init_2d_bdy_ncdf
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_2d_bdy_ncdf -
!
! !INTERFACE:
   subroutine do_2d_bdy_ncdf(loop)
!
! !DESCRIPTION:
!  kurt,kurt
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)	:: loop
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ncdf_2d_bdy.F90,v $
!  Revision 1.1  2002-05-02 14:01:46  gotm
!  Initial revision
!
!  Revision 1.8  2001/10/22 11:43:12  bbh
!  Proper check of offset time
!
!  Revision 1.7  2001/10/17 14:57:38  bbh
!  Force offset to _ZERO_ - needs fix
!
!  Revision 1.6  2001/09/19 14:21:13  bbh
!  Cleaning
!
!  Revision 1.5  2001/07/27 06:41:35  bbh
!  Added ncdf_lon_lat.F90
!
!  Revision 1.4  2001/07/26 14:07:18  bbh
!  Typos
!
!  Revision 1.3  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.2  2001/05/18 13:04:39  bbh
!  Cosmetics
!
!  Revision 1.1  2001/05/14 12:45:56  bbh
!  Introduced module ncdf_2d_bdy
!
! !LOCAL VARIABLES:
   integer,save	:: i,n
   integer	:: err
   logical	:: first=.true.
   REALTYPE	:: t
   REALTYPE, save	:: t1,t2= -_ONE_,loop0
!
!EOP
!-------------------------------------------------------------------------
!BOC
   include "netcdf.inc"
#ifdef DEBUG
   write(debug,*) 'do_2d_bdy_ncdf (NetCDF)'
#endif

   start(1) = 1   ; edges(1) = bdy_len

   if (first) then
      loop0=loop-1
   endif
   t = (loop-loop0)*dtm

   if (first) then
      first = .false.

      n = size(bdy_times)
      do i=1,n
         if (bdy_times(i) .gt. real(t + offset)) then
	    EXIT
         end if
      end do
      t1 = bdy_times(i-1) - offset
      t2 = bdy_times(i) - offset


      start(2) = i-1 ; edges(2) = 1
      err = nf_get_vara_real(ncid,elev_id,start,edges,bdy_old)
      if(err .NE. NF_NOERR) go to 10

      start(2) = i ; edges(2) = 1
      err = nf_get_vara_real(ncid,elev_id,start,edges,bdy_new)
      if(err .NE. NF_NOERR) go to 10
   end if

   if(t .gt. t2) then
      do i=1,n
         if(bdy_times(i) .gt. real(t + offset)) then
	    EXIT
         end if
      end do
      t1 = bdy_times(i-1) - offset
      t2 = bdy_times(i) - offset
      start(2) = i-1 ; edges(2) = 1
      err = nf_get_vara_real(ncid,elev_id,start,edges,bdy_old)
      if(err .NE. NF_NOERR) go to 10

      start(2) = i ; edges(2) = 1
      err = nf_get_vara_real(ncid,elev_id,start,edges,bdy_new)
      if(err .NE. NF_NOERR) go to 10
   end if

   bdy_data = bdy_old + (bdy_new - bdy_old)*(t-t1)/(t2-t1)

#ifdef DEBUG
   write(debug,*) 'Leaving do_2d_bdy_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'do_2d_bdy_ncdf: ',nf_strerror(err)
   stop
   end subroutine do_2d_bdy_ncdf
!EOC

!-----------------------------------------------------------------------

   end module ncdf_2d_bdy

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
