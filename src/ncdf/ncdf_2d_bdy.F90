!$Id: ncdf_2d_bdy.F90,v 1.8 2008-12-12 06:06:50 kb Exp $
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
   use netcdf
!KB   use m2d, only: dtm,bdy_times,bdy_old,bdy_new,bdy_data
   use m2d, only: dtm,bdy_times,bdy_data,bdy_data_u,bdy_data_v
   use time, only: string_to_julsecs,time_diff,add_secs
   use time, only: julianday,secondsofday,juln,secsn
   use time, only: write_time_string,timestr
   use domain,  only: need_2d_bdy_elev,need_2d_bdy_u,need_2d_bdy_v
   IMPLICIT NONE
!
   private
!
   public                              :: init_2d_bdy_ncdf,do_2d_bdy_ncdf
!
! !PRIVATE DATA MEMBERS:
   integer                             :: ncid
   integer                             :: time_id,elev_id=-1,nsets,bdy_len
   integer                             :: u_id=-1, v_id=-1
   integer                             :: start(2),edges(2)
   REALTYPE                            :: offset

   REAL_4B                            :: bdy_old(1500)
   REAL_4B                            :: bdy_new(1500)
   REAL_4B                            :: bdy_old_u(1500)
   REAL_4B                            :: bdy_new_u(1500)
   REAL_4B                            :: bdy_old_v(1500)
   REAL_4B                            :: bdy_new_v(1500)
!
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ncdf_2d_bdy.F90,v $
!  Revision 1.8  2008-12-12 06:06:50  kb
!  fixed serious error in specified u boundary velocity
!
!  Revision 1.7  2008-12-09 00:31:58  kb
!  added new 2D open boundaries
!
!  Revision 1.6  2007-09-30 13:00:43  kbk
!  prints real time as part of progessoutput
!
!  Revision 1.5  2005-05-04 11:45:29  kbk
!  adding model time stamp on IO
!
!  Revision 1.4  2004/04/06 16:32:29  kbk
!  TimeDiff --> time_diff
!
!  Revision 1.3  2003/04/23 11:54:03  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 12:49:47  kbk
!  dont need variables_3d
!
!  Revision 1.1.1.1  2002/05/02 14:01:46  gotm
!  recovering after CVS crash
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
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fname
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
   integer                   :: err,rec_id,bdy_id
   character(len=256)        :: units
   character(len=19)         :: tbuf
   integer                   :: j1,s1,j2,s2
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   write(debug,*) 'init_2d_bdy_ncdf'
   write(debug,*) 'Reading from: ',trim(fname)
#endif
   LEVEL3 'init_2d_bdy_ncdf'

   err = nf90_open(fname,NF90_NOWRITE,ncid)
   if (err .NE. NF90_NOERR) go to 10

   err = nf90_inquire(ncid, unlimitedDimID = rec_id)
   if (err .ne. NF90_NOERR) go to 10

   err = nf90_inquire_dimension(ncid, rec_id, len = nsets)
   if (err .ne. NF90_NOERR) go to 10

   err = nf90_inq_dimid(ncid, 'nbdyp', bdy_id)
   if (err .ne. NF90_NOERR)  go to 10

   err = nf90_inquire_dimension(ncid, bdy_id, len = bdy_len)
   if (err .ne. NF90_NOERR) go to 10

   err = nf90_inq_varid(ncid,'time',time_id)
   if (err .NE. NF90_NOERR) go to 10

   err =  nf90_get_att(ncid,time_id,'units',units)
   if (err .NE. NF90_NOERR) go to 10

#if 0
   err = nf90_inq_varid(ncid,'julday',jul_id)
   if (err .NE. NF90_NOERR) go to 10

   err = nf90_inq_varid(ncid,'secs',secs_id)
   if (err .NE. NF90_NOERR) go to 10
#endif

   if ( need_2d_bdy_elev ) then
      err = nf90_inq_varid(ncid,'elev',elev_id)
      if (err .NE. NF90_NOERR) go to 10
   end if

   if ( need_2d_bdy_u ) then
      err = nf90_inq_varid(ncid,'u',u_id)
      if (err .NE. NF90_NOERR) go to 10
   end if

   if ( need_2d_bdy_v ) then
      err = nf90_inq_varid(ncid,'v',v_id)
      if (err .NE. NF90_NOERR) go to 10
   end if

   allocate(bdy_times(nsets),stat=err)
   if (err /= 0) stop 'init_2d_bdy_ncdf: Error allocating memory (bdy_times)'

   err = nf90_get_var(ncid,time_id,bdy_times)
   if (err .NE. NF90_NOERR) go to 10

   call string_to_julsecs(units,j1,s1)
   offset = time_diff(julianday,secondsofday,j1,s1)
   if( offset .lt. bdy_times(1) ) then
      FATAL 'Model simulation starts before available boundary data'
      call write_time_string(julianday,secondsofday,tbuf)
      FATAL 'Simulation starts: ',tbuf
      call add_secs(j1,s1,nint(bdy_times(1)),j2,s2)
      call write_time_string(j2,s2,tbuf)
      FATAL 'Datafile starts:   ',tbuf
      stop 'init_2d_bdy_ncdf'
   else
      LEVEL3 'Boundary offset time ',offset
   end if

!  check if the bdy data file is long enough
   if( time_diff(juln,secsn,j1,s1) .gt. bdy_times(nsets) ) then
      FATAL 'Not enough 2D boundary data in file'
      call write_time_string(juln,secsn,tbuf)
      FATAL 'Simulation ends: ',tbuf
      call add_secs(j1,s1,nint(bdy_times(nsets)),j2,s2)
      call write_time_string(j2,s2,tbuf)
      FATAL 'Datafile ends:   ',tbuf
      stop 'init_2d_bdy_ncdf'
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving init_2d_bdy_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'init_2d_bdy_ncdf: ',nf90_strerror(err)
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
   integer, intent(in)                 :: loop
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ncdf_2d_bdy.F90,v $
!  Revision 1.8  2008-12-12 06:06:50  kb
!  fixed serious error in specified u boundary velocity
!
!  Revision 1.7  2008-12-09 00:31:58  kb
!  added new 2D open boundaries
!
!  Revision 1.6  2007-09-30 13:00:43  kbk
!  prints real time as part of progessoutput
!
!  Revision 1.5  2005-05-04 11:45:29  kbk
!  adding model time stamp on IO
!
!  Revision 1.4  2004/04/06 16:32:29  kbk
!  TimeDiff --> time_diff
!
!  Revision 1.3  2003/04/23 11:54:03  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 12:49:47  kbk
!  dont need variables_3d
!
!  Revision 1.1.1.1  2002/05/02 14:01:46  gotm
!  recovering after CVS crash
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
   integer,save              :: i,n
   integer                   :: err
   logical                   :: first=.true.
   REALTYPE                  :: t
   REALTYPE, save            :: t1,t2= -_ONE_,loop0
!
!EOP
!-------------------------------------------------------------------------
!BOC
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

      if ( need_2d_bdy_elev ) then
         start(2) = i-1 ; edges(2) = 1
         err = nf90_get_var(ncid,elev_id,bdy_old,start,edges)
         if(err .NE. NF90_NOERR) go to 10

         start(2) = i ; edges(2) = 1
         err = nf90_get_var(ncid,elev_id,bdy_new,start,edges)
         if(err .NE. NF90_NOERR) go to 10
      end if

      if ( need_2d_bdy_u ) then
         start(2) = i-1 ; edges(2) = 1
         err = nf90_get_var(ncid,u_id,bdy_old_u,start,edges)
         if(err .NE. NF90_NOERR) go to 10

         start(2) = i ; edges(2) = 1
         err = nf90_get_var(ncid,u_id,bdy_new_u,start,edges)
         if(err .NE. NF90_NOERR) go to 10
      end if

      if ( need_2d_bdy_v ) then
         start(2) = i-1 ; edges(2) = 1
         err = nf90_get_var(ncid,v_id,bdy_old_v,start,edges)
         if(err .NE. NF90_NOERR) go to 10

         start(2) = i ; edges(2) = 1
         err = nf90_get_var(ncid,v_id,bdy_new_v,start,edges)
         if(err .NE. NF90_NOERR) go to 10
      end if

   end if

   if(t .gt. t2) then
      do i=1,n
         if(bdy_times(i) .gt. real(t + offset)) then
            EXIT
         end if
      end do

      call write_time_string()
      LEVEL3 timestr,': reading 2D boundary data ...'

      t1 = bdy_times(i-1) - offset
      t2 = bdy_times(i) - offset

      if ( need_2d_bdy_elev ) then
         start(2) = i-1 ; edges(2) = 1
         err = nf90_get_var(ncid,elev_id,bdy_old,start,edges)
         if(err .NE. NF90_NOERR) go to 10

         start(2) = i ; edges(2) = 1
         err = nf90_get_var(ncid,elev_id,bdy_new,start,edges)
         if(err .NE. NF90_NOERR) go to 10
      end if

      if ( need_2d_bdy_u ) then
         start(2) = i-1 ; edges(2) = 1
         err = nf90_get_var(ncid,u_id,bdy_old_u,start,edges)
         if(err .NE. NF90_NOERR) go to 10

         start(2) = i ; edges(2) = 1
         err = nf90_get_var(ncid,u_id,bdy_new_u,start,edges)
         if(err .NE. NF90_NOERR) go to 10
      end if

      if ( need_2d_bdy_v ) then
         start(2) = i-1 ; edges(2) = 1
         err = nf90_get_var(ncid,v_id,bdy_old_v,start,edges)
         if(err .NE. NF90_NOERR) go to 10

         start(2) = i ; edges(2) = 1
         err = nf90_get_var(ncid,v_id,bdy_new_v,start,edges)
         if(err .NE. NF90_NOERR) go to 10
      end if

   end if

   if ( need_2d_bdy_elev ) then
      bdy_data = bdy_old + (bdy_new - bdy_old)*(t-t1)/(t2-t1)
   end if
   if ( need_2d_bdy_u ) then
      bdy_data_u = bdy_old_u + (bdy_new_u - bdy_old_u)*(t-t1)/(t2-t1)
   end if
   if ( need_2d_bdy_v ) then
      bdy_data_v = bdy_old_v + (bdy_new_v - bdy_old_v)*(t-t1)/(t2-t1)
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving do_2d_bdy_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'do_2d_bdy_ncdf: ',nf90_strerror(err)
   stop
   end subroutine do_2d_bdy_ncdf
!EOC

!-----------------------------------------------------------------------

   end module ncdf_2d_bdy

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
