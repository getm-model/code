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
   use variables_2d, only: dtm
   use bdy_2d, only: bdy_data,bdy_data_u,bdy_data_v
   use time, only: string_to_julsecs,time_diff,add_secs
   use time, only: julianday,secondsofday,juln,secsn
   use time, only: write_time_string,timestr
   use domain,  only: nsbv,need_2d_bdy_elev,need_2d_bdy_u,need_2d_bdy_v
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

   REAL_4B,dimension(:),allocatable   :: bdy_times
   REAL_4B,dimension(:),allocatable   :: bdy_old,bdy_new
   REAL_4B,dimension(:),allocatable   :: bdy_old_u,bdy_new_u
   REAL_4B,dimension(:),allocatable   :: bdy_old_v,bdy_new_v
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
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

   if (bdy_len .lt. nsbv) then
      stop 'init_2d_bdy_ncdf: netcdf file does not contain enough bdy points'
   else if (bdy_len .gt. nsbv) then
      LEVEL4 'WARNING: netcdf file contains data for more bdy points'
   end if

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
      allocate(bdy_old(bdy_len),stat=err)
      if (err /= 0) stop 'init_2d_bdy_ncdf: Error allocating memory (bdy_old)'
      allocate(bdy_new(bdy_len),stat=err)
      if (err /= 0) stop 'init_2d_bdy_ncdf: Error allocating memory (bdy_new)'
   end if

   if ( need_2d_bdy_u ) then
      err = nf90_inq_varid(ncid,'u',u_id)
      if (err .NE. NF90_NOERR) go to 10
      allocate(bdy_old_u(bdy_len),stat=err)
      if (err /= 0) stop 'init_2d_bdy_ncdf: Error allocating memory (bdy_old_u)'
      allocate(bdy_new_u(bdy_len),stat=err)
      if (err /= 0) stop 'init_2d_bdy_ncdf: Error allocating memory (bdy_new_u)'
   end if

   if ( need_2d_bdy_v ) then
      err = nf90_inq_varid(ncid,'v',v_id)
      if (err .NE. NF90_NOERR) go to 10
      allocate(bdy_old_v(bdy_len),stat=err)
      if (err /= 0) stop 'init_2d_bdy_ncdf: Error allocating memory (bdy_old_v)'
      allocate(bdy_new_v(bdy_len),stat=err)
      if (err /= 0) stop 'init_2d_bdy_ncdf: Error allocating memory (bdy_new_v)'
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
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer,save              :: i,n
   integer                   :: err
   logical                   :: first=.true.
   REALTYPE                  :: t
   REALTYPE, save            :: t1,t2= -_ONE_,loop0
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
