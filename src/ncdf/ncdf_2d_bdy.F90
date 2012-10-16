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
   use domain,  only: need_2d_bdy_elev,need_2d_bdy_u,need_2d_bdy_v
   use domain, only: nsbv,nsbvl,nbdy,NWB,NNB,NEB,NSB,bdy_index,bdy_index_l
   use domain, only: wi,wfj,wlj,nj,nfi,nli,ei,efj,elj,sj,sfi,sli
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
   integer                             :: loop0
   REALTYPE                            :: offset=_ZERO_

   REALTYPE,dimension(:),allocatable   :: bdy_times
   REALTYPE,dimension(:),pointer       :: bdy_data_new,d_bdy_data
   REALTYPE,dimension(:),pointer       :: bdy_data_u_new,d_bdy_data_u
   REALTYPE,dimension(:),pointer       :: bdy_data_v_new,d_bdy_data_v
   REALTYPE,dimension(:),allocatable   :: wrk
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
   subroutine init_2d_bdy_ncdf(fname,loop)
!
! !DESCRIPTION:
!  kurt,kurt
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fname
   integer, intent(in)                 :: loop
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
      bdy_len = nsbv
   end if

#if 0
   err = nf90_inq_varid(ncid,'julday',jul_id)
   if (err .NE. NF90_NOERR) go to 10

   err = nf90_inq_varid(ncid,'secs',secs_id)
   if (err .NE. NF90_NOERR) go to 10
#endif

   err = nf90_inq_varid(ncid,'time',time_id)
   if (err .NE. NF90_NOERR) go to 10
   err =  nf90_get_att(ncid,time_id,'units',units)
   if (err .NE. NF90_NOERR) go to 10

   allocate(bdy_times(nsets),stat=err)
   if (err /= 0) stop 'init_2d_bdy_ncdf: Error allocating memory (bdy_times)'
   err = nf90_get_var(ncid,time_id,bdy_times)
   if (err .NE. NF90_NOERR) go to 10

   loop0 = loop-1
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

   if ( need_2d_bdy_elev ) then
      err = nf90_inq_varid(ncid,'elev',elev_id)
      if (err .NE. NF90_NOERR) go to 10
      allocate(bdy_data_new(nsbvl),stat=err)
      if (err /= 0) stop 'init_2d_bdy_ncdf: Error allocating memory (bdy_data_new)'
      allocate(d_bdy_data(nsbvl),stat=err)
      if (err /= 0) stop 'init_2d_bdy_ncdf: Error allocating memory (d_bdy_data)'
   end if

   if ( need_2d_bdy_u ) then
      err = nf90_inq_varid(ncid,'u',u_id)
      if (err .NE. NF90_NOERR) go to 10
      allocate(bdy_data_u_new(nsbvl),stat=err)
      if (err /= 0) stop 'init_2d_bdy_ncdf: Error allocating memory (bdy_data_u_new)'
      allocate(d_bdy_data_u(nsbvl),stat=err)
      if (err /= 0) stop 'init_2d_bdy_ncdf: Error allocating memory (d_bdy_data_u)'
   end if

   if ( need_2d_bdy_v ) then
      err = nf90_inq_varid(ncid,'v',v_id)
      if (err .NE. NF90_NOERR) go to 10
      allocate(bdy_data_v_new(nsbvl),stat=err)
      if (err /= 0) stop 'init_2d_bdy_ncdf: Error allocating memory (bdy_data_v_new)'
      allocate(d_bdy_data_v(nsbvl),stat=err)
      if (err /= 0) stop 'init_2d_bdy_ncdf: Error allocating memory (d_bdy_data_v)'
   end if

   if (err .NE. NF90_NOERR) go to 10
   allocate(wrk(bdy_len),stat=err)
   wrk = _ZERO_

   call do_2d_bdy_ncdf(loop0)

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
   integer, intent(in)           :: loop
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer,save                  :: indx=1,start(2),edges(2)
   integer                       :: i,err
   logical                       :: first=.true.
   REALTYPE                      :: t,t_minus_t2
   REALTYPE, save                :: t1,t2= -_ONE_,deltm1
   REALTYPE,dimension(:),pointer :: bdy_data_old,bdy_data_u_old,bdy_data_v_old
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   write(debug,*) 'do_2d_bdy_ncdf (NetCDF)'
#endif

   t = (loop-loop0)*dtm


   if(t .gt. t2) then

      call write_time_string()
      LEVEL3 timestr,': reading 2D boundary data ...'
      t1 = t2
      do i=indx+1,nsets
         t2 = bdy_times(i) - offset
         if(t2 .gt. t) then
            EXIT
         end if
      end do

      if (first) then
         indx = i-1
         t2 = bdy_times(indx) - offset
         start(1) = 1 ; edges(1) = bdy_len
         edges(2) = 1
         first = .false.
      else
         indx = i
      end if
      start(2) = indx

      if ( elev_id .ne. -1 ) then
         err = nf90_get_var(ncid,elev_id,wrk,start,edges)
         if(err .NE. NF90_NOERR) go to 10
         call grid_2d_bdy_data_ncdf(bdy_len,wrk,nsbvl,bdy_data)
         bdy_data_old=>bdy_data_new;bdy_data_new=>bdy_data;bdy_data=>d_bdy_data;d_bdy_data=>bdy_data_old
         d_bdy_data = bdy_data_new - bdy_data_old
      end if
      if ( u_id .ne. -1 ) then
         err = nf90_get_var(ncid,u_id,wrk,start,edges)
         if(err .NE. NF90_NOERR) go to 10
         call grid_2d_bdy_data_ncdf(bdy_len,wrk,nsbvl,bdy_data_u)
         bdy_data_u_old=>bdy_data_u_new;bdy_data_u_new=>bdy_data_u;bdy_data_u=>d_bdy_data_u;d_bdy_data_u=>bdy_data_u_old
         d_bdy_data_u = bdy_data_u_new - bdy_data_u_old
      end if
      if ( v_id .ne. -1 ) then
         err = nf90_get_var(ncid,v_id,wrk,start,edges)
         if(err .NE. NF90_NOERR) go to 10
         call grid_2d_bdy_data_ncdf(bdy_len,wrk,nsbvl,bdy_data_v)
         bdy_data_v_old=>bdy_data_v_new;bdy_data_v_new=>bdy_data_v;bdy_data_v=>d_bdy_data_v;d_bdy_data_v=>bdy_data_v_old
         d_bdy_data_v = bdy_data_v_new - bdy_data_v_old
      end if

      deltm1 = _ONE_ / (t2 - t1)

   end if


   t_minus_t2 = t - t2

   if ( need_2d_bdy_elev ) then 
      bdy_data = bdy_data_new + d_bdy_data*deltm1*t_minus_t2
   end if
   if ( need_2d_bdy_u ) then
      bdy_data_u = bdy_data_u_new + d_bdy_data_u*deltm1*t_minus_t2
   end if
   if ( need_2d_bdy_v ) then
      bdy_data_v = bdy_data_v_new + d_bdy_data_v*deltm1*t_minus_t2
   end if


#ifdef DEBUG
   write(debug,*) 'Leaving do_2d_bdy_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'do_2d_bdy_data_ncdf: ',nf90_strerror(err)
   stop
   end subroutine do_2d_bdy_ncdf
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: grid_2d_bdy_data_ncdf -
!
! !INTERFACE:
   subroutine grid_2d_bdy_data_ncdf(nsbv,wrk,nsbvl,col)
!
! !DESCRIPTION:
!  kurt,kurt
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)  :: nsbv,nsbvl
   REALTYPE,intent(in) :: wrk(nsbv)
!
! !OUTPUT PARAMETERS:
   REALTYPE,intent(out) :: col(nsbvl)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer             :: i,j,k,kl,l,n
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   write(debug,*) 'grid_2d_bdy_data_ncdf'
#endif

   l = 0
   do n=1,NWB
      l = l+1
      k = bdy_index(l)
      kl = bdy_index_l(l)
      i = wi(n)
      do j=wfj(n),wlj(n)
         col(kl) = wrk(k)
         k = k+1
         kl = kl + 1
      end do
   end do
   do n = 1,NNB
      l = l+1
      k = bdy_index(l)
      kl = bdy_index_l(l)
      j = nj(n)
      do i = nfi(n),nli(n)
         col(kl) = wrk(k)
         k = k+1
         kl = kl + 1
      end do
   end do
   do n=1,NEB
      l = l+1
      k = bdy_index(l)
      kl = bdy_index_l(l)
      i = ei(n)
      do j=efj(n),elj(n)
         col(kl) = wrk(k)
         k = k+1
         kl = kl + 1
      end do
   end do
   do n = 1,NSB
      l = l+1
      k = bdy_index(l)
      kl = bdy_index_l(l)
      j = sj(n)
      do i = sfi(n),sli(n)
         col(kl) = wrk(k)
         k = k+1
         kl = kl + 1
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving grid_2d_bdy_data_ncdf()'
   write(debug,*)
#endif
   return
   end subroutine grid_2d_bdy_data_ncdf
!EOC

!-----------------------------------------------------------------------

   end module ncdf_2d_bdy

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
