!$Id: ncdf_3d_bdy.F90,v 1.3 2003-04-07 16:19:52 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  ncdf_3d_bdy - input in NetCDF format
!
! !INTERFACE:
   module ncdf_3d_bdy
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax
   use domain, only: nsbv,NWB,NNB,NEB,NSB
   use domain, only: wi,wfj,wlj,nj,nfi,nli,ei,efj,elj,sj,sfi,sli
   use domain, only: H
   use m2d, only: dtm
   use variables_3d, only: hn
   use bdy_3d, only: T_bdy,S_bdy
   use time, only: string_to_julsecs,TimeDiff,julianday,secondsofday
   IMPLICIT NONE
!
   private
!
   public :: init_3d_bdy_ncdf,do_3d_bdy_ncdf
!
! !PRIVATE DATA MEMBERS:
   integer 	:: ncid
   integer	:: time_id,temp_id,salt_id
   integer	:: start(4),edges(4)
   integer	:: time_dim,time_len,bdy_len
   logical	:: climatology=.false.,from_3d_fields=.false.
   REALTYPE	:: offset
   REAL_4B, allocatable	:: bdy_times(:),wrk(:)
   REALTYPE, allocatable, dimension(:,:)	:: T_old, T_new
   REALTYPE, allocatable, dimension(:,:)	:: S_old, S_new
   REALTYPE, allocatable, dimension(:,:,:)	:: T_bdy_clim,S_bdy_clim
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ncdf_3d_bdy.F90,v $
!  Revision 1.3  2003-04-07 16:19:52  kbk
!  parallel support
!
!  Revision 1.1.1.1  2002/05/02 14:01:49  gotm
!  recovering after CVS crash
!
!  Revision 1.1  2001/10/17 13:28:27  bbh
!  Initial import
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_3d_bdy_ncdf -
!
! !INTERFACE:
   subroutine init_3d_bdy_ncdf(fname)
!
! !DESCRIPTION:
!  kurt,kurt
!
! !USES:
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
   character(len=256)	:: units
   integer	:: j1,s1
   integer	:: ndims
   integer, allocatable, dimension(:)	:: dim_ids,dim_len
   character(len=16), allocatable 	:: dim_name(:)

   REAL_4B, allocatable, dimension(:,:,:)	:: sdum,tdum
   REAL_4B, allocatable, dimension(:)		:: zlev

   integer	:: rc,err
   integer	:: i,j,k,l,n,id
!EOP
!-------------------------------------------------------------------------
!BOC
   include "netcdf.inc"
#ifdef DEBUG
   write(debug,*) 'ncdf_init_3d_bdy (NetCDF)'
   write(debug,*) 'Reading from: ',trim(fname)
#endif

   LEVEL3 'init_3d_bdy_ncdf'

   err = nf_open(fname,NCNOWRIT,ncid)
   if (err .NE. NF_NOERR) go to 10

   err = nf_inq_ndims(ncid,ndims)
   if (err .NE. NF_NOERR) go to 10

   allocate(dim_ids(ndims),stat=rc)
   if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (dim_ids)'

   allocate(dim_len(ndims),stat=rc)
   if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (dim_len)'

   allocate(dim_name(ndims),stat=rc)
   if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (dim_name)'

   do n=1,ndims
      err = nf_inq_dim(ncid, n, dim_name(n), dim_len(n))
      if (err .NE. NF_NOERR) go to 10
      STDERR n,dim_name(n), dim_len(n)
   end do

   if(ndims .eq. 4) then
!     We are reading boundary values from a full 3D field
!     We assume COARDS conventions
!     1 -> lon,x-axis
!     2 -> lat,y-axis
!     3 -> zax,levels
!     4 -> time
      LEVEL4 'boundary data from 3D fields'
      from_3d_fields=.true.
      time_dim = 4

      allocate(zlev(dim_len(3)),stat=rc)
      if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (zlev)'
      err = nf_inq_varid(ncid, dim_name(3), id)
      if (err .ne. NF_NOERR) go to 10
      err = nf_get_var_real(ncid,id,zlev)
      if (err .ne. NF_NOERR) go to 10
      zlev = -_ONE_*zlev

      allocate(wrk(dim_len(3)),stat=rc)
      if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (wrk)'

      allocate(tdum(imin:imax,jmin:jmax,dim_len(3)),stat=rc)
      if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (tdum)'

      allocate(sdum(imin:imax,jmin:jmax,dim_len(3)),stat=rc)
      if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (sdum)'
   else
      LEVEL4 'special boundary data file'
!     We are reading boundary values from a special boundary data file
!     1 -> zax,levels
!     2 -> bdy_points
!     3 -> time
      bdy_len = dim_len(2)
      time_dim = 3
   end if
   time_len = dim_len(time_dim)
   if( time_len .le. 12) then
      climatology=.true.
      LEVEL4 'Assuming climatolgical 3D boundary conditions'
      LEVEL4 '# of times = ',time_len
   end if

   err = nf_inq_varid(ncid,'temp',temp_id)
   if (err .NE. NF_NOERR) go to 10

   err = nf_inq_varid(ncid,'salt',salt_id)
   if (err .NE. NF_NOERR) go to 10
   
   if (climatology) then
      allocate(T_bdy_clim(time_len,nsbv,0:kmax),stat=rc)
      if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (T_bdy_clim)'

      allocate(S_bdy_clim(time_len,nsbv,0:kmax),stat=rc)
      if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (S_bdy_clim)'

      if (from_3d_fields) then
         start(1) = imin; edges(1) = imax-imin+1;
         start(2) = jmin; edges(2) = jmax-jmin+1;
         start(3) = 1; edges(3) = dim_len(3);
         edges(4) = 1
         do l=1,time_len
            start(4) = l
            err = nf_get_vara_real(ncid,temp_id,start,edges,tdum)
            if (err .ne. NF_NOERR) go to 10
            err = nf_get_vara_real(ncid,salt_id,start,edges,sdum)
            if (err .ne. NF_NOERR) go to 10
            k = 0
            do n=1,NWB
               i = wi(n)
               do j=wfj(1),wlj(1)
                  k = k+1
                  wrk(:) = sdum(i,j,:)
                  call interpol(zlev,wrk,H(i,j),kmax,hn(i,j,:),S_bdy_clim(l,k,:))
                  wrk(:) = tdum(i,j,:)
                  call interpol(zlev,wrk,H(i,j),kmax,hn(i,j,:),T_bdy_clim(l,k,:))
               end do
            end do
            do n = 1,NNB
               j = nj(n)
               do i = nfi(n),nli(n)
                  k = k+1
                  wrk(:) = sdum(i,j,:)
                  call interpol(zlev,wrk,H(i,j),kmax,hn(i,j,:),S_bdy_clim(l,k,:))
                  wrk(:) = tdum(i,j,:)
                  call interpol(zlev,wrk,H(i,j),kmax,hn(i,j,:),T_bdy_clim(l,k,:))
               end do
            end do
            do n=1,NEB
               i = ei(n)
               do j=efj(1),elj(1)
                  k = k+1
                  wrk(:) = sdum(i,j,:)
                  call interpol(zlev,wrk,H(i,j),kmax,hn(i,j,:),S_bdy_clim(l,k,:))
                  wrk(:) = tdum(i,j,:)
                  call interpol(zlev,wrk,H(i,j),kmax,hn(i,j,:),T_bdy_clim(l,k,:))
               end do
            end do
            do n = 1,NSB
               j = sj(n)
               do i = sfi(n),sli(n)
                  k = k+1
                  wrk(:) = sdum(i,j,:)
                  call interpol(zlev,wrk,H(i,j),kmax,hn(i,j,:),S_bdy_clim(l,k,:))
                  wrk(:) = tdum(i,j,:)
                  call interpol(zlev,wrk,H(i,j),kmax,hn(i,j,:),T_bdy_clim(l,k,:))
               end do
            end do
         end do
      else
         start(1) = 1; edges(1) = kmax+1;
         start(2) = 1; edges(2) = nsbv;
         edges(3) = 1
	 STDERR 'ncdf_init_3d_bdy - not finished yet'
	 stop
      end if
      err = nf_close(ncid)
   else
      err = nf_inq_varid(ncid,'time',time_id)
      if (err .NE. NF_NOERR) go to 10
   
      err =  nf_get_att_text(ncid,time_id,'units',units)
      if (err .NE. NF_NOERR) go to 10
 
      allocate(bdy_times(time_len),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (bdy_times)'
   
      err = nf_get_var_real(ncid,time_id,bdy_times)
      if (err .NE. NF_NOERR) go to 10
   
      call string_to_julsecs(units,j1,s1)
      offset = TimeDiff(julianday,secondsofday,j1,s1)
      if( offset .lt. _ZERO_ ) then
         FATAL 'Model simulation starts before available boundary data'
         stop 'init_3d_bdy_ncdf'
      else
         LEVEL3 'Boundary offset time ',offset
      end if

      allocate(T_old(nsbv,0:kmax),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (T_old)'
      allocate(T_new(nsbv,0:kmax),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (T_new)'

      allocate(S_old(nsbv,0:kmax),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (S_old)'
      allocate(S_new(nsbv,0:kmax),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (S_new)'

      allocate(wrk(nsbv*(kmax+1)),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (wrk)'
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving init_3d_bdy_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'init_3d_bdy_ncdf: ',nf_strerror(err)
   stop
   end subroutine init_3d_bdy_ncdf
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: do_3d_bdy_ncdf -
!
! !INTERFACE:
   subroutine do_3d_bdy_ncdf(loop)
!
! !DESCRIPTION:
!  kurt,kurt
!
! !USES:
   use time, only: day,month,secondsofday,days_in_mon,leapyear,secsprday
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
!  $Log: ncdf_3d_bdy.F90,v $
!  Revision 1.3  2003-04-07 16:19:52  kbk
!  parallel support
!
!  Revision 1.1.1.1  2002/05/02 14:01:49  gotm
!  recovering after CVS crash
!
!  Revision 1.1  2001/10/17 13:28:27  bbh
!  Initial import
!
!  Revision 1.1  2001/05/14 12:45:56  bbh
!  Introduced module ncdf_2d_bdy
!
! !LOCAL VARIABLES:
#if 0
   integer,save	:: i,n
   integer,save	:: j,k,indx
   logical	:: first=.true.
   REALTYPE	:: t
   REALTYPE, save	:: t1,t2= -_ONE_,loop0
#endif
   integer	:: err
   REALTYPE	:: rat
   integer	:: monthsecs,prev,this,next
!EOP
!-------------------------------------------------------------------------
!BOC
   include "netcdf.inc"
#ifdef DEBUG
   write(debug,*) 'do_3d_bdy_ncdf (NetCDF)'
#endif

   if ( climatology ) then
      if (time_len .eq. 12) then
         this = month
         monthsecs = secsprday*days_in_mon(leapyear,month)
         rat=((day-1)*secsprday+secondsofday)/float(monthsecs)
         next=this+1
         if (next .gt. time_len) next=1
         prev=this-1
         if (prev .eq. 0) prev=time_len
      else
         STDERR 'do_3d_bdy_ncdf: climatology time_len .ne. 12'
	 stop
      end if
      S_bdy=(1.-rat)*0.5*(S_bdy_clim(prev,:,:)+S_bdy_clim(this,:,:))  &
         +     rat*0.5*(S_bdy_clim(next,:,:)+S_bdy_clim(this,:,:))
      T_bdy=(1.-rat)*0.5*(T_bdy_clim(prev,:,:)+T_bdy_clim(this,:,:))  &
         +     rat*0.5*(T_bdy_clim(next,:,:)+T_bdy_clim(this,:,:))
   else
   end if
#if 0
   start(1) = 1   ; edges(1) = kmax+1
   start(2) = 1   ; edges(2) = bdy_len

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
      start(3) = i-1 ; edges(3) = 1
      err = nf_get_vara_real(ncid,temp_id,start,edges,wrk)
      if(err .NE. NF_NOERR) go to 10
indx = 0
do k=1,89
do j=0,kmax
indx = indx+1
T_old(k,j) = wrk(indx)
end do
end do
      start(3) = i ; edges(3) = 1
      err = nf_get_vara_real(ncid,temp_id,start,edges,wrk)
      if(err .NE. NF_NOERR) go to 10
indx = 0
do k=1,89
do j=0,kmax
indx = indx+1
T_new(k,j) = wrk(indx)
end do
end do
      start(3) = i-1 ; edges(3) = 1
      err = nf_get_vara_real(ncid,salt_id,start,edges,wrk)
      if(err .NE. NF_NOERR) go to 10
indx = 0
do k=1,89
do j=0,kmax
indx = indx+1
S_old(k,j) = wrk(indx)
end do
end do

      start(3) = i ; edges(3) = 1
      err = nf_get_vara_real(ncid,salt_id,start,edges,wrk)
      if(err .NE. NF_NOERR) go to 10
indx = 0
do k=1,89
do j=0,kmax
indx = indx+1
S_new(k,j) = wrk(indx)
end do
end do
   end if

   if(t .gt. t2) then
STDERR 'New 3D boundary data'
      do i=1,n
         if(bdy_times(i) .gt. real(t + offset)) then
            EXIT
         end if
      end do
      t1 = bdy_times(i-1) - offset
      t2 = bdy_times(i) - offset

      T_old = T_new
      start(3) = i-1 ; edges(3) = 1
      err = nf_get_vara_real(ncid,temp_id,start,edges,wrk)
      if(err .NE. NF_NOERR) go to 10
indx = 0
do k=1,89
do j=0,kmax
indx = indx+1
T_new(k,j) = wrk(indx)
end do
end do

      S_old = S_new
      start(3) = i ; edges(3) = 1
      err = nf_get_vara_real(ncid,salt_id,start,edges,wrk)
      if(err .NE. NF_NOERR) go to 10
indx = 0
do k=1,89
do j=0,kmax
indx = indx+1
S_new(k,j) = wrk(indx)
end do
end do
   end if

   T_bdy = T_old + (T_new - T_old)*(t-t1)/(t2-t1)
   S_bdy = S_old + (S_new - S_old)*(t-t1)/(t2-t1)
#endif
#ifdef DEBUG
   write(debug,*) 'Leaving do_3d_bdy_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'do_3d_bdy_ncdf: ',nf_strerror(err)
   stop
   end subroutine do_3d_bdy_ncdf
!EOC

!-----------------------------------------------------------------------

! quick and dirty - should be merged with kbk_interpol.F90 and
! grid_interpol.F90

subroutine interpol(zlev,wrk,depth,kmax,zm,col)

REAL_4B		:: zlev(18),wrk(18)
REALTYPE	:: depth
integer		:: kmax
REALTYPE	:: zm(0:kmax),col(0:kmax)

REALTYPE	:: zmodel(kmax),rat
integer		:: k,n,nn

   zmodel(1) = -depth + 0.5*zm(1)
   do k=2,kmax
      zmodel(k) = zmodel(k-1) + 0.5*(zm(k-1)+zm(k))
   end do

   do k=kmax,1,-1
      if (zmodel(k) .ge. zlev(1)) col(k) = wrk(1)
   end do

   do k=1,kmax
      if (zmodel(k) .le. zlev(18)) col(k) = wrk(18)
   end do

   do k=1,kmax
      if (zmodel(k) .gt. zlev(18) .and. zmodel(k) .lt. zlev(1)) then
         nn=19
224      nn=nn-1
         if(zlev(nn) .le. zmodel(k)) goto 224
         rat = (zmodel(k)-zlev(nn+1))/(zlev(nn)-zlev(nn+1))
         col(k) = (_ONE_-rat)*wrk(nn+1)+rat*wrk(nn)
      end if
   end do
   col(0)=col(1)
end subroutine interpol

   end module ncdf_3d_bdy

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
