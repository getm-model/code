!$Id: ncdf_3d_bdy.F90,v 1.6 2003-08-03 09:19:41 kbk Exp $
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
   use domain, only: imin,imax,jmin,jmax,kmax,ioff,joff
   use domain, only: nsbv,NWB,NNB,NEB,NSB,bdy_index
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
   public                              :: init_3d_bdy_ncdf,do_3d_bdy_ncdf
!
! !PRIVATE DATA MEMBERS:
   integer                             :: ncid
   integer                             :: time_id,temp_id,salt_id
   integer                             :: start(4),edges(4)
   integer                             :: zax_dim,zax_len
   integer                             :: time_dim,time_len
   logical                             :: climatology=.false.
   logical                             :: from_3d_fields=.false.
   REALTYPE                            :: offset
   REAL_4B, allocatable                :: bdy_times(:),wrk(:)
   REALTYPE, allocatable, dimension(:,:)   :: T_old, T_new
   REALTYPE, allocatable, dimension(:,:)   :: S_old, S_new
   REALTYPE, allocatable, dimension(:,:,:) :: T_bdy_clim,S_bdy_clim
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ncdf_3d_bdy.F90,v $
!  Revision 1.6  2003-08-03 09:19:41  kbk
!  optimised reading of climatological boundary data
!
!  Revision 1.5  2003/05/05 15:44:20  kbk
!  reads boundary values from 3D fields as individual columns
!
!  Revision 1.4  2003/04/23 11:54:03  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.3  2003/04/07 16:19:52  kbk
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
   character(len=256)        :: units
   integer                   :: j1,s1
   integer                   :: ndims
   integer, allocatable, dimension(:):: dim_ids,dim_len
   character(len=16), allocatable :: dim_name(:)

   REAL_4B, allocatable, dimension(:):: zlev

   integer                   :: rc,err
   integer                   :: i,j,k,l,m,n,id
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
      LEVEL4 n,dim_name(n), dim_len(n)
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
      zax_dim = 3
      time_dim = 4
   else
!     We are reading boundary values from a special boundary data file
!     1 -> zax,levels
!     2 -> bdy_points
!     3 -> time
      LEVEL4 'special boundary data file'
      zax_dim = 1
      time_dim = 3
   end if
   zax_len = dim_len(zax_dim)
   time_len = dim_len(time_dim)

   allocate(zlev(zax_len),stat=rc)
   if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (zlev)'

   err = nf_inq_varid(ncid, dim_name(zax_dim), id)
   if (err .ne. NF_NOERR) go to 10

   err = nf_get_var_real(ncid,id,zlev)
   if (err .ne. NF_NOERR) go to 10
   zlev = -_ONE_*zlev

   allocate(wrk(zax_len),stat=rc)
   if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (wrk)'

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

!     we read each boundary column individually
!     here we can read from both a 3D field and from a
!     special ondary data file - only the arguments 'start' and 'edges'
!     varies in the calls to 'nf_get_vara_real()'
!     m counts the time
!     l counts the boundary number
!     k counts the number of the specific point
!     MUST cover the same area as in topo.nc 
      if (from_3d_fields) then
         edges(1) = 1;
         edges(2) = 1;
         start(3) = 1; edges(3) = dim_len(3);
         edges(4) = 1
      else
         start(1) = 1; edges(1) = kmax+1;
         edges(2) = 1;
         edges(3) = 1
      end if

      do m=1,time_len
         start(time_dim) = m

         l = 0
         do n=1,NWB
            l = l+1
            k = bdy_index(l)
            i = wi(n)
            do j=wfj(n),wlj(n)
               if (from_3d_fields) then
                  start(1) = i+ioff ; start(2) = j+joff
               else
                  start(2) = k
               end if

               err = nf_get_vara_real(ncid,salt_id,start,edges,wrk)
               if (err .ne. NF_NOERR) go to 10
               call interpol(zlev,wrk,H(i,j),kmax,hn(i,j,:), &
                             S_bdy_clim(m,k,:))
               err = nf_get_vara_real(ncid,temp_id,start,edges,wrk)
               if (err .ne. NF_NOERR) go to 10
               call interpol(zlev,wrk,H(i,j),kmax,hn(i,j,:), &
                             T_bdy_clim(m,k,:))
               k = k+1
            end do
         end do

         do n = 1,NNB
            l = l+1
            k = bdy_index(l)
            j = nj(n)
            do i = nfi(n),nli(n)
               if (from_3d_fields) then
                  start(1) = i+ioff ; start(2) = j+joff
               else
                  start(2) = k
               end if
               err = nf_get_vara_real(ncid,salt_id,start,edges,wrk)
               if (err .ne. NF_NOERR) go to 10
               call interpol(zlev,wrk,H(i,j),kmax,hn(i,j,:), &
                          S_bdy_clim(m,k,:))
               err = nf_get_vara_real(ncid,temp_id,start,edges,wrk)
               if (err .ne. NF_NOERR) go to 10
               call interpol(zlev,wrk,H(i,j),kmax,hn(i,j,:), &
                             T_bdy_clim(m,k,:))
               k = k+1
            end do
         end do

         do n=1,NEB
            l = l+1
            k = bdy_index(l)
            i = ei(n)
            do j=efj(1),elj(1)
               if (from_3d_fields) then
                  start(1) = i+ioff ; start(2) = j+joff
               else
                  start(2) = k
               end if
               err = nf_get_vara_real(ncid,salt_id,start,edges,wrk)
               if (err .ne. NF_NOERR) go to 10
               call interpol(zlev,wrk,H(i,j),kmax,hn(i,j,:), &
                             S_bdy_clim(m,k,:))
               err = nf_get_vara_real(ncid,temp_id,start,edges,wrk)
               if (err .ne. NF_NOERR) go to 10
               call interpol(zlev,wrk,H(i,j),kmax,hn(i,j,:), &
                             T_bdy_clim(m,k,:))
               k = k+1
            end do
         end do

         do n = 1,NSB
            l = l+1
            k = bdy_index(l)
            j = sj(n)
            do i = sfi(n),sli(n)
               if (from_3d_fields) then
                  start(1) = i+ioff ; start(2) = j+joff
               else
                  start(2) = k
               end if
               err = nf_get_vara_real(ncid,salt_id,start,edges,wrk)
               if (err .ne. NF_NOERR) go to 10
               call interpol(zlev,wrk,H(i,j),kmax,hn(i,j,:), &
                             S_bdy_clim(m,k,:))
               err = nf_get_vara_real(ncid,temp_id,start,edges,wrk)
               if (err .ne. NF_NOERR) go to 10
               call interpol(zlev,wrk,H(i,j),kmax,hn(i,j,:), &
                             T_bdy_clim(m,k,:))
               k = k+1
            end do
         end do
      end do
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
   integer, intent(in)                 :: loop
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ncdf_3d_bdy.F90,v $
!  Revision 1.6  2003-08-03 09:19:41  kbk
!  optimised reading of climatological boundary data
!
!  Revision 1.5  2003/05/05 15:44:20  kbk
!  reads boundary values from 3D fields as individual columns
!
!  Revision 1.4  2003/04/23 11:54:03  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.3  2003/04/07 16:19:52  kbk
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
   integer                   :: err
   REALTYPE                  :: rat
   integer                   :: monthsecs,prev,this,next
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
