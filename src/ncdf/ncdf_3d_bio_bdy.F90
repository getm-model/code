#ifdef _FABM_
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  ncdf_3d_bio_bdy - input in NetCDF format
!
! !INTERFACE:
   module ncdf_3d_bio_bdy
!
! !DESCRIPTION:
!
! !USES:
   use netcdf
   use domain, only: imin,imax,jmin,jmax,kmax,ioff,joff
   use domain, only: nsbv,NWB,NNB,NEB,NSB,bdy_index
   use domain, only: wi,wfj,wlj,nj,nfi,nli,ei,efj,elj,sj,sfi,sli
   use domain, only: H
   use m2d, only: dtm
   use variables_3d, only: hn
   use getm_fabm, only: fabm_calc,model
   use bdy_3d, only: bio_bdy,have_bio_bdy_values
   use time, only: string_to_julsecs,time_diff,julianday,secondsofday
   use time, only: write_time_string,timestr
   IMPLICIT NONE
!
   private
!
   public                              :: init_3d_bio_bdy_ncdf
   public                              :: do_3d_bio_bdy_ncdf
!
! !PRIVATE DATA MEMBERS:
   integer                             :: ncid
   integer                             :: time_id
   integer, allocatable, dimension(:)  :: bio_ids
   integer                             :: start(3),edges(3)
   integer                             :: zax_dim,zax_len
   integer                             :: time_dim,time_len
   logical                             :: climatology=.false.
   REALTYPE                            :: offset
   REAL_4B, allocatable                :: bdy_times(:),wrk(:)
   REAL_4B,  allocatable, dimension(:) :: zlev
   REALTYPE, allocatable, dimension(:,:,:)   :: bio_old, bio_new
   REAL_4B,  allocatable, dimension(:,:,:)   :: bio_wrk
   REALTYPE, allocatable, dimension(:,:,:,:) :: bio_bdy_clim
   integer                             :: npel=-1
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log$
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_3d_bio_bdy_ncdf -
!
! !INTERFACE:
   subroutine init_3d_bio_bdy_ncdf(fname)
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
   integer                   :: rc,err
   integer                   :: i,j,k,l,m,n,o,id
#ifdef _FABM_
   character(len=256)        :: varname
#endif
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   write(debug,*) 'ncdf_init_3d_bdy (NetCDF)'
   write(debug,*) 'Reading from: ',trim(fname)
#endif

   LEVEL3 'init_3d_bio_bdy_ncdf'

   err = nf90_open(fname,NF90_NOWRITE,ncid)
   if (err .NE. NF90_NOERR) go to 10

   err = nf90_inquire(ncid,nDimensions=ndims)
   if (err .NE. NF90_NOERR) go to 10

   allocate(dim_ids(ndims),stat=rc)
   if (rc /= 0) stop 'init_3d_bio_bdy_ncdf: Error allocating memory (dim_ids)'

   allocate(dim_len(ndims),stat=rc)
   if (rc /= 0) stop 'init_3d_bio_bdy_ncdf: Error allocating memory (dim_len)'

   allocate(dim_name(ndims),stat=rc)
   if (rc /= 0) stop 'init_3d_bio_bdy_ncdf: Error allocating memory (dim_name)'

   do n=1,ndims
      err = nf90_inquire_dimension(ncid,n,name=dim_name(n),len=dim_len(n))
      if (err .NE. NF90_NOERR) go to 10
      LEVEL4 n,dim_name(n), dim_len(n)
   end do

!  We are reading boundary values from a special boundary data file
!  The bio variables must all exist and be spanned by dimensions as:
!    1 -> zax,levels
!    2 -> bdy_points
!    3 -> time

!  made to work with the Bodden simulations
#if 0
   zax_dim  = 2
   time_dim = 1
#else
   time_dim = 1
   zax_dim  = 3
#endif

   zax_len = dim_len(zax_dim)
   time_len = dim_len(time_dim)

   allocate(zlev(zax_len),stat=rc)
   if (rc /= 0) stop 'init_3d_bio_bdy_ncdf: Error allocating memory (zlev)'

   err = nf90_inq_varid(ncid, trim(dim_name(zax_dim)), id)
   if (err .ne. NF90_NOERR) go to 10

   err = nf90_get_var(ncid,id,zlev)
   if (err .ne. NF90_NOERR) go to 10
   zlev = -_ONE_*zlev

   allocate(wrk(zax_len),stat=rc)
   if (rc /= 0) stop 'init_3d_bio_bdy_ncdf: Error allocating memory (wrk)'

   if( time_len .eq. 12) then
      climatology=.true.
      LEVEL4 'Assuming climatolgical BIO boundary conditions'
      LEVEL4 '# of times = ',time_len
   end if

   npel = size(model%info%state_variables)

!  npel is known and we can allocate memory for boundary conditions
   allocate(bio_ids(npel),stat=rc)
   if (rc /= 0) stop 'init_3d_bio_bdy_ncdf: Error allocating memory (bio_ids)'
   bio_ids = -1

   allocate(bio_bdy(0:kmax,nsbv,npel),stat=rc)
   if (rc /= 0) stop 'init_bdy_3d: Error allocating memory (bio_bdy)'
   allocate(have_bio_bdy_values(npel),stat=rc)
   if (rc /= 0) stop 'init_bdy_3d: Error allocating memory (have_bio_bdy_values)'
   have_bio_bdy_values = -1

   LEVEL4 'checking available boundary variables:'
   do n=1,npel
      varname = trim(model%info%state_variables(n)%name)
      err = nf90_inq_varid(ncid,trim(varname),bio_ids(n))
      if (err .NE. NF90_NOERR) then
         have_bio_bdy_values(n) = -1
         LEVEL4 trim(varname),': no'
      else
         have_bio_bdy_values(n) = 1
         LEVEL4 trim(varname),': yes'
      end if
   end do

   if (climatology) then

      allocate(bio_bdy_clim(0:kmax,nsbv,time_len,npel),stat=rc)
      if (rc /= 0) stop 'init_3d_bio_bdy_ncdf: Error allocating memory (bio_bdy_clim)'

!     we read each boundary column individually
!     m counts the time
!     l counts the boundary number
!     k counts the number of the specific point
!     MUST cover the same area as in topo.nc
      start(1) = 1; edges(1) = dim_len(zax_dim);
      edges(2) = 1
      edges(3) = 1

      do m=1,time_len
         start(time_dim) = m
         l = 0

         do n=1,NWB
            l = l+1
            i = wi(n)
            do o=1,npel
               k = bdy_index(l)
               if (have_bio_bdy_values(o) .eq. 1 ) then
                  do j=wfj(n),wlj(n)
                     start(2) = k
                     err = nf90_get_var(ncid,bio_ids(o),wrk,start,edges)
                     if (err .ne. NF90_NOERR) go to 10
                     call interpol(zax_len,zlev,wrk,H(i,j),kmax, &
                                   hn(i,j,:),bio_bdy_clim(:,k,m,o))
                     k = k+1
                  end do
               else
                  k = k + (wlj(n)-wfj(n)+1)
               end if
            end do
         end do

         do n = 1,NNB
            l = l+1
            j = nj(n)
            do o=1,npel
               k = bdy_index(l)
               if (have_bio_bdy_values(o) .eq. 1 ) then
                  do i = nfi(n),nli(n)
                     start(2) = k
                     err = nf90_get_var(ncid,bio_ids(o),wrk,start,edges)
                     if (err .ne. NF90_NOERR) go to 10
                     call interpol(zax_len,zlev,wrk,H(i,j),kmax, &
                                   hn(i,j,:),bio_bdy_clim(:,k,m,o))
                     k = k+1
                  end do
               else
                  k = k + (nli(n)-nfi(n)+1)
               end if
            end do
         end do

         do n=1,NEB
            l = l+1
            i = ei(n)
            do o=1,npel
               k = bdy_index(l)
               if (have_bio_bdy_values(o) .eq. 1 ) then
                  do j=efj(1),elj(1)
                     start(2) = k
                     err = nf90_get_var(ncid,bio_ids(o),wrk,start,edges)
                     if (err .ne. NF90_NOERR) go to 10
                     call interpol(zax_len,zlev,wrk,H(i,j),kmax, &
                                   hn(i,j,:),bio_bdy_clim(:,k,o,m))
                     k = k+1
                  end do
               else
                  k = k + (elj(n)-efj(n)+1)
               end if
            end do
         end do

         do n = 1,NSB
            l = l+1
            j = sj(n)
            do o=1,npel
               k = bdy_index(l)
               if (have_bio_bdy_values(o) .eq. 1 ) then
                  do i = sfi(n),sli(n)
                     start(2) = k
                     err = nf90_get_var(ncid,bio_ids(o),wrk,start,edges)
                     if (err .ne. NF90_NOERR) go to 10
                     call interpol(zax_len,zlev,wrk,H(i,j),kmax, &
                                   hn(i,j,:),bio_bdy_clim(:,k,o,m))
                     k = k+1
                  end do
               else
                  k = k + (sli(n)-sfi(n)+1)
               end if
            end do
         end do
      end do
      err = nf90_close(ncid)

   else

      err = nf90_inq_varid(ncid,'time',time_id)
      if (err .NE. NF90_NOERR) go to 10

      err =  nf90_get_att(ncid,time_id,'units',units)
      if (err .NE. NF90_NOERR) go to 10

      allocate(bdy_times(time_len),stat=err)
      if (err /= 0) stop 'init_3d_bio_bdy_ncdf: Error allocating memory (bdy_times)'

      err = nf90_get_var(ncid,time_id,bdy_times)
      if (err .NE. NF90_NOERR) go to 10

      call string_to_julsecs(units,j1,s1)
      offset = time_diff(julianday,secondsofday,j1,s1)
      if( offset .lt. _ZERO_ ) then
         FATAL 'Model simulation starts before available BIO boundary data'
         stop 'init_3d_bio_bdy_ncdf'
      else
         LEVEL3 'Boundary offset time ',offset
      end if

      allocate(bio_old(0:kmax,nsbv,npel),stat=err)
      if (err /= 0) stop 'init_3d_bio_bdy_ncdf: Error allocating memory (bio_old)'
      allocate(bio_new(0:kmax,nsbv,npel),stat=err)
      if (err /= 0) stop 'init_3d_bio_bdy_ncdf: Error allocating memory (bio_new)'
      allocate(bio_wrk(zax_len,nsbv,npel),stat=err)
      if (err /= 0) stop 'init_3d_bio_bdy_ncdf: Error allocating memory (bio_wrk)'
      bio_wrk = _ZERO_

      n = size(bdy_times)
      do i=1,n
         if(bdy_times(i) .ge. real(offset)) then
            EXIT
         end if
      end do

      if(i .gt. 1 .and. bdy_times(i) .gt. real(offset)) then
         i = i-1
      end if

      start(1) = 1; edges(1) = dim_len(zax_dim);
      start(2) = 1; edges(2) = nsbv;
      start(3) = i; edges(3) = 1

      do m=1,npel
         if (have_bio_bdy_values(m) .eq. 1 ) then
            err = nf90_get_var(ncid,bio_ids(m),bio_wrk(:,:,m),start,edges)
            if (err .ne. NF90_NOERR) go to 10
         end if
      end do

      l = 0
      do n=1,NWB
         l = l+1
         i = wi(n)
         do m=1,npel
            k = bdy_index(l)
            do j=wfj(n),wlj(n)
               call interpol(zax_len,zlev,bio_wrk(:,k,m),H(i,j), &
                             kmax,hn(i,j,:),bio_new(:,k,m))
               k = k+1
            end do
         end do
      end do

      do n = 1,NNB
         l = l+1
         j = nj(n)
         do m=1,npel
            k = bdy_index(l)
            do i = nfi(n),nli(n)
               call interpol(zax_len,zlev,bio_wrk(:,k,m),H(i,j), &
                             kmax,hn(i,j,:),bio_new(:,k,m))
               k = k+1
            end do
         end do
      end do

      do n=1,NEB
         l = l+1
         i = ei(n)
         do m=1,npel
            k = bdy_index(l)
            do j=efj(1),elj(1)
               call interpol(zax_len,zlev,bio_wrk(:,k,m),H(i,j), &
                             kmax,hn(i,j,:),bio_new(:,k,m))
               k = k+1
            end do
         end do
      end do

      do n = 1,NSB
         l = l+1
         j = sj(n)
         do m=1,npel
            k = bdy_index(l)
            do i = sfi(n),sli(n)
               call interpol(zax_len,zlev,bio_wrk(:,k,m),H(i,j), &
                             kmax,hn(i,j,:),bio_new(:,k,m))
               k = k+1
            end do
         end do
      end do
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving init_3d_bio_bdy_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'init_3d_bio_bdy_ncdf: ',nf90_strerror(err)
   stop
   end subroutine init_3d_bio_bdy_ncdf
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: do_3d_bio_bdy_ncdf -
!
! !INTERFACE:
   subroutine do_3d_bio_bdy_ncdf(loop)
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
!  $Log$
!
! !LOCAL VARIABLES:
   integer         :: err
   REALTYPE        :: rat
   integer         :: monthsecs,prev,this,next
   logical, save   :: first=.true.
   integer, save   :: loop0
   REALTYPE        :: t
   REALTYPE, save  :: t1=_ZERO_,t2=-_ONE_
   integer         :: i,j,k,l,n,o
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   write(debug,*) 'do_3d_bio_bdy_ncdf (NetCDF)'
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
         STDERR 'do_3d_bio_bdy_ncdf: climatology time_len .ne. 12'
         stop
      end if

      do o=1,npel
         bio_bdy(:,:,o)=(1.-rat)*0.5 &
                        *(bio_bdy_clim(prev,:,:,o)+bio_bdy_clim(this,:,:,o)) &
                        +rat*0.5     &
                        *(bio_bdy_clim(next,:,:,o)+bio_bdy_clim(this,:,:,o))
      end do
   else

      if (first) then
         loop0=loop-1
      endif
      t = (loop-loop0)*dtm

      if(t .gt. t2 .or. first) then

         if (first) then
            first = .false.
            t2=t
         else
            call write_time_string()
            LEVEL3 timestr,': reading BIO boundary data ...'
         end if

         n = size(bdy_times)
         do i=1,n
            if(bdy_times(i) .ge. real(t + offset)) then
               EXIT
            end if
         end do
         start(1) = 1; edges(1) = zax_len;
         start(2) = 1; edges(2) = nsbv;
         start(3) = i; edges(3) = 1

         t1=t2
         t2 = bdy_times(i) - offset

         bio_old = bio_new

         do o=1,npel
            if (have_bio_bdy_values(o) .eq. 1 ) then
               err = nf90_get_var(ncid,bio_ids(o),bio_wrk(:,:,o),start,edges)
               if (err .ne. NF90_NOERR) go to 10
            end if
         end do

         l = 0
         do n=1,NWB
            l = l+1
            i = wi(n)
            do o=1,npel
               k = bdy_index(l)
               do j=wfj(n),wlj(n)
                  call interpol(zax_len,zlev,bio_wrk(:,k,o),H(i,j),kmax, &
                                hn(i,j,:),bio_new(:,k,o))
                  k = k+1
               end do
            end do
         end do

         do n = 1,NNB
            l = l+1
            j = nj(n)
            do o=1,npel
               k = bdy_index(l)
               do i = nfi(n),nli(n)
                  call interpol(zax_len,zlev,bio_wrk(:,k,o),H(i,j),kmax, &
                                hn(i,j,:),bio_new(:,k,o))
                  k = k+1
               end do
            end do
         end do

         do n=1,NEB
            l = l+1
            i = ei(n)
            do o=1,npel
               k = bdy_index(l)
               do j=efj(1),elj(1)
                  call interpol(zax_len,zlev,bio_wrk(:,k,o),H(i,j),kmax, &
                                hn(i,j,:),bio_new(:,k,o))
                  k = k+1
               end do
            end do
         end do

         do n = 1,NSB
            l = l+1
            j = sj(n)
            do o=1,npel
               k = bdy_index(l)
               do i = sfi(n),sli(n)
                  call interpol(zax_len,zlev,bio_wrk(:,k,o),H(i,j),kmax, &
                                hn(i,j,:),bio_new(:,k,o))
                  k = k+1
               end do
            end do
         end do
      end if

      bio_bdy = bio_old + (bio_new - bio_old)*(t-t1)/(t2-t1)

   end if

#ifdef DEBUG
   write(debug,*) 'Leaving do_3d_bio_bdy_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'do_3d_bio_bdy_ncdf: ',nf90_strerror(err)
   stop
   end subroutine do_3d_bio_bdy_ncdf
!EOC

!-----------------------------------------------------------------------

! quick and dirty - should be merged with kbk_interpol.F90 and
! grid_interpol.F90

   subroutine interpol(nlev,zlev,wrk,depth,kmax,zm,col)

! !INPUT PARAMETERS:
   integer, intent(in)       :: nlev
   REAL_4B, intent(in)       :: zlev(nlev),wrk(nlev)
   REALTYPE, intent(in)      :: depth
   integer, intent(in)       :: kmax
   REALTYPE, intent(in)      :: zm(0:kmax)

! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)     :: col(0:kmax)

! !LOCAL VARIABLES:

   REALTYPE                  :: zmodel(kmax),rat
   integer                   :: k,li,n,nn

   zmodel(1) = -depth + 0.5*zm(1)
   do k=2,kmax
      zmodel(k) = zmodel(k-1) + 0.5*(zm(k-1)+zm(k))
   end do

   do k=kmax,1,-1
      if (zmodel(k) .ge. zlev(1)) col(k) = wrk(1)
   end do

!  find largest index with valid value in wrk
   do li=1,nlev
      if (wrk(li) .lt. -999. ) EXIT
   end do
   if (li .ne. nlev .or. wrk(li) .lt. -999.) li=li-1

   do k=1,kmax
      if (zmodel(k) .le. zlev(li)) col(k) = wrk(li)
   end do

   do k=1,kmax
      if (zmodel(k) .gt. zlev(li) .and. zmodel(k) .lt. zlev(1)) then
         nn=nlev+1
224      nn=nn-1
         if(zlev(nn) .le. zmodel(k)) goto 224
         rat = (zmodel(k)-zlev(nn+1))/(zlev(nn)-zlev(nn+1))
         col(k) = (_ONE_-rat)*wrk(nn+1)+rat*wrk(nn)
      end if
   end do
   col(0)=col(1)
   end subroutine interpol
!-----------------------------------------------------------------------

   end module ncdf_3d_bio_bdy

!-----------------------------------------------------------------------
! Copyright (C) 2012 - Karsten Bolding and Jorn Bruggeman (BB)         !
!-----------------------------------------------------------------------

#endif
