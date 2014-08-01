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
   use netcdf
   use domain, only: imin,imax,jmin,jmax,kmax,ioff,joff
   use domain, only: nsbv,nsbvl,nbdy,NWB,NNB,NEB,NSB
   use domain, only: bdy_index,bdy_index_l,bdy_index_stop
   use domain, only: wi,wfj,wlj,nj,nfi,nli,ei,efj,elj,sj,sfi,sli
   use variables_2d, only: dtm
   use m3d, only: update_salt,update_temp
   use bdy_3d, only: bdy3d_vel
   use bdy_3d, only: bdy_data_uu,bdy_data_vv,bdy_data_S,bdy_data_T
   use time, only: string_to_julsecs,time_diff,add_secs
   use time, only: julianday,secondsofday,juln,secsn
   use time, only: write_time_string,timestr
   IMPLICIT NONE
!
   private
!
   public                              :: init_3d_bdy_ncdf,do_3d_bdy_ncdf
!
! !PRIVATE DATA MEMBERS:
   integer                             :: ncid
   integer                             :: time_id,temp_id=-1,salt_id=-1
   integer                             :: uu_id=-1,vv_id=-1
   integer                             :: bdy_dim,bdy_len,bdy_pos
   integer                             :: zax_dim=-1,zax_len,zax_pos
   integer                             :: time_dim=-1,time_len,time_pos
   logical                             :: climatology=.false.
   logical                             :: from_3d_fields
   REALTYPE,dimension(:),allocatable   :: zlev
   REALTYPE,dimension(:,:),pointer     :: uu_bdy,uu_bdy_new,d_uu_bdy
   REALTYPE,dimension(:,:),pointer     :: vv_bdy,vv_bdy_new,d_vv_bdy
   REALTYPE,dimension(:,:),pointer     :: S_bdy,S_bdy_new,d_S_bdy
   REALTYPE,dimension(:,:),pointer     :: T_bdy,T_bdy_new,d_T_bdy
!  the following is used for climatology=.true.
   REALTYPE,dimension(:,:,:),allocatable :: uu_bdy_clim,vv_bdy_clim
   REALTYPE,dimension(:,:,:),allocatable :: S_bdy_clim,T_bdy_clim
!  the following is used for climatology=.false.
   integer                             :: loop0
   REALTYPE                            :: offset
   REALTYPE,dimension(:),allocatable   :: bdy_times
   REALTYPE,dimension(:,:),allocatable :: wrk
   integer                             :: wrk_len
   integer,parameter                   :: climatology_len=12
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
! !ROUTINE: init_3d_bdy_ncdf -
!
! !INTERFACE:
   subroutine init_3d_bdy_ncdf(fname,loop)
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
   character(len=256)        :: units
   character(len=19)         :: tbuf
   integer                   :: j1,s1,j2,s2
   integer                   :: ndims, nvardims
   integer                   :: vardim_ids(4)
   integer, allocatable, dimension(:):: dim_ids,dim_len
   character(len=16), allocatable :: dim_name(:)
   integer                   :: rc,err
   integer                   :: n,id
   logical                   :: dims_ready
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   write(debug,*) 'ncdf_init_3d_bdy (NetCDF)'
   write(debug,*) 'Reading from: ',trim(fname)
#endif

   LEVEL3 'init_3d_bdy_ncdf'

   err = nf90_open(fname,NF90_NOWRITE,ncid)
   if (err .NE. NF90_NOERR) go to 10

   err = nf90_inquire(ncid, nDimensions = nDims)
   if (err .NE. NF90_NOERR) go to 10

   allocate(dim_ids(ndims),stat=rc)
   if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (dim_ids)'
   allocate(dim_len(ndims),stat=rc)
   if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (dim_len)'
   allocate(dim_name(ndims),stat=rc)
   if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (dim_name)'

   do n=1,ndims
      err = nf90_inquire_dimension(ncid,n,name=dim_name(n),len=dim_len(n))
      if (err .NE. NF90_NOERR) go to 10
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
      zax_pos = 3
      time_pos = 4
   else
!     We are reading boundary values from a special boundary data file
!     The variables 'salt' and 'temp' must both exist and be spanned by
!     dimensions as:
!       1 -> zax,levels
!       2 -> bdy_points
!       3 -> time
      LEVEL4 'special boundary data file'
      from_3d_fields=.false.
      zax_pos = 1
      bdy_pos = 2
      time_pos = 3

!     Note(BJB): This test may break backward compatibility,
!                so I leave it out for now:
      !if (ndims .NE. 3) stop 'init_3d_bdy_ncdf: Wrong number of dims in file (must be 3)'
   end if

!  We will use this information to actually find the dimension
!  index numbers in the data set.
!  Some of the tests will be repeated later (fixing is possible but not
!  high priority, BJB 2007-04-25).

   dims_ready = .false.
   if (bdy3d_vel) then
      LEVEL4 ' ... checking variable "uu"'
      err = nf90_inq_varid(ncid,'uu',uu_id)
      if (err .NE. NF90_NOERR) go to 10
      err = nf90_inquire_variable(ncid,uu_id,ndims=nvardims)
      if (err .NE. NF90_NOERR) go to 10
      if (nvardims .NE. ndims) &
           stop 'init_3d_bdy_ncdf: Wrong number of dims in uu'
      err = nf90_inquire_variable(ncid,uu_id,dimids=vardim_ids)
      if (err .NE. NF90_NOERR) go to 10
      zax_dim  = vardim_ids(zax_pos)
      time_dim = vardim_ids(time_pos)
      dims_ready = .true.
      LEVEL4 ' ... checking variable "vv"'
      err = nf90_inq_varid(ncid,'vv',vv_id)
      if (err .NE. NF90_NOERR) go to 10
      err = nf90_inquire_variable(ncid,vv_id,ndims=nvardims)
      if (err .NE. NF90_NOERR) go to 10
      if (nvardims .NE. ndims) &
           stop 'init_3d_bdy_ncdf: Wrong number of dims in vv'
      err = nf90_inquire_variable(ncid,vv_id,dimids=vardim_ids)
      if (err .NE. NF90_NOERR) go to 10
      if (zax_dim /= vardim_ids(zax_pos)) &
           stop 'init_3d_bdy_ncdf: Position of zax dimension of vv'
      if (time_dim /= vardim_ids(time_pos)) &
           stop 'init_3d_bdy_ncdf: Position of time dimension of vv'
   end if
   if (update_salt) then
      LEVEL4 ' ... checking variable "salt"'
      err = nf90_inq_varid(ncid,'salt',salt_id)
      if (err .NE. NF90_NOERR) go to 10
      err = nf90_inquire_variable(ncid,salt_id,ndims=nvardims)
      if (err .NE. NF90_NOERR) go to 10
      if (nvardims .NE. ndims) &
           stop 'init_3d_bdy_ncdf: Wrong number of dims in salt'
      err = nf90_inquire_variable(ncid,salt_id,dimids=vardim_ids)
      if (err .NE. NF90_NOERR) go to 10
      if (dims_ready) then
         if (zax_dim /= vardim_ids(zax_pos)) &
              stop 'init_3d_bdy_ncdf: Position of zax dimension of salt'
         if (time_dim /= vardim_ids(time_pos)) &
              stop 'init_3d_bdy_ncdf: Position of time dimension of salt'
      else
         zax_dim  = vardim_ids(zax_pos)
         time_dim = vardim_ids(time_pos)
         dims_ready = .true.
      end if
   end if
   if (update_temp) then
      LEVEL4 ' ... checking variable "temp"'
      err = nf90_inq_varid(ncid,'temp',temp_id)
      if (err .NE. NF90_NOERR) go to 10
      err = nf90_inquire_variable(ncid,temp_id,ndims=nvardims)
      if (err .NE. NF90_NOERR) go to 10
      if (nvardims .NE. ndims) &
           stop 'init_3d_bdy_ncdf: Wrong number of dims in temp'
      err = nf90_inquire_variable(ncid,temp_id,dimids=vardim_ids)
      if (err .NE. NF90_NOERR) go to 10
      if (dims_ready) then
         if (zax_dim /= vardim_ids(zax_pos)) &
              stop 'init_3d_bdy_ncdf: Position of zax dimension of temp'
         if (time_dim /= vardim_ids(time_pos)) &
              stop 'init_3d_bdy_ncdf: Position of time dimension of temp'
      else
         zax_dim  = vardim_ids(zax_pos)
         time_dim = vardim_ids(time_pos)
         dims_ready = .true.
      end if
   end if

   if (.not. from_3d_fields) then
      bdy_dim = vardim_ids(bdy_pos)
      bdy_len = dim_len(bdy_dim)
      if (bdy_len .lt. nsbv) then
         stop 'init_3d_bdy_ncdf: netcdf file does not contain enough bdy points'
      else if (bdy_len .gt. nsbv) then
         LEVEL4 'WARNING: netcdf file contains data for more bdy points'
      end if
   end if

   zax_len = dim_len(zax_dim)
   time_len = dim_len(time_dim)

   allocate(zlev(zax_len),stat=rc)
   if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (zlev)'
   zlev = _ZERO_

   err = nf90_inq_varid(ncid, dim_name(zax_dim), id)
   if (err .ne. NF90_NOERR) go to 10
   err = nf90_get_var(ncid,id,zlev)
   if (err .ne. NF90_NOERR) go to 10

!  a few sanity checks on the vertical axis for the 3D boundaries
   do n=1,zax_len
      if (zlev(n) .eq. NF90_FILL_REAL) then
         FATAL '3D boundary z-axis contains NF90_FILL_REAL values'
         FATAL 'proper interpolation cant be done'
         stop 'init_3d_bdy_ncdf'
      end if
   end do
!  not sure if this check is safe - kb
   if ( zlev(1) .ge. _ZERO_ .and. zlev(zax_len) .gt. _ZERO_ ) then
      LEVEL4 'converting positive z-axis (depth) values to negative'
      zlev = -_ONE_*zlev
   end if
!  check strict monotonicity
   do n=1,zax_len-1
      if ( .not. zlev(n) .gt. zlev(n+1) ) then
         FATAL '3D boundary z-axis not strict monotone: ',zlev(n),zlev(n+1)
         stop 'init_3d_bdy_ncdf'
      end if
   end do

   if (bdy3d_vel) then
      allocate(uu_bdy(zax_len,nsbvl),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (uu_bdy)'
      allocate(uu_bdy_new(zax_len,nsbvl),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (uu_bdy_new)'
      allocate(d_uu_bdy(zax_len,nsbvl),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (d_uu_bdy)'
      allocate(vv_bdy(zax_len,nsbvl),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (vv_bdy)'
      allocate(vv_bdy_new(zax_len,nsbvl),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (vv_bdy_new)'
      allocate(d_vv_bdy(zax_len,nsbvl),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (d_vv_bdy)'
   end if
   if (update_salt) then
      allocate(S_bdy(zax_len,nsbvl),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (S_bdy)'
      allocate(S_bdy_new(zax_len,nsbvl),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (S_bdy_new)'
      allocate(d_S_bdy(zax_len,nsbvl),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (d_S_bdy)'
   end if
   if (update_temp) then
      allocate(T_bdy(zax_len,nsbvl),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (T_bdy)'
      allocate(T_bdy_new(zax_len,nsbvl),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (T_bdy_new)'
      allocate(d_T_bdy(zax_len,nsbvl),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (d_T_bdy)'
   end if

   if( time_len .eq. climatology_len) then
      climatology=.true.
      LEVEL4 'Assuming climatolgical 3D boundary conditions'
      LEVEL4 '# of times = ',time_len
   end if

   if (climatology) then

      if (bdy3d_vel) then
         allocate(uu_bdy_clim(climatology_len,zax_len,nsbvl),stat=rc)
         if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (uu_bdy_clim)'
         allocate(vv_bdy_clim(climatology_len,zax_len,nsbvl),stat=rc)
         if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (vv_bdy_clim)'
      end if
      if (update_salt) then
         allocate(S_bdy_clim(climatology_len,zax_len,nsbvl),stat=rc)
         if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (S_bdy_clim)'
      end if
      if (update_temp) then
         allocate(T_bdy_clim(climatology_len,zax_len,nsbvl),stat=rc)
         if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (T_bdy_clim)'
      end if

      call read_3d_bdy_climatology_ncdf()

      err = nf90_close(ncid)

   else

      if (from_3d_fields) then
         FATAL 'non-climatology bdy data only support special bdy data file'
         stop 'init_3d_bdy_ncdf'
      end if

      err = nf90_inq_varid(ncid,'time',time_id)
      if (err .NE. NF90_NOERR) go to 10
      err =  nf90_get_att(ncid,time_id,'units',units)
      if (err .NE. NF90_NOERR) go to 10

      allocate(bdy_times(time_len),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (bdy_times)'
      bdy_times = _ZERO_
      err = nf90_get_var(ncid,time_id,bdy_times)
      if (err .NE. NF90_NOERR) go to 10

      loop0 = loop - 1
      call string_to_julsecs(units,j1,s1)
      offset = time_diff(julianday,secondsofday,j1,s1)
      if( offset .lt. bdy_times(1) ) then
         FATAL 'Model simulation starts before available boundary data'
         call write_time_string(julianday,secondsofday,tbuf)
         FATAL 'Simulation starts: ',tbuf
         call add_secs(j1,s1,nint(bdy_times(1)),j2,s2)
         call write_time_string(j2,s2,tbuf)
         FATAL 'Datafile starts:   ',tbuf
         stop 'init_3d_bdy_ncdf'
      else
         LEVEL3 'Boundary offset time ',offset
      end if

!     check if the bdy data file is long enough
      if( time_diff(juln,secsn,j1,s1) .gt. bdy_times(time_len) ) then
         FATAL 'Not enough 3D boundary data in file'
         call write_time_string(juln,secsn,tbuf)
         FATAL 'Simulation ends: ',tbuf
         call add_secs(j1,s1,nint(bdy_times(time_len)),j2,s2)
         call write_time_string(j2,s2,tbuf)
         FATAL 'Datafile ends:   ',tbuf
         stop 'init_3d_bdy_ncdf'
      end if

      wrk_len = bdy_index_stop - bdy_index(1) + 1
      allocate(wrk(zax_len,bdy_index(1):bdy_index_stop),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (wrk)'
      wrk = _ZERO_

      call do_3d_bdy_ncdf(loop0)

   end if


#ifdef DEBUG
   write(debug,*) 'Leaving init_3d_bdy_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'init_3d_bdy_ncdf: ',nf90_strerror(err)
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
   use domain, only: H,HU,HV
   use variables_3d, only: hn,hun,hvn
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: loop
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer,save    :: this,indx=1,start(3),edges(3)
   integer         :: prev,next,i,err
   logical, save   :: first=.true.
   logical         :: new_set
   REALTYPE        :: t,t_minus_t2
   REALTYPE, save  :: t1=_ZERO_,t2=-_ONE_,deltm1
   REALTYPE,dimension(:,:),pointer :: uu_bdy_old,vv_bdy_old,S_bdy_old,T_bdy_old
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   write(debug,*) 'do_3d_bdy_ncdf (NetCDF)'
#endif

   new_set = .false.


   if ( climatology ) then


!     Note(KK): We already read in all data and only need to interpolate in time

      if (first) then
!        needed to initialise prev
         this = month-1
         if (this .eq. 0) this=time_len
      end if

      t = (day-1)*secsprday+secondsofday

      if (month .ne. this) then

         new_set = .true.

         prev = this
         this = month
         next = mod(month,climatology_len)+1

!        t1 initialised to 0
         t2 = secsprday*days_in_mon(leapyear,month)

!        Note (KK): in principle the average can be calculated only once
         if (first) then
            if (bdy3d_vel) then
               uu_bdy_new = _HALF_ * ( uu_bdy_clim(prev,:,:) + uu_bdy_clim(this,:,:) )
               vv_bdy_new = _HALF_ * ( vv_bdy_clim(prev,:,:) + vv_bdy_clim(this,:,:) )
            end if
            if (update_salt) then
               S_bdy_new = _HALF_ * ( S_bdy_clim(prev,:,:) + S_bdy_clim(this,:,:) )
            end if
            if (update_temp) then
               T_bdy_new = _HALF_ * ( T_bdy_clim(prev,:,:) + T_bdy_clim(this,:,:) )
            end if
            first = .false.
         end if

         if (bdy3d_vel) then
            uu_bdy = _HALF_ * ( uu_bdy_clim(this,:,:) + uu_bdy_clim(next,:,:) )
            vv_bdy = _HALF_ * ( vv_bdy_clim(this,:,:) + vv_bdy_clim(next,:,:) )
         end if
         if (update_salt) then
            S_bdy = _HALF_ * ( S_bdy_clim(this,:,:) + S_bdy_clim(next,:,:) )
         end if
         if (update_temp) then
            T_bdy = _HALF_ * ( T_bdy_clim(this,:,:) + T_bdy_clim(next,:,:) )
         end if

      end if


   else


      t = (loop-loop0)*dtm

      if(t .gt. t2) then

         new_set = .true.

         t1 = t2
         do i=indx+1,time_len
            t2 = bdy_times(i) - offset
            if(t2 .gt. t) then
               EXIT
            end if
         end do

         if (first) then
            indx = i-1
            t2 = bdy_times(indx) - offset
            start(1) = 1
            edges(1) = zax_len
            start(2) = bdy_index(1)
            edges(2) = wrk_len
            edges(3) = 1
            first = .false.
         else
            indx = i
         end if
         start(3) = indx

         call write_time_string()
         LEVEL3 timestr,': reading 3D boundary data ...',indx

!        Note(KK): We read in at once the data of all global bdy cells
!                  but only for the current time stage.
!                  Interpolation extracts all local bdy cells.

         if (uu_id .ne. -1) then
            err = nf90_get_var(ncid,uu_id,wrk,start,edges)
            if (err .ne. NF90_NOERR) go to 10
            call grid_3d_bdy_data_ncdf(wrk,uu_bdy)
         end if
         if (vv_id .ne. -1) then
            err = nf90_get_var(ncid,vv_id,wrk,start,edges)
            if (err .ne. NF90_NOERR) go to 10
            call grid_3d_bdy_data_ncdf(wrk,vv_bdy)
         end if
         if (salt_id .ne. -1) then
            err = nf90_get_var(ncid,salt_id,wrk,start,edges)
            if (err .ne. NF90_NOERR) go to 10
            call grid_3d_bdy_data_ncdf(wrk,S_bdy)
         end if
         if (temp_id .ne. -1) then
            err = nf90_get_var(ncid,temp_id,wrk,start,edges)
            if (err .ne. NF90_NOERR) go to 10
            call grid_3d_bdy_data_ncdf(wrk,T_bdy)
         end if

      end if


   end if


   if (new_set) then
      if (bdy3d_vel) then
         uu_bdy_old=>uu_bdy_new;uu_bdy_new=>uu_bdy;uu_bdy=>d_uu_bdy;d_uu_bdy=>uu_bdy_old
         d_uu_bdy = uu_bdy_new - uu_bdy_old
         vv_bdy_old=>vv_bdy_new;vv_bdy_new=>vv_bdy;vv_bdy=>d_vv_bdy;d_vv_bdy=>vv_bdy_old
         d_vv_bdy = vv_bdy_new - vv_bdy_old
      end if
      if (update_salt) then
         S_bdy_old=>S_bdy_new;S_bdy_new=>S_bdy;S_bdy=>d_S_bdy;d_S_bdy=>S_bdy_old
         d_S_bdy = S_bdy_new - S_bdy_old
      end if
      if (update_temp) then
         T_bdy_old=>T_bdy_new;T_bdy_new=>T_bdy;T_bdy=>d_T_bdy;d_T_bdy=>T_bdy_old
         d_T_bdy = T_bdy_new - T_bdy_old
      end if
      deltm1 = _ONE_ / (t2 - t1)
   end if

   t_minus_t2 = t - t2

   if (bdy3d_vel) then
      uu_bdy = uu_bdy_new + d_uu_bdy*deltm1*t_minus_t2
      call interpolate_3d_bdy_ncdf(nsbvl,zax_len,uu_bdy,HU,kmax,hun,bdy_data_uu)
      vv_bdy = vv_bdy_new + d_vv_bdy*deltm1*t_minus_t2
      call interpolate_3d_bdy_ncdf(nsbvl,zax_len,vv_bdy,HV,kmax,hvn,bdy_data_vv)
   end if
   if (update_salt) then
      S_bdy = S_bdy_new + d_S_bdy*deltm1*t_minus_t2
      call interpolate_3d_bdy_ncdf(nsbvl,zax_len,S_bdy,H,kmax,hn,bdy_data_S)
   end if
   if (update_temp) then
      T_bdy = T_bdy_new + d_T_bdy*deltm1*t_minus_t2
      call interpolate_3d_bdy_ncdf(nsbvl,zax_len,T_bdy,H,kmax,hn,bdy_data_T)
   end if


#ifdef DEBUG
   write(debug,*) 'Leaving do_3d_bdy_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'do_3d_bdy_ncdf: ',nf90_strerror(err)
   stop
   end subroutine do_3d_bdy_ncdf
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: grid_3d_bdy_data_ncdf -
!
! !INTERFACE:
   subroutine grid_3d_bdy_data_ncdf(wrk,col)
!
! !DESCRIPTION:
!  kurt,kurt
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in) :: wrk(zax_len,bdy_index(1):bdy_index_stop)
!
! !OUTPUT PARAMETERS:
   REALTYPE,intent(out) :: col(zax_len,nsbvl)
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
   write(debug,*) 'grid_3d_bdy_data_ncdf'
#endif

   l = 0
   do n=1,NWB
      l = l+1
      k = bdy_index(l)
      kl = bdy_index_l(l)
      i = wi(n)
      do j=wfj(n),wlj(n)
         col(:,kl) = wrk(:,k)
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
         col(:,kl) = wrk(:,k)
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
         col(:,kl) = wrk(:,k)
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
         col(:,kl) = wrk(:,k)
         k = k+1
         kl = kl + 1
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving grid_3d_bdy_data_ncdf()'
   write(debug,*)
#endif
   return
   end subroutine grid_3d_bdy_data_ncdf
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: read_3d_bdy_climatology_ncdf -
!
! !INTERFACE:
   subroutine read_3d_bdy_climatology_ncdf()
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer         :: start(4),edges(4)
   integer         :: i,j,k,kl,l,m,n,err
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   write(debug,*) 'read_3d_bdy_climatology_ncdf'
#endif

!  Note(KK): We read in the data columnwise for all time stages
!  here we can read from both a 3D field and from a
!  special boundary data file - only the arguments 'start' and 'edges'
!  varies in the calls to 'nf90_get_var()'
!  m counts the time
!  l counts the boundary number
!  k counts the global index of the bdy point
!  kl counts the local index of the bdy point
!  MUST cover the same area as in topo.nc

   edges = 1
   edges(zax_pos) = zax_len
   start(zax_pos) = 1

   do m=1,climatology_len
      start(time_pos) = m
      l = 0
      do n=1,NWB
         l = l+1
         k = bdy_index(l)
         kl = bdy_index_l(l)
         i = wi(n)
         do j=wfj(n),wlj(n)
            if (from_3d_fields) then
               start(1) = i+ioff ; start(2) = j+joff
            else
               start(2) = k
            end if
            if (uu_id .ne. -1) then
               err = nf90_get_var(ncid,uu_id,uu_bdy_clim(m,:,kl),start,edges)
               if (err .ne. NF90_NOERR) go to 10
            end if
            if (vv_id .ne. -1) then
               err = nf90_get_var(ncid,vv_id,vv_bdy_clim(m,:,kl),start,edges)
               if (err .ne. NF90_NOERR) go to 10
            end if
            if (salt_id .ne. -1) then
               err = nf90_get_var(ncid,salt_id,S_bdy_clim(m,:,kl),start,edges)
               if (err .ne. NF90_NOERR) go to 10
            end if
            if (temp_id .ne. -1) then
               err = nf90_get_var(ncid,temp_id,T_bdy_clim(m,:,kl),start,edges)
               if (err .ne. NF90_NOERR) go to 10
            end if
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
            if (from_3d_fields) then
               start(1) = i+ioff ; start(2) = j+joff
            else
               start(2) = k
            end if
            if (uu_id .ne. -1) then
               err = nf90_get_var(ncid,uu_id,uu_bdy_clim(m,:,kl),start,edges)
               if (err .ne. NF90_NOERR) go to 10
            end if
            if (vv_id .ne. -1) then
               err = nf90_get_var(ncid,vv_id,vv_bdy_clim(m,:,kl),start,edges)
               if (err .ne. NF90_NOERR) go to 10
            end if
            if (salt_id .ne. -1) then
               err = nf90_get_var(ncid,salt_id,S_bdy_clim(m,:,kl),start,edges)
               if (err .ne. NF90_NOERR) go to 10
            end if
            if (temp_id .ne. -1) then
               err = nf90_get_var(ncid,temp_id,T_bdy_clim(m,:,kl),start,edges)
               if (err .ne. NF90_NOERR) go to 10
            end if
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
            if (from_3d_fields) then
               start(1) = i+ioff ; start(2) = j+joff
            else
               start(2) = k
            end if
            if (uu_id .ne. -1) then
               err = nf90_get_var(ncid,uu_id,uu_bdy_clim(m,:,kl),start,edges)
               if (err .ne. NF90_NOERR) go to 10
            end if
            if (vv_id .ne. -1) then
               err = nf90_get_var(ncid,vv_id,vv_bdy_clim(m,:,kl),start,edges)
               if (err .ne. NF90_NOERR) go to 10
            end if
            if (salt_id .ne. -1) then
               err = nf90_get_var(ncid,salt_id,S_bdy_clim(m,:,kl),start,edges)
               if (err .ne. NF90_NOERR) go to 10
            end if
            if (temp_id .ne. -1) then
               err = nf90_get_var(ncid,temp_id,T_bdy_clim(m,:,kl),start,edges)
               if (err .ne. NF90_NOERR) go to 10
            end if
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
            if (from_3d_fields) then
               start(1) = i+ioff ; start(2) = j+joff
            else
               start(2) = k
            end if
            if (uu_id .ne. -1) then
               err = nf90_get_var(ncid,uu_id,uu_bdy_clim(m,:,kl),start,edges)
               if (err .ne. NF90_NOERR) go to 10
            end if
            if (vv_id .ne. -1) then
               err = nf90_get_var(ncid,vv_id,vv_bdy_clim(m,:,kl),start,edges)
               if (err .ne. NF90_NOERR) go to 10
            end if
            if (salt_id .ne. -1) then
               err = nf90_get_var(ncid,salt_id,S_bdy_clim(m,:,kl),start,edges)
               if (err .ne. NF90_NOERR) go to 10
            end if
            if (temp_id .ne. -1) then
               err = nf90_get_var(ncid,temp_id,T_bdy_clim(m,:,kl),start,edges)
               if (err .ne. NF90_NOERR) go to 10
            end if
            k = k+1
            kl = kl + 1
         end do
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving read_3d_bdy_climatology_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'read_3d_bdy_climatology_ncdf: ',nf90_strerror(err)
   stop
   end subroutine read_3d_bdy_climatology_ncdf
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: interpolate_3d_bdy_ncdf -
!
! !INTERFACE:
   subroutine interpolate_3d_bdy_ncdf(nsbvl,nlev,data_zax,H,kmax,hn,data_gvc)
!
! !DESCRIPTION:
!  Here the interpolation is called for the locally active bdy columns.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)   :: nsbvl,nlev,kmax
   REALTYPE,intent(in)  :: data_zax(nlev,nsbvl)
   REALTYPE,dimension(E2DFIELD),intent(in) :: H
   REALTYPE,dimension(I3DFIELD),intent(in) :: hn

! !OUTPUT PARAMETERS:
   REALTYPE,intent(out) :: data_gvc(0:kmax,nsbvl)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer         :: i,j,k,kl,l,n
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   write(debug,*) 'interpolate_3d_bdy_ncdf'
#endif

   l = 0
   do n=1,NWB
      l = l+1
      k = bdy_index(l)
      kl = bdy_index_l(l)
      i = wi(n)
      do j=wfj(n),wlj(n)
         call interpol(nlev,zlev,data_zax(:,kl),H(i,j),kmax,hn(i,j,:), &
                       data_gvc(:,kl))
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
         call interpol(nlev,zlev,data_zax(:,kl),H(i,j),kmax,hn(i,j,:), &
                       data_gvc(:,kl))
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
         call interpol(nlev,zlev,data_zax(:,kl),H(i,j),kmax,hn(i,j,:), &
                       data_gvc(:,kl))
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
         call interpol(nlev,zlev,data_zax(:,kl),H(i,j),kmax,hn(i,j,:), &
                       data_gvc(:,kl))
         k = k+1
         kl = kl + 1
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving interpolate_3d_bdy_ncdf()'
   write(debug,*)
#endif
   return
   end subroutine interpolate_3d_bdy_ncdf
!EOC
!-----------------------------------------------------------------------

! quick and dirty - should be merged with kbk_interpol.F90 and
! grid_interpol.F90

   subroutine interpol(nlev,zlev,wrk,depth,kmax,zm,col)

! !INPUT PARAMETERS:
   integer,intent(in)                    :: nlev,kmax
   REALTYPE,dimension(nlev),intent(in)   :: zlev,wrk
   REALTYPE,intent(in)                   :: depth
   REALTYPE,dimension(0:kmax),intent(in) :: zm

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
   ! BJB-NOTE: Typically, li will end up as nlev+1, so the first
   !   of the following tests gets false. However, during debug
   !   compilation the second condition *MAY* evaulate wrk(li),
   !   which will result in a "forrtl: severe".
   !if (li .ne. nlev .or. wrk(li) .lt. -999.) li=li-1
   if (li .ne. nlev) then
      li=li-1
   elseif (wrk(li) .lt. -999.) then
      li=li-1
   end if

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

   end module ncdf_3d_bdy

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
