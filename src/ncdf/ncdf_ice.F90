#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdf_ice -
!
! !INTERFACE:
   module ncdf_ice
!
! !DESCRIPTION:
!
! !USES:
   use netcdf
   use exceptions
   use time         ,only: string_to_julsecs,time_diff,add_secs,in_interval
   use time         ,only: jul0,secs0,julianday,secondsofday,timestep
   use time         ,only: write_time_string,timestr
   use domain       ,only: imin,imax,jmin,jmax,iextr,jextr
   use domain       ,only: ill,ihl,jll,jhl,ilg,ihg,jlg,jhg
   use domain       ,only: az,lonc,latc
   use grid_interpol,only: init_grid_interpol,do_grid_interpol
   use grid_interpol,only: to_rotated_lat_lon
   use getm_ice     ,only: ice_file,ice_hi
!   use variables_ice  ,only: ice_hi
   IMPLICIT NONE
!
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_ice_input_ncdf,get_ice_data_ncdf
!
! !PRIVATE DATA MEMBERS:
   integer         :: ncid
   integer         :: ice_hi_id=-1
   integer         :: ilen,jlen,nlen=0
   integer         :: start(3),edges(3)
   REALTYPE        :: offset
   logical         :: stationary,on_grid

   integer         :: grid_scan=1
   logical         :: point_source=.false.
   logical         :: rotated_ice_grid=.false.

   integer, allocatable      :: ice_hi_mask(:,:)
   REALTYPE, allocatable     :: ice_lon(:),ice_lat(:)

!  For gridinterpolation
   REALTYPE, allocatable     :: beta(:,:)
   REALTYPE, allocatable     :: ti(:,:),ui(:,:)
   integer, allocatable      :: gridmap(:,:,:)
!
   REALTYPE, parameter       :: pi=3.1415926535897932384626433832795029
   REALTYPE, parameter       :: deg2rad=pi/180.,rad2deg=180./pi
   REALTYPE                  :: southpole(3) = (/0.0,-90.0,0.0/)


   REALTYPE,dimension(:),allocatable   :: ice_times(:)
   REALTYPE,dimension(:,:),pointer     :: ice_hi_new,d_ice_hi,ice_hi_input
   REALTYPE,dimension(:,:),allocatable :: wrk

   character(len=10)         :: name_ice_hi="ice_hi"
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_ice_input_ncdf -
!
! !INTERFACE:
   subroutine init_ice_input_ncdf(fn,nstart)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn
   integer, intent(in)                 :: nstart
!
! !LOCAL VARIABLES:
   integer         :: il,ih,jl,jh
   integer         :: rc
!EOP
!-------------------------------------------------------------------------
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_ice_input_ncdf() # ',Ncall
#endif

   call open_ice_file(ice_file)

   if (ilen.eq.iextr .and. jlen.eq.jextr) then
      LEVEL3 'Assuming On-Grid ice forcing'
      on_grid = .true.
      il = ilg ; jl = jlg ; ih = ihg ; jh = jhg
   else if (ilen.eq.1 .and. jlen.eq.1) then
      LEVEL3 'Assuming Point Source ice forcing'
      point_source = .true.
      on_grid = .true.
      il = 1 ; jl = 1 ; ih = 1 ; jh = 1
   else
      on_grid = .false.
      il = 1 ; jl = 1 ; ih = ilen ; jh = jlen

      allocate(ti(E2DFIELD),stat=rc)
      if (rc /= 0) &
          stop 'init_meteo_input_ncdf: Error allocating memory (ti)'
      ti = -999.

      allocate(ui(E2DFIELD),stat=rc)
      if (rc /= 0) stop &
              'init_meteo_input_ncdf: Error allocating memory (ui)'
      ui = -999.

      allocate(gridmap(E2DFIELD,1:2),stat=rc)
      if (rc /= 0) stop &
              'init_meteo_input_ncdf: Error allocating memory (gridmap)'
      gridmap(:,:,:) = -999

      allocate(beta(E2DFIELD),stat=rc)
      if (rc /= 0) &
          stop 'init_meteo_input_ncdf: Error allocating memory (beta)'
      beta = _ZERO_

!     do not call with ice_hi_mask, otherwise we have nearest neighbour
      call init_grid_interpol(imin,imax,jmin,jmax,az,  &
                lonc,latc,ice_lon,ice_lat,southpole,gridmap,beta,ti,ui, &
                break_on_missing=.false.)

   end if

   start(1) = il; start(2) = jl;
   edges(1) = ih-il+1; edges(2) = jh-jl+1;
   edges(3) = 1

   allocate(wrk(edges(1),edges(2)),stat=rc)
   if (rc /= 0) call getm_error('init_ice_input_ncdf()',               &
                                'Error allocating memory (wrk)')

   allocate(ice_hi_new(E2DFIELD),stat=rc)
   if (rc /= 0) call getm_error('init_ice_input_ncdf()',               &
                                'Error allocating memory (ice_hi_new)')


   if ( stationary ) then

      ice_hi_input => ice_hi_new
      call read_data(0)
      ice_hi = ice_hi_new

   else

      allocate(d_ice_hi(E2DFIELD),stat=rc)
      if (rc /= 0) call getm_error('init_ice_input_ncdf()',            &
                                   'Error allocating memory (d_ice_hi)')

      ice_hi_input => d_ice_hi
      call get_ice_data_ncdf(nstart-1)

   end if

#ifdef DEBUG
   write(debug,*) 'Leaving init_ice_input_ncdf()'
   write(debug,*)
#endif
   return
   end subroutine init_ice_input_ncdf
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ice_data_ncdf - .
!
! !INTERFACE:
   subroutine get_ice_data_ncdf(loop)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: loop
!
! !LOCAL VARIABLES:
   integer                         :: indx
   integer, save                   :: save_n=1
   REALTYPE                        :: t,t_minus_t2
   REALTYPE,save                   :: t1,t2=-_ONE_
   REALTYPE,save                   :: deltm1=_ZERO_
   REALTYPE,dimension(:,:),pointer :: ice_hi_old
   logical, save                   :: first=.true.
!EOP
!-------------------------------------------------------------------------
#ifdef DEBUG
   integer, save   :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'get_ice_data_ncdf() # ',Ncall
#endif

   if (stationary) then
      ice_hi = ice_hi_new
      return
   end if

!  find the right index

   t = loop*timestep

   if ( t .gt. t2 ) then

      t1 = t2

!     Note (KK): Even if the first time stage is the last entry of an
!                input file, we must not jump to the next file because
!                then indx-1 does not work!
!                Therefore .ge. and save_n must be initialised to 1!
      do indx=save_n,nlen
         t2 = ice_times(indx) - offset
         if ( t2 .ge. t ) then
            EXIT
         end if
      end do

!     end of simulation?
      if (indx .gt. nlen) then
!        Note (KK): here we are not in case of first
!                   (because of in_interval check in open_ice_file)
         LEVEL2 'Need new ice file'
         call open_ice_file(ice_file)
         do indx=1,nlen
            t2 = ice_times(indx) - offset
!           Note (KK): For ice there is no check for long enough
!                      data sets. Therefore .ge.!
            if ( t2 .ge. t ) then
               EXIT
            end if
         end do
      end if

      if (first) then
         if ( t2 .gt. t ) then
            indx = indx-1
         end if
         t2 = ice_times(indx) - offset
      end if

      call read_data(indx)
      save_n = indx+1

      ice_hi_old=>ice_hi_new;ice_hi_new=>d_ice_hi;d_ice_hi=>ice_hi_old;ice_hi_input=>d_ice_hi

      if ( .not. first ) then
         d_ice_hi = ice_hi_new - ice_hi_old
         deltm1 = _ONE_ / (t2 - t1)
      end if

   end if


   t_minus_t2 = t - t2
   ice_hi = ice_hi_new + d_ice_hi*deltm1*t_minus_t2

   first = .false.

#ifdef DEBUG
   write(debug,*) 'Leaving get_ice_data_ncdf()'
   write(debug,*)
#endif
   return
   end subroutine get_ice_data_ncdf
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: open_ice_file - .
!
! !INTERFACE:
   subroutine open_ice_file(ice_file)
!
! !DESCRIPTION:
!  Instead of specifying the name of one single ncdf file directly - a list
!  of names can be specified in \emph{ice\_file}. The rationale for this
!  approach is that output from operational meteorological models are of
!  typically 2-5 days length. Collecting a number of these files allows for
!  longer model integrations without have to reformat the data.
!  It is assumed that the different files contains the same variables
!  and that they are of the same shape.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: ice_file
!
! !LOCAL VARIABLES:
   integer, parameter        :: iunit=62
   character(len=256)        :: fn,time_units
   integer         :: junit,sunit,j1,s1,j2,s2
   integer         :: n,err
   logical,save       :: first=.true.
   logical,save       :: first_open=.true.
   logical,save       :: found=.false.

   integer            :: ndims,nvardims
   integer            :: lon_dim,lat_dim,time_dim=-1
   integer            :: lon_id,lat_id,time_id=-1
   integer            :: dim_len(3),vardim_ids(3)
   character(len=16)  :: dim_name(3),name_ice_hi_mask
   integer            :: id
   logical            :: have_southpole
!
!EOP
!-------------------------------------------------------------------------
#ifdef DEBUG
   integer, save   :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'open_ice_file() # ',Ncall
#endif

   if (first) open(iunit,file=ice_file,status='old',action='read',err=80)

   found = .false.


   do

      if (.not. first_open) then
         err = nf90_close(ncid)
         if (err .NE. NF90_NOERR) go to 10
      end if

      read(iunit,*,err=85,end=90) fn
      LEVEL3 'Trying ice from:'
      LEVEL4 trim(fn)

      err = nf90_open(fn,NF90_NOWRITE,ncid)
      if (err .NE. NF90_NOERR) go to 10

      err = nf90_inquire(ncid, nDimensions = ndims)
      if (err .NE. NF90_NOERR) go to 10

      LEVEL4 'dimensions'
      do n=1,ndims
         err = nf90_inquire_dimension(ncid,n,name=dim_name(n),len=dim_len(n))
         if (err .NE. NF90_NOERR) go to 10
         LEVEL4 n,dim_name(n), dim_len(n)
      end do

      select case (ndims)
         case(2)
            LEVEL4 'stationary ice fields'
            stationary=.true.
            if (.not. first_open) call getm_error('open_ice_file()',   &
                        'stationary fields only possible in first file')
         case(3)
            LEVEL4 'non-stationary ice fields'
            stationary=.false.
         case default
           call getm_error('open_ice_file()', &
                           'invalid number of dimensions')
      end select

      LEVEL4 ' ... checking variable ',name_ice_hi
      err = nf90_inq_varid(ncid,name_ice_hi,ice_hi_id)
      if (err .NE. NF90_NOERR) go to 10
      err = nf90_inquire_variable(ncid,ice_hi_id,ndims=nvardims)
      if (err .NE. NF90_NOERR) go to 10
      if (nvardims .NE. ndims) call getm_error('open_ice_file()',      &
                                'Wrong number of dims in '//name_ice_hi)
      err = nf90_inquire_variable(ncid,ice_hi_id,dimids=vardim_ids)
      if (err .NE. NF90_NOERR) go to 10

      lon_dim = vardim_ids(1)
      lat_dim = vardim_ids(2)

      if (stationary) exit

      time_dim = vardim_ids(3)

      err = nf90_inq_varid(ncid,dim_name(time_dim),time_id)
      if (err .ne. NF90_NOERR) go to 10

      if (dim_len(time_dim) > nlen) then
         if (.not. first) then
            deallocate(ice_times,stat=err)
            if (err /= 0) call getm_error('open_ice_file()',           &
                               'Error de-allocating memory (ice_times)')
         end if
         allocate(ice_times(dim_len(time_dim)),stat=err)
         if (err /= 0) call getm_error('open_ice_file()',              &
                                  'Error allocating memory (ice_times)')
      end if
      nlen = dim_len(time_dim)
      err = nf90_get_var(ncid,time_id,ice_times(1:nlen))
      if (err .ne. NF90_NOERR) go to 10
      err =  nf90_get_att(ncid,time_id,'units',time_units)
      if (err .NE. NF90_NOERR) go to 10
      call string_to_julsecs(time_units,junit,sunit)

      offset = time_diff(jul0,secs0,junit,sunit)

      call add_secs(junit,sunit,nint(ice_times(1   )),j1,s1)
      call add_secs(junit,sunit,nint(ice_times(nlen)),j2,s2)

      if (first) then
         if (in_interval(j1,s1,julianday,secondsofday,j2,s2)) then
            found = .true.
         end if
      else
         if (time_diff(j2,s2,julianday,secondsofday) > _ZERO_) then
            found = .true.
         else
            LEVEL0 'WARNING: skipping ice file ',trim(fn)
         end if
      end if

      if (found) exit
      first_open = .false.

   end do

   if ( .not. found ) call getm_error('open_ice_file()',               &
                  'Could not find valid iceforcing in '//trim(ice_file))

   LEVEL3 'Using ice from:'
   LEVEL4 trim(fn)
   if ( .not. stationary ) then
      LEVEL3 'ice offset time ',offset
   end if
   ilen = dim_len(lon_dim)
   jlen = dim_len(lat_dim)

   if (first) then

      err = nf90_inq_varid(ncid,dim_name(lon_dim),lon_id)
      if (err .ne. NF90_NOERR) go to 10
      allocate(ice_lon(ilen),stat=err)
      if (err /= 0) call getm_error('open_ice_file()',                 &
                                    'Error allocating memory (ice_lon)')
      err = nf90_get_var(ncid,lon_id,ice_lon(1:ilen))
      if (err .ne. NF90_NOERR) go to 10

      err = nf90_inq_varid(ncid,dim_name(lat_dim),lat_id)
      if (err .ne. NF90_NOERR) go to 10
      allocate(ice_lat(ilen),stat=err)
      if (err /= 0) call getm_error('open_ice_file()',                 &
                                    'Error allocating memory (ice_lat)')
      err = nf90_get_var(ncid,lat_id,ice_lat(1:jlen))
      if (err .ne. NF90_NOERR) go to 10

!           first we check for CF compatible grid_mapping_name
            err = nf90_inq_varid(ncid,'rotated_pole',id)
            if (err .eq. NF90_NOERR) then
               LEVEL4 'Reading CF-compliant rotated grid specification'
               err = nf90_get_att(ncid,id, &
                                  'grid_north_pole_latitude',southpole(1))
               if (err .ne. NF90_NOERR) go to 10
               err = nf90_get_att(ncid,id, &
                                  'grid_north_pole_longitude',southpole(2))
               if (err .ne. NF90_NOERR) go to 10
               err = nf90_get_att(ncid,id, &
                                  'north_pole_grid_longitude',southpole(3))
               if (err .ne. NF90_NOERR) then
                  southpole(3) = _ZERO_
               end if
!              Northpole ---> Southpole transformation
               LEVEL4 'Transforming North Pole to South Pole specification'
               if (southpole(2) .ge. 0) then
                  southpole(2) = southpole(2) - 180.
               else
                  southpole(2) = southpole(2) + 180.
               end if 
               southpole(1) = -southpole(1)
               southpole(3) = _ZERO_
               have_southpole = .true.
               rotated_ice_grid = .true.
            else
               have_southpole = .false.
            end if

!           and then we revert to the old way - checking 'southpole' directly
            if (.not. have_southpole) then
               err = nf90_inq_varid(ncid,'southpole',id)
               if (err .ne. NF90_NOERR) then
                  LEVEL4 'Setting southpole to (0,-90,0)'
               else
                  err = nf90_get_var(ncid,id,southpole)
                  if (err .ne. NF90_NOERR) go to 10
                  rotated_ice_grid = .true.
               end if
            end if
            if (rotated_ice_grid) then
               LEVEL4 'south pole:'
!              changed indices - kb 2014-12-15
               LEVEL4 '      lon ',southpole(2)
               LEVEL4 '      lat ',southpole(1)
            end if

      allocate(ice_hi_mask(1:ilen,1:jlen),stat=err)
      if (err /= 0) &
         stop 'init_meteo_input_ncdf: Error allocating memory (ice_hi_mask)'
      err =  nf90_get_att(ncid,ice_hi_id,'mask',name_ice_hi_mask)
      if (err .ne. NF90_NOERR) name_ice_hi_mask='ice_mask'
      err = nf90_inq_varid(ncid,trim(name_ice_hi_mask),id)
      if (err .eq. NF90_NOERR) then
         LEVEL4 'taking variable' // trim(name_ice_hi_mask) // ' as ice_hi_mask'
         err = nf90_get_var(ncid,id,ice_hi_mask)
         if (err .ne. NF90_NOERR) go to 10
      else
         ice_hi_mask = 1
      end if

   end if


   first = .false.

#ifdef DEBUG
   write(debug,*) 'Leaving open_ice_file()'
   write(debug,*)
#endif
   return

10 FATAL 'open_ice_file: ',nf90_strerror(err)
   stop 'open_ice_file()'
80 FATAL 'I could not open: ',trim(ice_file)
   stop 'open_ice_file()'
85 FATAL 'Error reading: ',trim(ice_file)
   stop 'open_ice_file()'
90 FATAL 'Reached eof in: ',trim(ice_file)
   stop 'open_ice_file()'

   end subroutine open_ice_file
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_data -
!
! !INTERFACE:
   subroutine read_data(indx)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in) :: indx
!
! !LOCAL VARIABLES:
   integer                      :: err
!EOP
!-----------------------------------------------------------------------

   call write_time_string()
   LEVEL3 timestr,': reading ice data ...',indx
   start(3) = indx

   err = nf90_get_var(ncid,ice_hi_id,wrk,start,edges)
   if (err .ne. NF90_NOERR) go to 10
   if (on_grid) then
      if (point_source) then
         ice_hi_input = wrk(1,1)
      else
         ice_hi_input(ill:ihl,jll:jhl) = wrk
      end if
   else
      call do_grid_interpol(az,wrk,gridmap,ti,ui,ice_hi_input,         &
                            imask=ice_hi_mask,fillvalue=_ZERO_)
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving read_data()'
   write(debug,*)
#endif
   return

10 FATAL 'read_data: ',nf90_strerror(err)
   stop

   end subroutine read_data
!EOC
!-----------------------------------------------------------------------

   end module ncdf_ice

!-----------------------------------------------------------------------
! Copyright (C) 2014 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
