#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdf_waves -
!
! !INTERFACE:
   module ncdf_waves
!
! !DESCRIPTION:
!
! !USES:
   use netcdf
   use exceptions
   use time           ,only: string_to_julsecs,time_diff,add_secs,in_interval
   use time           ,only: jul0,secs0,julianday,secondsofday,timestep
   use time           ,only: write_time_string,timestr
   use domain         ,only: imin,imax,jmin,jmax,iextr,jextr
   use domain         ,only: ill,ihl,jll,jhl,ilg,ihg,jlg,jhg
   use domain         ,only: az,lonc,latc,cosconv,sinconv
   use grid_interpol  ,only: init_grid_interpol,do_grid_interpol
   use grid_interpol  ,only: to_rotated_lat_lon
   use waves          ,only: waves_file
   use variables_waves,only: waveH,waveL,coswavedir,sinwavedir
   IMPLICIT NONE
!
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_waves_input_ncdf,get_waves_data_ncdf
!
! !PRIVATE DATA MEMBERS:
   integer         :: ncid
   integer         :: waveH_id=-1
   integer         :: waveL_id=-1
   integer         :: waveDir_id=-1
   integer         :: ilen,jlen,nlen=0
   integer         :: start(3),edges(3)
   REALTYPE        :: offset
   logical         :: stationary,on_grid

   integer         :: grid_scan=1
   logical         :: point_source=.false.
   logical         :: rotated_waves_grid=.false.

   integer, allocatable      :: waves_mask(:,:)
   REALTYPE, allocatable     :: waves_lon(:),waves_lat(:)

!  For gridinterpolation
   REALTYPE, allocatable     :: beta(:,:)
   REALTYPE, allocatable     :: ti(:,:),ui(:,:)
   integer, allocatable      :: gridmap(:,:,:)
!
   REALTYPE, parameter       :: pi=3.1415926535897932384626433832795029
   REALTYPE, parameter       :: deg2rad=pi/180.,rad2deg=180./pi
   REALTYPE                  :: southpole(3) = (/0.0,-90.0,0.0/)

   REALTYPE,dimension(:),allocatable   :: wave_times(:)
   REALTYPE,dimension(:,:),pointer     :: waveH_new,d_waveH,waveH_input
   REALTYPE,dimension(:,:),pointer     :: waveL_new,d_waveL,waveL_input
   REALTYPE,dimension(:,:),pointer     :: coswavedir_new,d_coswavedir
   REALTYPE,dimension(:,:),pointer     :: sinwavedir_new,d_sinwavedir

   character(len=10)         :: name_waveH="hs"
   character(len=10)         :: name_waveL="L"
   character(len=10)         :: name_waveDir="theta0"
!
! !REVISION HISTORY:
!  Original author(s): Ulf Gräwe
!                      Saeed Moghimi
!                      Knut Klingbeil
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_waves_input_ncdf -
!
! !INTERFACE:
   subroutine init_waves_input_ncdf(fn,nstart)
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
   write(debug,*) 'init_waves_input_ncdf() # ',Ncall
#endif

   call open_waves_file(waves_file)

   if (ilen.eq.iextr .and. jlen.eq.jextr) then
      LEVEL3 'Assuming On-Grid waves forcing'
      on_grid = .true.
      il = ilg ; jl = jlg ; ih = ihg ; jh = jhg
   else if (ilen.eq.1 .and. jlen.eq.1) then
      LEVEL3 'Assuming Point Source waves forcing'
      point_source = .true.
      on_grid = .true.
      il = 1 ; jl = 1 ; ih = 1 ; jh = 1
   else
      on_grid = .false.
      il = 1 ; jl = 1 ; ih = ilen ; jh = jlen

      allocate(ti(E2DFIELD),stat=rc)
      if (rc /= 0) &
          stop 'init_waves_input_ncdf: Error allocating memory (ti)'
      ti = -999.

      allocate(ui(E2DFIELD),stat=rc)
      if (rc /= 0) stop &
              'init_waves_input_ncdf: Error allocating memory (ui)'
      ui = -999.

      allocate(gridmap(E2DFIELD,1:2),stat=rc)
      if (rc /= 0) stop &
              'init_waves_input_ncdf: Error allocating memory (gridmap)'
      gridmap(:,:,:) = -999

      allocate(beta(E2DFIELD),stat=rc)
      if (rc /= 0) &
          stop 'init_waves_input_ncdf: Error allocating memory (beta)'
      beta = _ZERO_

      call init_grid_interpol(imin,imax,jmin,jmax,az,  &
                lonc,latc,waves_lon,waves_lat,southpole,gridmap,beta,ti,ui, &
                break_on_missing=.false.)

   end if

   start(1) = il; start(2) = jl;
   edges(1) = ih-il+1; edges(2) = jh-jl+1;
   edges(3) = 1

   allocate(waveH_new(E2DFIELD),stat=rc)
   if (rc /= 0) call getm_error('init_waves_input_ncdf()',             &
                                'Error allocating memory (waveH_new)')

   if ( stationary ) then

      waveH_input => waveH_new
      waveL_input => waveL
      call read_data(0)
      waveH = waveH_new

   else

      allocate(d_waveH(E2DFIELD),stat=rc)
      if (rc /= 0) call getm_error('init_waves_input_ncdf()',          &
                             'Error allocating memory (d_waveH)')
      allocate(waveL_new(E2DFIELD),stat=rc)
      if (rc /= 0) call getm_error('init_waves_input_ncdf()',          &
                             'Error allocating memory (waveL_new)')
      allocate(d_waveL(E2DFIELD),stat=rc)
      if (rc /= 0) call getm_error('init_waves_input_ncdf()',          &
                             'Error allocating memory (d_waveL)')
      allocate(coswavedir_new(E2DFIELD),stat=rc)
      if (rc /= 0) call getm_error('init_waves_input_ncdf()',          &
                             'Error allocating memory (coswavedir_new)')
      allocate(d_coswavedir(E2DFIELD),stat=rc)
      if (rc /= 0) call getm_error('init_waves_input_ncdf()',          &
                             'Error allocating memory (d_coswavedir)')
      allocate(sinwavedir_new(E2DFIELD),stat=rc)
      if (rc /= 0) call getm_error('init_waves_input_ncdf()',          &
                             'Error allocating memory (sinwavedir_new)')
      allocate(d_sinwavedir(E2DFIELD),stat=rc)
      if (rc /= 0) call getm_error('init_waves_input_ncdf()',          &
                             'Error allocating memory (d_sinwavedir)')

      waveH_input => d_waveH
      waveL_input => d_waveL
      call get_waves_data_ncdf(nstart-1)

   end if

#ifdef DEBUG
   write(debug,*) 'Leaving init_waves_input_ncdf()'
   write(debug,*)
#endif
   return
   end subroutine init_waves_input_ncdf
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_waves_data_ncdf - .
!
! !INTERFACE:
   subroutine get_waves_data_ncdf(loop)
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
   REALTYPE,dimension(:,:),pointer :: waveH_old,waveL_old
   REALTYPE,dimension(:,:),pointer :: coswavedir_old,sinwavedir_old
   logical, save                   :: first=.true.
!EOP
!-------------------------------------------------------------------------
#ifdef DEBUG
   integer, save   :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'get_waves_data_ncdf() # ',Ncall
#endif

   if (stationary) then
      waveH = waveH_new
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
         t2 = wave_times(indx) - offset
         if ( t2 .ge. t ) then
            EXIT
         end if
      end do

!     end of simulation?
      if (indx .gt. nlen) then
!        Note (KK): here we are not in case of first
!                   (because of in_interval check in open_waves_file)
         LEVEL2 'Need new waves file'
         call open_waves_file(waves_file)
         do indx=1,nlen
            t2 = wave_times(indx) - offset
!           Note (KK): For waves there is no check for long enough
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
         t2 = wave_times(indx) - offset
      end if

      call read_data(indx)
      save_n = indx+1

      waveH_old=>waveH_new;waveH_new=>d_waveH;d_waveH=>waveH_old;waveH_input=>d_waveH
      waveL_old=>waveL_new;waveL_new=>d_waveL;d_waveL=>waveL_old;waveL_input=>d_waveL
      coswavedir_old=>coswavedir_new;coswavedir_new=>coswavedir;coswavedir=>d_coswavedir;d_coswavedir=>coswavedir_old
      sinwavedir_old=>sinwavedir_new;sinwavedir_new=>sinwavedir;sinwavedir=>d_sinwavedir;d_sinwavedir=>sinwavedir_old

      if ( .not. first ) then
         d_waveH = waveH_new - waveH_old
         d_waveL = waveL_new - waveL_old
         d_coswavedir = coswavedir_new - coswavedir_old
         d_sinwavedir = sinwavedir_new - sinwavedir_old
         deltm1 = _ONE_ / (t2 - t1)
      end if

   end if


   t_minus_t2 = t - t2
   waveH      = waveH_new      + d_waveH     *deltm1*t_minus_t2
   waveL      = waveL_new      + d_waveL     *deltm1*t_minus_t2
   coswavedir = coswavedir_new + d_coswavedir*deltm1*t_minus_t2
   sinwavedir = sinwavedir_new + d_sinwavedir*deltm1*t_minus_t2

   first = .false.

#ifdef DEBUG
   write(debug,*) 'Leaving get_waves_data_ncdf()'
   write(debug,*)
#endif
   return
   end subroutine get_waves_data_ncdf
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: open_waves_file - .
!
! !INTERFACE:
   subroutine open_waves_file(waves_file)
!
! !DESCRIPTION:
!  Instead of specifying the name of the wavesrological file directly - a list
!  of names can be specified in \emph{waves\_file}. The rationale for this
!  approach is that output from operational wavesrological models are of
!  typically 2-5 days length. Collecting a number of these files allows for
!  longer model integrations without have to reformat the data.
!  It is assumed that the different files contains the same variables
!  and that they are of the same shape.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: waves_file
!
! !LOCAL VARIABLES:
   integer, parameter        :: iunit=60
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
   character(len=16)  :: dim_name(3),name_waves_mask
   integer            :: id
   logical            :: have_southpole
!
!EOP
!-------------------------------------------------------------------------
#ifdef DEBUG
   integer, save   :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'open_waves_file() # ',Ncall
#endif

   if (first) open(iunit,file=waves_file,status='old',action='read',err=80)

   found = .false.


   do

      if (.not. first_open) then
         err = nf90_close(ncid)
         if (err .NE. NF90_NOERR) go to 10
      end if

      read(iunit,*,err=85,end=90) fn
      LEVEL3 'Trying waves from:'
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
            LEVEL4 'stationary wave fields'
            stationary=.true.
            if (.not. first_open) call getm_error('open_waves_file()', &
                        'stationary fields only possible in first file')
         case(3)
            LEVEL4 'non-stationary wave fields'
            stationary=.false.
         case default
           call getm_error('open_waves_file()', &
                           'invalid number of dimensions')
      end select

      LEVEL4 ' ... checking variable ',name_waveH
      err = nf90_inq_varid(ncid,name_waveH,waveH_id)
      if (err .NE. NF90_NOERR) go to 10
      err = nf90_inquire_variable(ncid,waveH_id,ndims=nvardims)
      if (err .NE. NF90_NOERR) go to 10
      if (nvardims .NE. ndims) call getm_error('open_waves_file()',    &
                                 'Wrong number of dims in '//name_waveH)
      err = nf90_inquire_variable(ncid,waveH_id,dimids=vardim_ids)
      if (err .NE. NF90_NOERR) go to 10

      lon_dim = vardim_ids(1)
      lat_dim = vardim_ids(2)

      if (stationary) exit

      time_dim = vardim_ids(3)

      err = nf90_inq_varid(ncid,dim_name(time_dim),time_id)
      if (err .ne. NF90_NOERR) go to 10

      if (dim_len(time_dim) > nlen) then
         if (.not. first) then
            deallocate(wave_times,stat=err)
            if (err /= 0) call getm_error('open_waves_file()',         &
                              'Error de-allocating memory (wave_times)')
         end if
         allocate(wave_times(dim_len(time_dim)),stat=err)
         if (err /= 0) call getm_error('open_waves_file()',            &
                                 'Error allocating memory (wave_times)')
      end if
      nlen = dim_len(time_dim)
      err = nf90_get_var(ncid,time_id,wave_times(1:nlen))
      if (err .ne. NF90_NOERR) go to 10
      err =  nf90_get_att(ncid,time_id,'units',time_units)
      if (err .NE. NF90_NOERR) go to 10
      call string_to_julsecs(time_units,junit,sunit)

      offset = time_diff(jul0,secs0,junit,sunit)

      call add_secs(junit,sunit,nint(wave_times(1   )),j1,s1)
      call add_secs(junit,sunit,nint(wave_times(nlen)),j2,s2)

      if (first) then
         if (in_interval(j1,s1,julianday,secondsofday,j2,s2)) then
            found = .true.
         end if
      else
         if (time_diff(j2,s2,julianday,secondsofday) > _ZERO_) then
            found = .true.
         else
            LEVEL0 'WARNING: skipping waves file ',trim(fn)
         end if
      end if

      if (found) exit
      first_open = .false.

   end do

   if ( .not. found ) call getm_error('open_waves_file()',             &
               'Could not find valid waveforcing in '//trim(waves_file))

   LEVEL4 ' ... checking variable ',name_waveL
   err = nf90_inq_varid(ncid,name_waveL,waveL_id)
   if (err .NE. NF90_NOERR) go to 10
   err = nf90_inquire_variable(ncid,waveL_id,ndims=nvardims)
   if (err .NE. NF90_NOERR) go to 10
   if (nvardims .NE. ndims) call getm_error('open_waves_file()',       &
                                 'Wrong number of dims in '//name_waveL)

   LEVEL4 ' ... checking variable ',name_waveDir
   err = nf90_inq_varid(ncid,name_waveDir,waveDir_id)
   if (err .NE. NF90_NOERR) go to 10
   err = nf90_inquire_variable(ncid,waveDir_id,ndims=nvardims)
   if (err .NE. NF90_NOERR) go to 10
   if (nvardims .NE. ndims) call getm_error('open_waves_file()',       &
                               'Wrong number of dims in '//name_waveDir)

   LEVEL3 'Using waves from:'
   LEVEL4 trim(fn)
   if ( .not. stationary ) then
      LEVEL3 'waves offset time ',offset
   end if
   ilen = dim_len(lon_dim)
   jlen = dim_len(lat_dim)

   if (first) then

      err = nf90_inq_varid(ncid,dim_name(lon_dim),lon_id)
      if (err .ne. NF90_NOERR) go to 10
      allocate(waves_lon(ilen),stat=err)
      if (err /= 0) call getm_error('open_waves_file()',                 &
                                    'Error allocating memory (waves_lon)')
      err = nf90_get_var(ncid,lon_id,waves_lon(1:ilen))
      if (err .ne. NF90_NOERR) go to 10

      err = nf90_inq_varid(ncid,dim_name(lat_dim),lat_id)
      if (err .ne. NF90_NOERR) go to 10
      allocate(waves_lat(ilen),stat=err)
      if (err /= 0) call getm_error('open_waves_file()',                 &
                                    'Error allocating memory (waves_lat)')
      err = nf90_get_var(ncid,lat_id,waves_lat(1:jlen))
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
               rotated_waves_grid = .true.
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
                  rotated_waves_grid = .true.
               end if
            end if
            if (rotated_waves_grid) then
               LEVEL4 'south pole:'
!              changed indices - kb 2014-12-15
               LEVEL4 '      lon ',southpole(2)
               LEVEL4 '      lat ',southpole(1)
            end if

      allocate(waves_mask(1:ilen,1:jlen),stat=err)
      if (err /= 0) &
         stop 'open_waves_file: Error allocating memory (waves_mask)'
      err =  nf90_get_att(ncid,waveH_id,'mask',name_waves_mask)
      if (err .ne. NF90_NOERR) name_waves_mask='waves_mask'
      err = nf90_inq_varid(ncid,trim(name_waves_mask),id)
      if (err .eq. NF90_NOERR) then
         LEVEL4 'taking variable' // trim(name_waves_mask) // ' as waves_mask'
         err = nf90_get_var(ncid,id,waves_mask)
         if (err .ne. NF90_NOERR) go to 10
      else
         waves_mask = 1
      end if

   end if


   first = .false.

#ifdef DEBUG
   write(debug,*) 'Leaving open_waves_file()'
   write(debug,*)
#endif
   return

10 FATAL 'open_waves_file: ',nf90_strerror(err)
   stop 'open_waves_file()'
80 FATAL 'I could not open: ',trim(waves_file)
   stop 'open_waves_file()'
85 FATAL 'Error reading: ',trim(waves_file)
   stop 'open_waves_file()'
90 FATAL 'Reached eof in: ',trim(waves_file)
   stop 'open_waves_file()'

   end subroutine open_waves_file
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
   REALTYPE,dimension(edges(1),edges(2)) :: wrk,wd
   REALTYPE,dimension(E2DFIELD) :: cwd,swd
   integer                      :: err
   REALTYPE, parameter :: pi=3.1415926535897932384626433832795029d0
   REALTYPE, parameter :: deg2rad=pi/180.,rad2deg=180./pi
!EOP
!-----------------------------------------------------------------------

   call write_time_string()
   LEVEL3 timestr,': reading wave data ...',indx
   start(3) = indx

   err = nf90_get_var(ncid,waveH_id,wrk,start,edges)
   if (err .ne. NF90_NOERR) go to 10
   if (on_grid) then
      if (point_source) then
         waveH_input = wrk(1,1)
      else
         waveH_input(ill:ihl,jll:jhl) = wrk
      end if
   else
      call do_grid_interpol(az,wrk,gridmap,ti,ui,waveH_input,          &
                            imask=waves_mask,fillvalue=_ZERO_)
   end if

   err = nf90_get_var(ncid,waveL_id,wrk,start,edges)
   if (err .ne. NF90_NOERR) go to 10
   if (on_grid) then
      if (point_source) then
         waveL_input = wrk(1,1)
      else
         waveL_input(ill:ihl,jll:jhl) = wrk
      end if
   else
      call do_grid_interpol(az,wrk,gridmap,ti,ui,waveL_input,          &
                            imask=waves_mask,fillvalue=_ZERO_)
   end if

   err = nf90_get_var(ncid,waveDir_id,wrk,start,edges)
   if (err .ne. NF90_NOERR) go to 10
   if (on_grid) then
      if (point_source) then
         cwd = cos(wrk(1,1))
         swd = sin(wrk(1,1))
      else
         cwd(ill:ihl,jll:jhl) = cos(wrk)
         swd(ill:ihl,jll:jhl) = sin(wrk)
      end if
   else
      wd = wrk*deg2rad
      wrk = cos(wd)
      call do_grid_interpol(az,wd,gridmap,ti,ui,cwd,                   &
                            imask=waves_mask,fillvalue=_ZERO_)
      wrk = sin(wd)
      call do_grid_interpol(az,wd,gridmap,ti,ui,swd,                   &
                            imask=waves_mask,fillvalue=_ZERO_)
   end if
   coswavedir = cwd*cosconv - swd*sinconv
   sinwavedir = swd*cosconv + cwd*sinconv

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

   end module ncdf_waves

!-----------------------------------------------------------------------
! Copyright (C) 2014 - Knut Klingbeil                                  !
!-----------------------------------------------------------------------
