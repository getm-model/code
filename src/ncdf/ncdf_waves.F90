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
   use domain         ,only: az,convc
   use waves          ,only: waves_file,on_grid
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
   logical         :: stationary

   REALTYPE,dimension(:),allocatable   :: wave_times(:)
   REALTYPE,dimension(:,:),pointer     :: waveH_new,d_waveH
   REALTYPE,dimension(:,:),pointer     :: waveL_new,d_waveL
   REALTYPE,dimension(:,:),pointer     :: coswavedir_new,d_coswavedir
   REALTYPE,dimension(:,:),pointer     :: sinwavedir_new,d_sinwavedir
   REALTYPE,dimension(:,:),allocatable :: wrk

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

   if ( .not. on_grid ) then
      call getm_error("init_waves_input_ncdf()", &
                      "currently only on_grid=T supported")
   end if

   call open_waves_file(waves_file)

   if (on_grid) then
      if (ilen.ne.iextr .or. jlen.ne.jextr) then
         call getm_error("init_waves_input_ncdf()", &
                         "dimensions do not match")
      end if
      il = ilg ; jl = jlg ; ih = ihg ; jh = jhg
   else
      il = 1 ; jl = 1 ; ih = ilen ; jh = jlen
   end if

   start(1) = il; start(2) = jl;
   edges(1) = ih-il+1; edges(2) = jh-jl+1;
   edges(3) = 1

   allocate(wrk(edges(1),edges(2)),stat=rc)
   if (rc /= 0) call getm_error('init_waves_input_ncdf()',             &
                                'Error allocating memory (wrk)')

   allocate(waveH_new(E2DFIELD),stat=rc)
   if (rc /= 0) call getm_error('init_waves_input_ncdf()',             &
                                'Error allocating memory (waveH_new)')


   if ( stationary ) then

      call read_data(0)
      waveH_new = waveH

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

      waveH_old=>waveH_new;waveH_new=>waveH;waveH=>d_waveH;d_waveH=>waveH_old
      waveL_old=>waveL_new;waveL_new=>waveL;waveL=>d_waveL;d_waveL=>waveL_old
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
   integer            :: x_dim,y_dim,time_dim=-1,time_id=-1
   integer            :: dim_len(3),vardim_ids(3)
   character(len=16)  :: dim_name(3)
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
      x_dim = vardim_ids(1)
      y_dim = vardim_ids(2)
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
         allocate(wave_times(nlen),stat=err)
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
            LEVEL0 'WARNING: skipping meteo file ',trim(fn)
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
   ilen = dim_len(x_dim)
   jlen = dim_len(y_dim)

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
   REALTYPE,dimension(E2DFIELD) :: waveDir
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
      waveH(ill:ihl,jll:jhl) = wrk
   else
   end if

   err = nf90_get_var(ncid,waveL_id,wrk,start,edges)
   if (err .ne. NF90_NOERR) go to 10
   if (on_grid) then
      waveL(ill:ihl,jll:jhl) = wrk
   else
   end if

   err = nf90_get_var(ncid,waveDir_id,wrk,start,edges)
   if (err .ne. NF90_NOERR) go to 10
   if (on_grid) then
      waveDir(ill:ihl,jll:jhl) = wrk
   else
   end if
   waveDir = (waveDir - convc) * deg2rad
   coswavedir = cos(waveDir)
   sinwavedir = sin(waveDir)

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
! Copyright (C) 2014 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
