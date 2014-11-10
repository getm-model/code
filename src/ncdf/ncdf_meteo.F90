#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdf_meteo -
!
! !INTERFACE:
   module ncdf_meteo
!
! !DESCRIPTION:
!
! !USES:
   use netcdf
   use time, only: string_to_julsecs,time_diff,add_secs,in_interval
   use time, only: jul0,secs0,julianday,secondsofday,timestep,simtime
   use time, only: write_time_string,timestr
   use domain, only: imin,imax,jmin,jmax,az,lonc,latc,convc
   use domain, only: iextr_domain=>iextr,jextr_domain=>jextr
   use domain, only: ill,ihl,jll,jhl,ilg,ihg,jlg,jhg
   use grid_interpol, only: init_grid_interpol,do_grid_interpol
   use grid_interpol, only: to_rotated_lat_lon
   use meteo, only: meteo_file,on_grid,calc_met,met_method,hum_method
   use meteo, only: airp,u10,v10,t2,hum,tcc
   use meteo, only: fwf_method,evap,precip
   use meteo, only: tausx,tausy,swr,shf
   use meteo, only: new_meteo,t_1,t_2
   use meteo, only: wind_factor,evap_factor,precip_factor
   use meteo, only: nudge_sst,sst
   use exceptions
   IMPLICIT NONE
!
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_meteo_input_ncdf,get_meteo_data_ncdf
!
! !PRIVATE DATA MEMBERS:
   REALTYPE        :: offset
   integer         :: ncid,ndims,dims(3)
   integer         :: start(3),edges(3)
   integer         :: u10_id,v10_id,airp_id,t2_id
   integer         :: hum_id,convp_id,largep_id,tcc_id
   integer         :: evap_id=-1,precip_id=-1
   integer         :: tausx_id,tausy_id,swr_id,shf_id
   integer         :: sst_id=-1
   integer         :: iextr,jextr,textr,tmax=-1
   integer         :: grid_scan=1
   logical         :: point_source=.false.
   logical         :: rotated_meteo_grid=.false.

   REALTYPE, allocatable     :: met_lon(:),met_lat(:)
   REALTYPE, allocatable      :: met_times(:)
   REALTYPE, allocatable      :: wrk(:,:)
   REALTYPE, allocatable     :: wrk_dp(:,:)

!  For gridinterpolation
   REALTYPE, allocatable     :: beta(:,:)
   REALTYPE, allocatable     :: ti(:,:),ui(:,:)
   integer, allocatable      :: gridmap(:,:,:)
!
   REALTYPE, parameter       :: pi=3.1415926535897932384626433832795029
   REALTYPE, parameter       :: deg2rad=pi/180.,rad2deg=180./pi
   REALTYPE                  :: southpole(3) = (/0.0,-90.0,0.0/)
   character(len=10)         :: name_lon="lon"
   character(len=10)         :: name_lat="lat"
   character(len=10)         :: name_time="time"
   character(len=10)         :: name_u10="u10"
   character(len=10)         :: name_v10="v10"
   character(len=10)         :: name_airp="slp"
   character(len=10)         :: name_t2="t2"
   character(len=10)         :: name_hum1="sh"
   character(len=10)         :: name_hum2="rh"
   character(len=10)         :: name_hum3="dev2"
   character(len=10)         :: name_hum4="twet"
   character(len=10)         :: name_tcc="tcc"
   character(len=10)         :: name_evap="evap"
   character(len=10)         :: name_precip="precip"
   integer, parameter        :: SPECIFIC_HUM=1
   integer, parameter        :: RELATIVE_HUM=2
   integer, parameter        :: DEW_POINT=3
   integer, parameter        :: WET_BULB=4

   character(len=10)         :: name_tausx="tausx"
   character(len=10)         :: name_tausy="tausy"
   character(len=10)         :: name_swr="swr"
   character(len=10)         :: name_shf="shf"
   character(len=10)         :: name_sst="sst"
   character(len=128)        :: model_time
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
! !IROUTINE: init_meteo_input_ncdf -
!
! !INTERFACE:
   subroutine init_meteo_input_ncdf(fn,nstart)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Prepares reading meteorological forcing from a NetCDF formatted file.
!  Based on names of various variables the corresponding variable ids
!  are obtained from the NetCDF file.
!  The dimensions of the meteological grid is read (x,y,t).
!  If the southpole is not (0,-90,0) a rotated grid is assumed and
!  coefficients for interpolation between the meteorological grid and
!  the model grid are calculated.
!  The arry \emph{met\_times} are filled with the times where forcing is
!  available.
!  Finally, meteorological fields are initialised by a call to
!  \emph{get\_meteo\_data\_ncdf}.
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn
   integer, intent(in)                 :: nstart
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   integer         :: i,j,il,ih,jl,jh
   integer         :: err
   logical         :: ok=.true.
   REALTYPE        :: olon,olat,rlon,rlat,x
!EOP
!-------------------------------------------------------------------------
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_meteo_input_ncdf() # ',Ncall
#endif

   call open_meteo_file(meteo_file)

   if (met_lat(1) .gt. met_lat(2)) then
      LEVEL3 'Reverting lat-axis and setting grid_scan to 0'
      grid_scan = 0
      x = met_lat(1)
      do j=1,jextr/2
         met_lat(j) = met_lat(jextr-j+1)
         met_lat(jextr-j+1) = x
         x = met_lat(j+1)
      end do
   end if

   if (on_grid) then
      if (iextr_domain.ne.iextr .or. jextr_domain.ne.jextr) then
         call getm_error("init_meteo_input_ncdf()", &
                         "dimensions do not match")
      end if
      il = ilg ; jl = jlg ; ih = ihg ; jh = jhg
   else
      il = 1 ; jl = 1 ; ih = iextr ; jh = jextr
   end if

   start(1) = il; start(2) = jl;
   edges(1) = ih-il+1; edges(2) = jh-jl+1;
   edges(3) = 1

   allocate(wrk(edges(1),edges(2)),stat=err)
   if (err /= 0) stop 'ncdf_meteo: Error allocating memory (wrk)'
   wrk = 0.

   allocate(beta(E2DFIELD),stat=err)
   if (err /= 0) &
       stop 'init_meteo_input_ncdf: Error allocating memory (beta)'
   beta = _ZERO_

   if(iextr .eq. 1 .and. jextr .eq. 1) then
      point_source = .true.
      LEVEL3 'Assuming Point Source meteo forcing'
      if ( .not. on_grid) then
         LEVEL3 'Setting on_grid to true'
         on_grid=.true.
      end if
      if (rotated_meteo_grid) then
         rlon=met_lon(1)
         rlat=met_lat(1)
         call to_rotated_lat_lon(southpole,olon,olat,rlon,rlat,x)
         beta = x
      end if
   end if


   if ( .not. on_grid ) then

      if (.not. calc_met) then
         stop 'init_meteo_input_ncdf: calc_met=false requires on_grid=true'
      end if

      allocate(wrk_dp(edges(1),edges(2)),stat=err)
      if (err /= 0) stop 'ncdf_meteo: Error allocating memory (wrk_dp)'
      wrk_dp = _ZERO_

      allocate(ti(E2DFIELD),stat=err)
      if (err /= 0) &
          stop 'init_meteo_input_ncdf: Error allocating memory (ti)'
      ti = -999.

      allocate(ui(E2DFIELD),stat=err)
      if (err /= 0) stop &
              'init_meteo_input_ncdf: Error allocating memory (ui)'
      ui = -999.

      allocate(gridmap(E2DFIELD,1:2),stat=err)
      if (err /= 0) stop &
              'init_meteo_input_ncdf: Error allocating memory (gridmap)'
      gridmap(:,:,:) = -999
#if 0
#ifdef MED_15X15MINS_TEST
      do i=1,iextr
         met_lon(i) = -10.125 + (i-1)*1.125
      end do
      do j=1,jextr
         met_lat(j) =  28.125 + (j-1)*1.125
      end do
#endif

#ifdef NS_06NM_TEST
      do i=1,iextr
         met_lon(i) = -21.0 + (i-1)*1.
      end do
      do j=1,jextr
         met_lat(j) =  48.0 + (j-1)*1.
      end do
      grid_scan=0
#endif
#endif

      call init_grid_interpol(imin,imax,jmin,jmax,az,  &
                lonc,latc,met_lon,met_lat,southpole,gridmap,beta,ti,ui)

      LEVEL2 "Checking interpolation coefficients"
      do j=jmin,jmax
         do i=imin,imax
            if ( az(i,j) .gt. 0 .and. &
                (ui(i,j) .lt. _ZERO_ .or. ti(i,j) .lt. _ZERO_ )) then
               ok=.false.
               LEVEL3 "error at (i,j) ",i,j
            end if
         end do
      end do
      if ( ok ) then
         LEVEL2 "done"
      else
         call getm_error("init_meteo_input_ncdf()", &
                          "Some interpolation coefficients are not valid")
      end if
   end if

   call get_meteo_data_ncdf(nstart-1)

#ifdef DEBUG
   write(debug,*) 'Leaving init_meteo_input_ncdf()'
   write(debug,*)
#endif
   return
   end subroutine init_meteo_input_ncdf
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_meteo_data_ncdf - .
!
! !INTERFACE:
   subroutine get_meteo_data_ncdf(loop)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Do book keeping about when new fields are to be read. Set variables
!  used by \emph{do\_meteo} and finally calls \emph{read\_data} if
!  necessary.
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: loop
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   integer         :: indx
   REALTYPE        :: t
   logical, save   :: first=.true.
   integer, save   :: save_n=1
!EOP
!-------------------------------------------------------------------------
#ifdef DEBUG
   integer, save   :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'get_meteo_data_ncdf() # ',Ncall
#endif

   if (met_method .eq. 2) then

!     find the right index

      t = loop*timestep
      new_meteo = .false.

      if ( t .gt. t_2 ) then

         t_1 = t_2
         new_meteo = .true.

!        Note (KK): Even if the first time stage is the last entry of an
!                   input file, we must not jump to the next file because
!                   then indx-1 does not work!
!                   Therefore .ge. and save_n must be initialised to 1!
         do indx=save_n,tmax
            t_2 = met_times(indx) - offset
            if ( t_2 .ge. t ) then
               EXIT
            end if
         end do

!        end of simulation?
         if (indx .gt. tmax) then
!           Note (KK): here we are not in case of first
!                      (because of in_interval check in open_meteo_file)
            LEVEL2 'Need new meteo file'
            call open_meteo_file(meteo_file)
            do indx=1,tmax
               t_2 = met_times(indx) - offset
!              Note (KK): For meteo there is no check for long enough
!                         data sets. Therefore .ge.!
               if ( t_2 .ge. t ) then
                  EXIT
               end if
            end do
         end if

         if (first) then
            first = .false.
!            if (indx .gt. 1) then
            if ( t_2 .gt. t ) then
               indx = indx-1
            end if
            t_2 = met_times(indx) - offset
         end if

         start(3) = indx
         call write_time_string()
         LEVEL3 timestr,': reading meteo data ...',indx
         call read_data()
         save_n = indx+1

      end if
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving get_meteo_data_ncdf()'
   write(debug,*)
#endif
   return
   end subroutine get_meteo_data_ncdf
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: open_meteo_file - .
!
! !INTERFACE:
   subroutine open_meteo_file(meteo_file)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Instead of specifying the name of the meteorological file directly - a list
!  of names can be specified in \emph{meteo\_file}. The rationale for this
!  approach is that output from operational meteorological models are of
!  typically 2-5 days length. Collecting a number of these files allows for
!  longer model integrations without have to reformat the data.
!  It is assumed that the different files contains the same variables
!  and that they are of the same shape.
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: meteo_file
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   integer, parameter        :: iunit=55
   character(len=256)        :: fn,time_units
   integer         :: junit,sunit,j1,s1,j2,s2
   integer         :: n,err,idum
   logical         :: first=.true.
   logical         :: found=.false.,first_open=.true.
   integer, save   :: lon_id=-1,lat_id=-1,time_id=-1,id=-1
   integer, save   :: time_var_id=-1
   character(len=256) :: dimname
!
!EOP
!-------------------------------------------------------------------------
#ifdef DEBUG
   integer, save   :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'open_meteo_file() # ',Ncall
#endif

   if (first) then
      first = .false.
      open(iunit,file=meteo_file,status='old',action='read',err=80)
      do
         if (found) EXIT
         read(iunit,*,err=85,end=90) fn
         LEVEL3 'Trying meteo from:'
         LEVEL4 trim(fn)
         err = nf90_open(fn,NF90_NOWRITE,ncid)
         if (err .ne. NF90_NOERR) go to 10

         if (first_open) then
            first_open = .false.
            err = nf90_inquire(ncid,nDimensions=ndims)
            if (err .NE. NF90_NOERR) go to 10

            LEVEL4 'dimensions'
            do n=1,ndims
               err = nf90_inquire_dimension(ncid,n,name=dimname)
               if (err .ne. NF90_NOERR) go to 10
               if( dimname .eq. name_lon ) then
                  lon_id = n
                  err = nf90_inquire_dimension(ncid,lon_id,len=iextr)
                  if (err .ne. NF90_NOERR) go to 10
                  LEVEL4 'lon_id  --> ',lon_id,', len = ',iextr
               end if
               if( dimname .eq. name_lat ) then
                  lat_id = n
                  err = nf90_inquire_dimension(ncid,lat_id,len=jextr)
                  if (err .ne. NF90_NOERR) go to 10
                  LEVEL4 'lat_id  --> ',lat_id,', len = ',jextr
               end if
               if( dimname .eq. name_time ) then
                  time_id = n
                  err = nf90_inquire_dimension(ncid,time_id,len=textr)
                  if (err .ne. NF90_NOERR) go to 10
                  LEVEL4 'time_id --> ',time_id,', len = ',textr
!                  if (tmax .lt. 0) tmax=textr
                  tmax=textr
               end if
            end do
            if(lon_id .eq. -1) then
               FATAL 'could not find longitude coordinate in meteo file'
               stop 'open_meteo_file()'
            end if
            if(lat_id .eq. -1) then
               FATAL 'could not find latitude coordinate in meteo file'
               stop 'open_meteo_file()'
            end if
            if(time_id .eq. -1) then
               FATAL 'could not find time coordinate in meteo file'
               stop 'open_meteo_file()'
            end if

            allocate(met_lon(iextr),stat=err)
            if (err /= 0) stop &
                  'open_meteo_file(): Error allocating memory (met_lon)'
            err = nf90_inq_varid(ncid,name_lon,id)
            if (err .NE. NF90_NOERR) go to 10
            err = nf90_get_var(ncid,id,met_lon)
            if (err .ne. NF90_NOERR) go to 10

            allocate(met_lat(jextr),stat=err)
            if (err /= 0) stop &
                  'open_meteo_file(): Error allocating memory (met_lat)'
            err = nf90_inq_varid(ncid,name_lat,id)
            if (err .NE. NF90_NOERR) go to 10
            err = nf90_get_var(ncid,id,met_lat)
            if (err .ne. NF90_NOERR) go to 10

            allocate(met_times(textr),stat=err)
            if (err /= 0) stop &
                  'open_meteo_file(): Error allocating memory (met_times)'

            err = nf90_inq_varid(ncid,'southpole',id)
            if (err .ne. NF90_NOERR) then
               LEVEL4 'Setting southpole to (0,-90,0)'
            else
               err = nf90_get_var(ncid,id,southpole)
               if (err .ne. NF90_NOERR) go to 10
               rotated_meteo_grid = .true.
            end if
            LEVEL4 'south pole:'
            LEVEL4 '      lon ',southpole(1)
            LEVEL4 '      lat ',southpole(2)

         end if

         err = nf90_inquire_dimension(ncid,time_id,len=idum)
         if (err .ne. NF90_NOERR) go to 10
         if(idum .gt. size(met_times)) then
            deallocate(met_times,stat=err)
            if (err /= 0) stop      &
               'open_meteo_file(): Error de-allocating memory (met_times)'
            allocate(met_times(idum),stat=err)
            if (err /= 0) stop &
               'open_meteo_file(): Error allocating memory (met_times)'
         end if
         textr = idum
         LEVEL3 'time_id --> ',time_id,', len = ',textr
!        if (tmax .lt. 0) tmax=textr
         tmax=textr

         err = nf90_inq_varid(ncid,name_time,time_var_id)
         if (err .NE. NF90_NOERR) go to 10
         err =  nf90_get_att(ncid,time_var_id,'units',time_units)
         if (err .NE. NF90_NOERR) go to 10
         call string_to_julsecs(time_units,junit,sunit)
         err = nf90_get_var(ncid,time_var_id,met_times(1:textr))
         if (err .ne. NF90_NOERR) go to 10

         call add_secs(junit,sunit,nint(met_times(1)),    j1,s1)
         call add_secs(junit,sunit,nint(met_times(textr)),j2,s2)

         if (in_interval(j1,s1,julianday,secondsofday,j2,s2)) then
            found = .true.
         else
            err = nf90_close(ncid)
            if (err .NE. NF90_NOERR) go to 10
         end if
      end do
   else
      err = nf90_close(ncid)
      if (err .NE. NF90_NOERR) go to 10
!     open next file
      read(iunit,*,err=85,end=90) fn
      err = nf90_open(fn,NF90_NOWRITE,ncid)
      if (err .ne. NF90_NOERR) go to 10

      err = nf90_inquire_dimension(ncid,time_id,len=idum)
      if (err .ne. NF90_NOERR) go to 10
      if(idum .gt. size(met_times)) then
         deallocate(met_times,stat=err)
         if (err /= 0) stop      &
            'open_meteo_file(): Error de-allocating memory (met_times)'
         allocate(met_times(idum),stat=err)
         if (err /= 0) stop &
            'open_meteo_file(): Error allocating memory (met_times)'
      end if
      textr = idum
      LEVEL3 'time_id --> ',time_id,', len = ',textr
!     if (tmax .lt. 0) tmax=textr
      tmax=textr

      err = nf90_inq_varid(ncid,name_time,time_var_id)
      if (err .NE. NF90_NOERR) go to 10
      err =  nf90_get_att(ncid,time_var_id,'units',time_units)
      if (err .NE. NF90_NOERR) go to 10
      call string_to_julsecs(time_units,junit,sunit)
      err = nf90_get_var(ncid,time_var_id,met_times(1:textr))
      if (err .ne. NF90_NOERR) go to 10

      call add_secs(junit,sunit,nint(met_times(1)),    j1,s1)
      call add_secs(junit,sunit,nint(met_times(textr)),j2,s2)
   end if

   LEVEL4 ' ... checking variable ',name_airp
   err = nf90_inq_varid(ncid,name_airp,airp_id)
   if (err .NE. NF90_NOERR) go to 10

   if (fwf_method .eq. 2) then
      LEVEL4 ' ... checking variable ',name_evap
      err = nf90_inq_varid(ncid,name_evap,evap_id)
      if (err .NE. NF90_NOERR) go to 10
   end if
   if (fwf_method .eq. 2 .or. fwf_method .eq. 3) then
      LEVEL4 ' ... checking variable ',name_precip
      err = nf90_inq_varid(ncid,name_precip,precip_id)
      if (err .NE. NF90_NOERR) go to 10
   end if

   if (calc_met) then
      LEVEL4 ' ... checking variable ',name_u10
      err = nf90_inq_varid(ncid,name_u10,u10_id)
      if (err .NE. NF90_NOERR) go to 10

      LEVEL4 ' ... checking variable ',name_v10
      err = nf90_inq_varid(ncid,name_v10,v10_id)
      if (err .NE. NF90_NOERR) go to 10

      LEVEL4 ' ... checking variable ',name_t2
      err = nf90_inq_varid(ncid,name_t2,t2_id)
      if (err .NE. NF90_NOERR) go to 10

      hum_id = -1
      err = nf90_inq_varid(ncid,name_hum1,hum_id)
      if (err .NE. NF90_NOERR) then
         err = nf90_inq_varid(ncid,name_hum2,hum_id)
         if (err .NE. NF90_NOERR) then
            err = nf90_inq_varid(ncid,name_hum3,hum_id)
            if (err .NE. NF90_NOERR) then
               FATAL 'Not able to find valid humidity parameter'
               stop 'init_meteo_input_ncdf()'
            else
               LEVEL2 'Taking hum as dew point temperature'
               hum_method = DEW_POINT
            end if
         else
            LEVEL2 'Taking hum as relative humidity'
            hum_method = RELATIVE_HUM
         end if
      else
         LEVEL2 'Taking hum as atmospheric specific humidity'
         hum_method = SPECIFIC_HUM
      end if
!KBKSTDERR 'Taking hum as wet bulb temperature'

      LEVEL4 ' ... checking variable ',name_tcc
      err = nf90_inq_varid(ncid,name_tcc,tcc_id)
      if (err .NE. NF90_NOERR) go to 10

   else

      LEVEL4 ' ... checking variable ',name_tausx
      err = nf90_inq_varid(ncid,name_tausx,tausx_id)
      if (err .NE. NF90_NOERR) go to 10

      LEVEL4 ' ... checking variable ',name_tausy
      err = nf90_inq_varid(ncid,name_tausy,tausy_id)
      if (err .NE. NF90_NOERR) go to 10

      LEVEL4 ' ... checking variable ',name_swr
      err = nf90_inq_varid(ncid,name_swr,swr_id)
      if (err .NE. NF90_NOERR) go to 10

      LEVEL4 ' ... checking variable ',name_shf
      err = nf90_inq_varid(ncid,name_shf,shf_id)
      if (err .NE. NF90_NOERR) go to 10

   end if

   if (nudge_sst) then
      LEVEL4 ' ... checking variable ',name_sst
      err = nf90_inq_varid(ncid,name_sst,sst_id)
      if (err .NE. NF90_NOERR) go to 10
   end if

   if (found) then
      offset = time_diff(jul0,secs0,junit,sunit)
      LEVEL3 'Using meteo from:'
      LEVEL4 trim(fn)
      LEVEL3 'Meteorological offset time ',offset
   else
      FATAL 'Could not find any valid meteo-files'
      stop 'open_meteo_file'
   end if

   return
10 FATAL 'open_meteo_file: ',nf90_strerror(err)
   stop 'open_meteo_file()'
80 FATAL 'I could not open: ',trim(meteo_file)
   stop 'open_meteo_file()'
85 FATAL 'Error reading: ',trim(meteo_file)
   stop 'open_meteo_file()'
90 FATAL 'Reached eof in: ',trim(meteo_file)
   stop 'open_meteo_file()'

#ifdef DEBUG
   write(debug,*) 'Leaving open_meteo_file()'
   write(debug,*)
#endif
   return
   end subroutine open_meteo_file
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_data -
!
! !INTERFACE:
   subroutine read_data()
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Reads the relevant variables from the NetCDF file. Interpolates to the
!  model grid if necessary. After a call to this routine updated versions
!  of either variables used for calculating stresses and fluxes or directly
!  the stresses/fluxes directly are available to \emph{do\_meteo}.
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   integer         :: i,j,err
   REALTYPE        :: angle,uu,vv,sinconv,cosconv
!EOP
!-----------------------------------------------------------------------

   err = nf90_get_var(ncid,airp_id,wrk,start,edges)
   if (err .ne. NF90_NOERR) go to 10
   if (on_grid) then
      if (point_source) then
         airp = wrk(1,1)
      else
         airp(ill:ihl,jll:jhl) = wrk
      end if
   else
      !KBKwrk_dp = _ZERO_
      call copy_var(grid_scan,wrk,wrk_dp)
      call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,airp)
   end if

   if (evap_id .ge. 0) then
      err = nf90_get_var(ncid,evap_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (on_grid) then
         if (point_source) then
            evap = wrk(1,1)
         else
            evap(ill:ihl,jll:jhl) = wrk
         end if
      else
         call copy_var(grid_scan,wrk,wrk_dp)
         call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,evap)
      end if
      if (evap_factor .ne. _ONE_) then
         evap = evap * evap_factor
      end if
   end if

   if (precip_id .gt. 0) then
      err = nf90_get_var(ncid,precip_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (on_grid) then
         if (point_source) then
            precip = wrk(1,1)
         else
            precip(ill:ihl,jll:jhl) = wrk
         end if
      else
         call copy_var(grid_scan,wrk,wrk_dp)
         call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,precip)
      end if
      if (precip_factor .ne. _ONE_) then
         precip = precip * precip_factor
      end if
   end if

   if (calc_met) then

      err = nf90_get_var(ncid,u10_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (on_grid) then
         if (point_source) then
            u10 = wrk(1,1)
         else
            u10(ill:ihl,jll:jhl) = wrk
         end if
      else
         !KBKwrk_dp = _ZERO_
         call copy_var(grid_scan,wrk,wrk_dp)
         call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,u10)
      end if

      err = nf90_get_var(ncid,v10_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (on_grid) then
         if (point_source) then
            v10 = wrk(1,1)
         else
            v10(ill:ihl,jll:jhl) = wrk
         end if
      else
         !KBKwrk_dp = _ZERO_
         call copy_var(grid_scan,wrk,wrk_dp)
         call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,v10)
      end if

!     Rotation of wind due to the combined effect of possible rotation of
!     meteorological grid and possible hydrodynamic grid convergence
!     (cartesian and curvi-linear grids where conv <> 0.)
      do j=jmin-1,jmax+1
         do i=imin-1,imax+1
!KBK            angle=-convc(i,j)*deg2rad
!KBK            angle=beta(i,j)
            angle=beta(i,j)-convc(i,j)*deg2rad
            if(angle .ne. _ZERO_) then
               sinconv=sin(angle)
               cosconv=cos(angle)
               uu=u10(i,j)
               vv=v10(i,j)
               u10(i,j)= uu*cosconv+vv*sinconv
               v10(i,j)=-uu*sinconv+vv*cosconv
            end if
         end do
      end do

      if (wind_factor .ne. _ONE_) then
         u10 = u10 * wind_factor
         v10 = v10 * wind_factor
      end if


      err = nf90_get_var(ncid,t2_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (on_grid) then
         if (point_source) then
            t2 = wrk(1,1)
         else
            t2(ill:ihl,jll:jhl) = wrk
         end if
      else
         !KBKwrk_dp = _ZERO_
         call copy_var(grid_scan,wrk,wrk_dp)
         call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,t2)
      end if

      err = nf90_get_var(ncid,hum_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (on_grid) then
         if (point_source) then
            hum = wrk(1,1)
         else
            hum(ill:ihl,jll:jhl) = wrk
         end if
      else
         !KBKwrk_dp = _ZERO_
         call copy_var(grid_scan,wrk,wrk_dp)
         call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,hum)
      end if

      err = nf90_get_var(ncid,tcc_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (on_grid) then
         if (point_source) then
            tcc = wrk(1,1)
         else
            tcc(ill:ihl,jll:jhl) = wrk
         end if
      else
         !KBKwrk_dp = _ZERO_
         call copy_var(grid_scan,wrk,wrk_dp)
         call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,tcc)
      end if

   else

      err = nf90_get_var(ncid,tausx_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (point_source) then
         tausx = wrk(1,1)
      else
         tausx(ill:ihl,jll:jhl) = wrk
      end if

      err = nf90_get_var(ncid,tausy_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (point_source) then
         tausy = wrk(1,1)
      else
         tausy(ill:ihl,jll:jhl) = wrk
      end if

      err = nf90_get_var(ncid,swr_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (point_source) then
         swr = wrk(1,1)
      else
         swr(ill:ihl,jll:jhl) = wrk
      end if

      err = nf90_get_var(ncid,shf_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (point_source) then
         shf = wrk(1,1)
      else
         shf(ill:ihl,jll:jhl) = wrk
      end if

   end if

   if (sst_id .gt. 0) then
      err = nf90_get_var(ncid,sst_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (on_grid) then
         if (point_source) then
            sst = wrk(1,1)
         else
            sst(ill:ihl,jll:jhl) = wrk
         end if
      else if (calc_met) then
         call copy_var(grid_scan,wrk,wrk_dp)
         call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,sst)
      end if
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
!BOP
!
! !IROUTINE: copy_var -
!
! !INTERFACE:
!   subroutine copy_var(grid_scan,var)
   subroutine copy_var(grid_scan,inf,outf)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Reads the relevant variables from the NetCDF file. Interpolates to the
!  model grid if necessary. After a call to this routine updated versions
!  of either variables used for calculating stresses and fluxes or directly
!  the stresses/fluxes directly are available to \emph{do\_meteo}.
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: grid_scan
   REALTYPE, intent(in)                :: inf(:,:)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: outf(:,:)
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   integer         :: i,j
!EOP
!-----------------------------------------------------------------------

   select case (grid_scan)
      case (0)
         do j=1,jextr
            do i=1,iextr
               outf(i,jextr-j+1) = inf(i,j)
            end do
         end do
      case (1) ! ?????
         do j=1,jextr
            do i=1,iextr
               outf(i,j) = inf(i,j)
            end do
         end do
      case default
         FATAL 'Do something here - copy_var'
   end select
   return
   end subroutine copy_var
!EOC

!-----------------------------------------------------------------------

   end module ncdf_meteo

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
