!$Id: ncdf_meteo.F90,v 1.15 2005-01-12 19:26:16 kbk Exp $
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
   use time, only: string_to_julsecs,time_diff,add_secs,in_interval
   use time, only: jul0,secs0,julianday,secondsofday,timestep
   use domain, only: imin,imax,jmin,jmax,az,lonc,latc,conv
   use grid_interpol, only: init_grid_interpol,do_grid_interpol
   use meteo, only: meteo_file,on_grid,calc_met,method,hum_method
   use meteo, only: airp,u10,v10,t2,hum,tcc
   use meteo, only: tausx,tausy,swr,shf
   use meteo, only: new_meteo,t_1,t_2
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
   integer         :: tausx_id,tausy_id,swr_id,shf_id
   integer         :: iextr,jextr,textr,tmax=-1
   integer         :: grid_scan=1
   logical         :: point_source=.false.

   REALTYPE, allocatable     :: met_lon(:),met_lat(:)
   REAL_4B, allocatable      :: met_times(:)
   REAL_4B, allocatable      :: wrk(:,:)
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
   integer, parameter        :: SPECIFIC_HUM=1
   integer, parameter        :: RELATIVE_HUM=2
   integer, parameter        :: DEW_POINT=3
   integer, parameter        :: WET_BULB=4

   character(len=10)         :: name_tausx="tausx"
   character(len=10)         :: name_tausy="tausy"
   character(len=10)         :: name_swr="swr"
   character(len=10)         :: name_shf="shf"
   character(len=128)        :: model_time
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ncdf_meteo.F90,v $
!  Revision 1.15  2005-01-12 19:26:16  kbk
!  fixed printing of south pole
!
!  Revision 1.14  2005/01/12 19:17:47  kbk
!  setting grid_scan depending on lat-axis - Stips
!
!  Revision 1.13  2004/08/09 10:43:59  kbk
!  correct length of met_times - Buchmann
!
!  Revision 1.12  2004/08/09 08:39:36  kbk
!  if SPHERICAL and rotated meteo grid fixed turning of wind - Carsten Hansen (FRV)
!
!  Revision 1.11  2004/04/06 16:32:29  kbk
!  TimeDiff --> time_diff
!
!  Revision 1.10  2004/01/15 11:45:01  kbk
!  meteo point source forcing - taus, swr and shf - implemented
!
!  Revision 1.9  2003/12/16 16:50:41  kbk
!  added support for Intel/IFORT compiler - expanded TABS, same types in subroutine calls
!
!  Revision 1.8  2003/11/03 14:34:54  kbk
!  use time_var_id in addition to time_id
!
!  Revision 1.7  2003/10/30 16:31:36  kbk
!  check validity of meteo interpolation coeffcients
!
!  Revision 1.6  2003/10/07 15:16:50  kbk
!  now works properly with varying length (time) files
!
!  Revision 1.5  2003/07/01 16:38:33  kbk
!  cleaned code - new methods
!
!  Revision 1.4  2003/06/17 14:53:29  kbk
!  default meteo variables names comply with Adolf Stips suggestion + southpole(3)
!
!  Revision 1.3  2003/04/07 15:34:15  kbk
!  updated to lonc,latc
!
!  Revision 1.1.1.1  2002/05/02 14:01:47  gotm
!  recovering after CVS crash
!
!  Revision 1.4  2001/10/17 14:27:39  bbh
!  Met-data can now be read from a series of .nc files
!
!  Revision 1.3  2001/07/26 13:57:14  bbh
!  Meteo working - needs some polishing
!
!  Revision 1.2  2001/06/04 13:15:12  bbh
!  Further steps towards full implementation of meteorological forcing
!
!  Revision 1.1  2001/05/25 19:26:22  bbh
!  ncdf_meteo.F90
!
! !TO DO:
!  Unified method of obtaining time info - needs some namelist variables.
!  Loop over met-files listed in meteo_file.
!  Make code independent of HIRLAM/ECMWF etc.
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
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   integer         :: i,j,n
   integer         :: err
   logical         :: ok=.true.
   REALTYPE        :: x
!EOP
!-------------------------------------------------------------------------
   include "netcdf.inc"
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_meteo_input_ncdf() # ',Ncall
#endif

   call open_meteo_file(meteo_file)

   if(iextr .eq. 1 .and. jextr .eq. 1) then
      point_source = .true.
      LEVEL3 'Assuming Point Source meteo forcing'
      if (on_grid .eq. .false. ) then
         LEVEL3 'Setting on_grid to true'
         on_grid=.true.
      end if
   end if

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

   allocate(wrk(iextr,jextr),stat=err)
   if (err /= 0) stop 'ncdf_meteo: Error allocating memory (wrk)'
   wrk = 0.

   allocate(wrk_dp(iextr,jextr),stat=err)
   if (err /= 0) stop 'ncdf_meteo: Error allocating memory (wrk_dp)'
   wrk_dp = _ZERO_

   if ( .not. on_grid ) then

      allocate(beta(E2DFIELD),stat=err)
      if (err /= 0) &
          stop 'init_meteo_input_ncdf: Error allocating memory (beta)'
      beta = _ZERO_

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

   if (calc_met) then

      err = nf_inq_varid(ncid,name_u10,u10_id)
      if (err .NE. NF_NOERR) go to 10

      err = nf_inq_varid(ncid,name_v10,v10_id)
      if (err .NE. NF_NOERR) go to 10

      err = nf_inq_varid(ncid,name_airp,airp_id)
      if (err .NE. NF_NOERR) go to 10

      err = nf_inq_varid(ncid,name_t2,t2_id)
      if (err .NE. NF_NOERR) go to 10

      hum_id = -1
      err = nf_inq_varid(ncid,name_hum1,hum_id)
      if (err .NE. NF_NOERR) then
         err = nf_inq_varid(ncid,name_hum2,hum_id)
         if (err .NE. NF_NOERR) then
            err = nf_inq_varid(ncid,name_hum3,hum_id)
            if (err .NE. NF_NOERR) then
               FATAL 'Not able to find valid humudity parameter'
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

      err = nf_inq_varid(ncid,name_tcc,tcc_id)
      if (err .NE. NF_NOERR) go to 10

   else

      err = nf_inq_varid(ncid,name_tausx,tausx_id)
      if (err .NE. NF_NOERR) go to 10

      err = nf_inq_varid(ncid,name_tausy,tausy_id)
      if (err .NE. NF_NOERR) go to 10

      err = nf_inq_varid(ncid,name_swr,swr_id)
      if (err .NE. NF_NOERR) go to 10

      err = nf_inq_varid(ncid,name_shf,shf_id)
      if (err .NE. NF_NOERR) go to 10

   end if

   if (method .eq. 2) then
      start(1) = 1; start(2) = 1;
      edges(1) = iextr; edges(2) = jextr;
      edges(3) = 1
      call get_meteo_data_ncdf(nstart)
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving init_meteo_input_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'init_meteo_input_ncdf: ',nf_strerror(err)
   stop
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
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   integer         :: i,indx
   REALTYPE        :: t
   logical, save   :: first=.true.
   integer, save   :: save_n=1
!
!EOP
!-------------------------------------------------------------------------
#ifdef DEBUG
   integer, save   :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'get_meteo_data_ncdf() # ',Ncall
#endif

   if (method .eq. 2) then

!     find the right index

      t = loop*timestep
      do indx=save_n,tmax
         if (met_times(indx) .gt. real(t + offset)) EXIT
      end do
      if (indx .gt. tmax) then
         LEVEL2 'Need new meteo file'
         call open_meteo_file(meteo_file)
         save_n = 1
         do indx=save_n,tmax
            if (met_times(indx) .gt. real(t + offset)) EXIT
         end do
! indx = 1
         save_n = 0
      end if
      start(3) = indx

      ! First time through we have to initialize t_1
      if (first) then
         first = .false.
         new_meteo = .true.
         if (indx .gt. 1) then
            indx = indx-1
         end if
         start(3) = indx
         t_1 = met_times(indx) - offset
         t_2 = t_1
         call read_data()
      else
         if (indx .gt. save_n) then
            new_meteo = .true.
            save_n = indx
            t_1 = t_2
            t_2 = met_times(indx) - offset
            call read_data()
         else
            new_meteo = .false.
         end if
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
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   integer, parameter        :: iunit=55
   character(len=256)        :: fn,time_units
   integer         :: j1,s1,j2,s2
   integer         :: n,err,idum
   logical         :: first=.true.
   logical         :: found=.false.,first_open=.true.
   integer, save   :: lon_id=-1,lat_id=-1,time_id=-1,id=-1
   integer, save   :: time_var_id=-1
   character(len=256) :: dimname
!
! !TO DO:
!  Need a variable to indicate homw much to read from each file.
!EOP
!-------------------------------------------------------------------------
   include "netcdf.inc"
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
         err = nf_open(fn,NCNOWRIT,ncid)
         if (err .ne. NF_NOERR) go to 10

         if (first_open) then
            first_open = .false.
            err = nf_inq_ndims(ncid,ndims)
            if (err .NE. NF_NOERR) go to 10

            LEVEL4 'dimensions'
            do n=1,ndims
               err = nf_inq_dimname(ncid,n,dimname)
               if (err .ne. NF_NOERR) go to 10
               if( dimname .eq. name_lon ) then
                  lon_id = n
                  err = nf_inq_dimlen(ncid,lon_id,iextr)
                  if (err .ne. NF_NOERR) go to 10
                  LEVEL4 'lon_id  --> ',lon_id,', len = ',iextr
               end if
               if( dimname .eq. name_lat ) then
                  lat_id = n
                  err = nf_inq_dimlen(ncid,lat_id,jextr)
                  if (err .ne. NF_NOERR) go to 10
                  LEVEL4 'lat_id  --> ',lat_id,', len = ',jextr
               end if
               if( dimname .eq. name_time ) then
                  time_id = n
                  err = nf_inq_dimlen(ncid,time_id,textr)
                  if (err .ne. NF_NOERR) go to 10
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
            err = nf_inq_varid(ncid,name_lon,id)
            if (err .NE. NF_NOERR) go to 10
            err = nf_get_var_double(ncid,id,met_lon)
            if (err .ne. NF_NOERR) go to 10

            allocate(met_lat(jextr),stat=err)
            if (err /= 0) stop &
                  'open_meteo_file(): Error allocating memory (met_lat)'
            err = nf_inq_varid(ncid,name_lat,id)
            if (err .NE. NF_NOERR) go to 10
            err = nf_get_var_double(ncid,id,met_lat)
            if (err .ne. NF_NOERR) go to 10

            allocate(met_times(textr),stat=err)
            if (err /= 0) stop &
                  'open_meteo_file(): Error allocating memory (met_times)'

            err = nf_inq_varid(ncid,'southpole',id)
            if (err .ne. NF_NOERR) then
               LEVEL4 'Setting southpole to (0,-90,0)'
            else
               err = nf_get_var_double(ncid,id,southpole)
               if (err .ne. NF_NOERR) go to 10
            end if
            LEVEL4 'south pole:'
            LEVEL4 '      lon ',southpole(1)
            LEVEL4 '      lat ',southpole(2)

         end if

         err = nf_inq_dimlen(ncid,time_id,idum)
         if (err .ne. NF_NOERR) go to 10
         if(idum .gt. textr) then
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

         err = nf_inq_varid(ncid,name_time,time_var_id)
         if (err .NE. NF_NOERR) go to 10
         err =  nf_get_att_text(ncid,time_var_id,'units',time_units)
         if (err .NE. NF_NOERR) go to 10
         call string_to_julsecs(time_units,j1,s1)
         err = nf_get_var_real(ncid,time_var_id,met_times)
         if (err .ne. NF_NOERR) go to 10
         call add_secs(j1,s1,nint(met_times(textr)),j2,s2)

         if (in_interval(j1,s1,julianday,secondsofday,j2,s2)) then
            found = .true.
         else
            err = nf_close(ncid)
            if (err .NE. NF_NOERR) go to 10
         end if
      end do
   else
      err = nf_close(ncid)
      if (err .NE. NF_NOERR) go to 10
!     open next file
      read(iunit,*,err=85,end=90) fn
      err = nf_open(fn,NCNOWRIT,ncid)
      if (err .ne. NF_NOERR) go to 10

      err = nf_inq_dimlen(ncid,time_id,idum)
      if (err .ne. NF_NOERR) go to 10
      if(idum .gt. textr) then
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

      err = nf_inq_varid(ncid,name_time,time_var_id)
      if (err .NE. NF_NOERR) go to 10
      err =  nf_get_att_text(ncid,time_var_id,'units',time_units)
      if (err .NE. NF_NOERR) go to 10
      call string_to_julsecs(time_units,j1,s1)
      err = nf_get_var_real(ncid,time_var_id,met_times)
      if (err .ne. NF_NOERR) go to 10

      call add_secs(j1,s1,nint(met_times(textr)),j2,s2)
   end if

   if (found) then
      offset = time_diff(jul0,secs0,j1,s1)
      LEVEL3 'Using meteo from:'
      LEVEL4 trim(fn)
      LEVEL3 'Meteorological offset time ',offset
   else
      FATAL 'Could not find any valid meteo-files'
      stop 'open_meteo_file'
   end if

   return
10 FATAL 'open_meteo_file: ',nf_strerror(err)
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
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   integer         :: i1,i2,istr,j1,j2,jstr
   integer         :: i,j,err
   REALTYPE        :: uu,vv,sinconv,cosconv
!EOP
!-----------------------------------------------------------------------
   include "netcdf.inc"

   if (calc_met) then

      err = nf_get_vara_real(ncid,u10_id,start,edges,wrk)
      if (err .ne. NF_NOERR) go to 10
      if (on_grid) then
         do j=jmin,jmax
            do i=imin,imax
               u10(i,j) = wrk(i,j)
            end do
         end do
      else
         !KBKwrk_dp = _ZERO_
         call copy_var(grid_scan,wrk,wrk_dp)
         call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,u10)
      end if

      err = nf_get_vara_real(ncid,v10_id,start,edges,wrk)
      if (err .ne. NF_NOERR) go to 10
      if (on_grid) then
         do j=jmin,jmax
            do i=imin,imax
               v10(i,j) = wrk(i,j)
            end do
         end do
      else
         !KBKwrk_dp = _ZERO_
         call copy_var(grid_scan,wrk,wrk_dp)
         call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,v10)
      end if

#ifdef SPHERICAL
!     Rotation of wind due to grid convergence
      do j=jmin,jmax
         do i=imin,imax
            if(beta(i,j) .ne. _ZERO_) then
               sinconv=sin(beta(i,j))
               cosconv=cos(beta(i,j))
               uu=u10(i,j)
               vv=v10(i,j)
               u10(i,j)= uu*cosconv+vv*sinconv
               v10(i,j)=-uu*sinconv+vv*cosconv
            end if
         end do
      end do
#else
      if (southpole(1) .ne. 0. .or. southpole(2) .ne. -90.) then
LEVEL0 "rotation of wind due to the combined effect of grid convergence"
LEVEL0 "and rotated meteorological grid is not implemented yet."
         stop "read_data()"
      end if
!     Rotation of wind due to grid convergence
      do j=jmin,jmax
         do i=imin,imax
            if(conv(i,j) .ne. _ZERO_) then
               sinconv=sin(-conv(i,j)*deg2rad)
               cosconv=cos(-conv(i,j)*deg2rad)
               uu=u10(i,j)
               vv=v10(i,j)
               u10(i,j)= uu*cosconv+vv*sinconv
               v10(i,j)=-uu*sinconv+vv*cosconv
            end if
         end do
      end do
#endif

      err = nf_get_vara_real(ncid,airp_id,start,edges,wrk)
      if (err .ne. NF_NOERR) go to 10
      if (on_grid) then
         do j=jmin,jmax
            do i=imin,imax
               airp(i,j) = wrk(i,j)
            end do
         end do
      else
         !KBKwrk_dp = _ZERO_
         call copy_var(grid_scan,wrk,wrk_dp)
         call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,airp)
      end if

      err = nf_get_vara_real(ncid,t2_id,start,edges,wrk)
      if (err .ne. NF_NOERR) go to 10
      if (on_grid) then
         do j=jmin,jmax
            do i=imin,imax
               t2(i,j) = wrk(i,j)
            end do
         end do
      else
         !KBKwrk_dp = _ZERO_
         call copy_var(grid_scan,wrk,wrk_dp)
         call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,t2)
      end if

      err = nf_get_vara_real(ncid,hum_id,start,edges,wrk)
      if (err .ne. NF_NOERR) go to 10
      if (on_grid) then
         do j=jmin,jmax
            do i=imin,imax
               hum(i,j) = wrk(i,j)
            end do
         end do
      else
         !KBKwrk_dp = _ZERO_
         call copy_var(grid_scan,wrk,wrk_dp)
         call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,hum)
      end if

      err = nf_get_vara_real(ncid,tcc_id,start,edges,wrk)
      if (err .ne. NF_NOERR) go to 10
      if (on_grid) then
         do j=jmin,jmax
            do i=imin,imax
               tcc(i,j) = wrk(i,j)
            end do
         end do
      else
         !KBKwrk_dp = _ZERO_
         call copy_var(grid_scan,wrk,wrk_dp)
         call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,tcc)
      end if

   else

      err = nf_get_vara_real(ncid,tausx_id,start,edges,wrk)
      if (err .ne. NF_NOERR) go to 10
      if (point_source) then
         tausx = wrk(1,1)
      else
         do j=jmin,jmax
            do i=imin,imax
               tausx(i,j) = wrk(i,j)
            end do
         end do
      end if

      err = nf_get_vara_real(ncid,tausy_id,start,edges,wrk)
      if (err .ne. NF_NOERR) go to 10
      if (point_source) then
         tausy = wrk(1,1)
      else
         do j=jmin,jmax
            do i=imin,imax
               tausy(i,j) = wrk(i,j)
            end do
         end do
      end if

      err = nf_get_vara_real(ncid,swr_id,start,edges,wrk)
      if (err .ne. NF_NOERR) go to 10
      if (point_source) then
         swr = wrk(1,1)
      else
         do j=jmin,jmax
            do i=imin,imax
               swr(i,j) = wrk(i,j)
            end do
         end do
      end if

      err = nf_get_vara_real(ncid,shf_id,start,edges,wrk)
      if (err .ne. NF_NOERR) go to 10
      if (point_source) then
         shf = wrk(1,1)
      else
         do j=jmin,jmax
            do i=imin,imax
               shf(i,j) = wrk(i,j)
            end do
         end do
      end if

   end if

#ifdef DEBUG
   write(debug,*) 'Leaving read_data()'
   write(debug,*)
#endif
   return
10 FATAL 'read_data: ',nf_strerror(err)
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
   REAL_4B, intent(in)                 :: inf(:,:)
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
   integer         :: i1,i2,istr,j1,j2,jstr
   integer         :: i,j,err
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
