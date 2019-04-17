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
   use domain, only: imin,imax,jmin,jmax,az,lonc,latc,convc,calc_points
   use domain, only: iextr_domain=>iextr,jextr_domain=>jextr
   use domain, only: ill,ihl,jll,jhl,ilg,ihg,jlg,jhg
   use grid_interpol, only: init_grid_interpol,do_grid_interpol
   use grid_interpol, only: to_rotated_lat_lon
   use meteo
   use halo_zones, only : H_TAG,update_2d_halo,wait_halo,periodic_domain
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
   integer         :: sss_id=-1
   integer         :: iextr,jextr,textr,tmax=-1
   integer         :: grid_scan=1
   logical         :: point_source=.false.
   logical         :: rotated_meteo_grid=.false.

   REALTYPE, allocatable     :: met_lon(:),met_lat(:)
   REALTYPE, allocatable      :: met_times(:)

   REALTYPE,dimension(:,:),pointer :: d_airp,airp_input
   REALTYPE,dimension(:,:),pointer :: d_u10,u10_input
   REALTYPE,dimension(:,:),pointer :: d_v10,v10_input
   REALTYPE,dimension(:,:),pointer :: d_tcc,tcc_input
   REALTYPE,dimension(:,:),pointer :: t2_new,d_t2,t2_input
   REALTYPE,dimension(:,:),pointer :: hum_new,d_hum,hum_input
   REALTYPE,dimension(:,:),pointer :: swr_new,d_swr,swr_input
   REALTYPE,dimension(:,:),pointer :: d_precip,precip_input
   REALTYPE,dimension(:,:),pointer :: sst_new,d_sst,sst_input
   REALTYPE,dimension(:,:),pointer :: sss_new,d_sss,sss_input

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

   character(len=10)         :: name_tausx="tausx"
   character(len=10)         :: name_tausy="tausy"
   character(len=10)         :: name_swr="swr"
   character(len=10)         :: name_shf="shf"
   character(len=10)         :: name_sst="sst"
   character(len=10)         :: name_sss="sss"
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
   integer         :: rc,err
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

   if (on_grid) then
      if (iextr_domain.ne.iextr .or. jextr_domain.ne.jextr) then
         call getm_error("init_meteo_input_ncdf()", &
                         "dimensions do not match")
      end if
      il = ilg ; jl = jlg ; ih = ihg ; jh = jhg
   end if

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
      il = 1 ; jl = 1 ; ih = iextr ; jh = jextr
   else
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
   end if


   if ( .not. on_grid ) then

      if (.not. calc_met) then
         stop 'init_meteo_input_ncdf: calc_met=false requires on_grid=true'
      end if

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

      call init_grid_interpol(ill,ihl,jll,jhl,az,  &
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

      il = minval(gridmap(:,:,1),mask=(gridmap(:,:,1).gt.0))
      jl = minval(gridmap(:,:,2),mask=(gridmap(:,:,2).gt.0))
      ih = min( maxval(gridmap(:,:,1))+1 , iextr )
      jh = min( maxval(gridmap(:,:,2))+1 , jextr )

      where( gridmap(:,:,1).gt.0 ) gridmap(:,:,1) = gridmap(:,:,1) - il + 1
      where( gridmap(:,:,2).gt.0 ) gridmap(:,:,2) = gridmap(:,:,2) - jl + 1

   end if

   start(1) = il; start(2) = jl;
   edges(1) = ih-il+1; edges(2) = jh-jl+1;
   edges(3) = 1

   allocate(d_airp(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_meteo_input_ncdf: Error allocating memory (d_airp)'
   d_airp = -9999*_ONE_
   airp_input => d_airp

   if (calc_met) then

      allocate(d_u10(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo_input_ncdf: Error allocating memory (d_u10)'
      d_u10 = -9999*_ONE_
      u10_input => d_u10

      allocate(d_v10(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo_input_ncdf: Error allocating memory (d_v10)'
      d_v10 = -9999*_ONE_
      v10_input => d_v10

      allocate(d_tcc(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo_input_ncdf: Error allocating memory (d_tcc)'
      d_tcc = -9999*_ONE_
      tcc_input => d_tcc

      if (interpolate_meteo) then

         allocate(t2_new(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo_input_ncdf: Error allocating memory (t2_new)'
         t2_new = -9999*_ONE_
         allocate(d_t2(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo_input_ncdf: Error allocating memory (d_t2)'
         d_t2 = -9999*_ONE_
         t2_input => d_t2

         allocate(hum_new(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo_input_ncdf: Error allocating memory (hum_new)'
         hum_new = -9999*_ONE_
         allocate(d_hum(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_meteo_input_ncdf: Error allocating memory (d_hum)'
         d_hum = -9999*_ONE_
         hum_input => d_hum

      else

         t2_input  => t2
         hum_input => hum

      end if

   else

      allocate(swr_new(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo_input_ncdf: Error allocating memory (swr_new)'
      swr_new = -9999*_ONE_
      allocate(d_swr(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo_input_ncdf: Error allocating memory (d_swr)'
      swr_input => d_swr

   end if

   if (fwf_method .eq. 2 .or. fwf_method .eq. 3) then
      allocate(d_precip(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_meteo_input_ncdf: Error allocating memory (d_precip)'
      d_precip = -9999*_ONE_
      precip_input => d_precip
   end if

   if (nudge_sst) then
      allocate(sst_new(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'do_meteo: Error allocating memory (sst_new)'
      sst_new = -9999*_ONE_
      allocate(d_sst(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'do_meteo: Error allocating memory (d_sst)'
      d_sst = -9999*_ONE_
      sst_input => d_sst
   end if

   if (nudge_sss) then
      allocate(sss_new(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'do_meteo: Error allocating memory (sss_new)'
      sss_new = -9999*_ONE_
      allocate(d_sss(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'do_meteo: Error allocating memory (d_sss)'
      d_sss = -9999*_ONE_
      sss_input => d_sss
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
   REALTYPE,dimension(:,:),pointer :: airp_old
   REALTYPE,dimension(:,:),pointer :: u10_old
   REALTYPE,dimension(:,:),pointer :: v10_old
   REALTYPE,dimension(:,:),pointer :: tcc_old
   REALTYPE,dimension(:,:),pointer :: t2_old
   REALTYPE,dimension(:,:),pointer :: hum_old
   REALTYPE,dimension(:,:),pointer :: swr_old
   REALTYPE,dimension(:,:),pointer :: precip_old
   REALTYPE,dimension(:,:),pointer :: sst_old
   REALTYPE,dimension(:,:),pointer :: sss_old

   integer         :: indx
   REALTYPE        :: t,t_minus_t2
   REALTYPE,save   :: deltm1=_ZERO_
   logical, save   :: first=.true.
   integer, save   :: save_n=1
!EOP
!-------------------------------------------------------------------------
#ifdef DEBUG
   integer, save   :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'get_meteo_data_ncdf() # ',Ncall
#endif

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
!            if (indx .gt. 1) then
            if ( t_2 .gt. t ) then
               indx = indx-1
            end if
            t_2 = met_times(indx) - offset
         end if

         if (calc_points.gt.0) call read_data(indx)
         save_n = indx+1

         if (periodic_domain) then
!           need airp(ihl+1,:) and/or airp(:,jhl+1) for periodic velocity points
            call update_2d_halo(airp_input,airp_input,az,imin,jmin,imax,jmax,H_TAG)
            call wait_halo(H_TAG)
         end if

         airp_old=>airp_new;airp_new=>d_airp;d_airp=>airp_old;airp_input=>d_airp
         if (calc_met) then
            u10_old=>u10_new;u10_new=>d_u10;d_u10=>u10_old;u10_input=>d_u10
            v10_old=>v10_new;v10_new=>d_v10;d_v10=>v10_old;v10_input=>d_v10
            tcc_old=>tcc_new;tcc_new=>d_tcc;d_tcc=>tcc_old;tcc_input=>d_tcc
            if (interpolate_meteo) then
               t2_old =>t2_new ;t2_new =>d_t2 ;d_t2 =>t2_old ;t2_input =>d_t2
               hum_old=>hum_new;hum_new=>d_hum;d_hum=>hum_old;hum_input=>d_hum
            end if
         else
            swr_old=>swr_new;swr_new=>d_swr;d_swr=>swr_old;swr_input=>d_swr
         end if
         if (fwf_method.eq.2 .or. fwf_method.eq.3) then
            precip_old=>precip_new;precip_new=>d_precip;d_precip=>precip_old;precip_input=>d_precip
         end if
         if (nudge_sst) then
            sst_old=>sst_new;sst_new=>d_sst;d_sst=>sst_old;sst_input=>d_sst
         end if
         if (nudge_sss) then
            sss_old=>sss_new;sss_new=>d_sss;d_sss=>sss_old;sss_input=>d_sss
         end if

         if (.not. first) then
            d_airp = airp_new - airp_old
            if (calc_met) then
               d_tcc = tcc_new - tcc_old
               d_u10 = u10_new - u10_old
               d_v10 = v10_new - v10_old
               if (interpolate_meteo) then
                  d_t2  = t2_new  - t2_old
                  d_hum = hum_new - hum_old
               end if
            else
               d_swr = swr_new - swr_old
            end if
            if (fwf_method.eq.2 .or. fwf_method.eq.3) then
               d_precip = precip_new - precip_old
            end if
            if (nudge_sst) then
               d_sst = sst_new - sst_old
            end if
            if (nudge_sss) then
               d_sss = sss_new - sss_old
            end if
            deltm1 = _ONE_ / (t_2 - t_1)
         end if

      end if

      t_minus_t2 = t - t_2

      airp = airp_new + d_airp*deltm1*t_minus_t2
      if (calc_met) then
         tcc = tcc_new + d_tcc*deltm1*t_minus_t2
         u10 = u10_new + d_u10*deltm1*t_minus_t2
         v10 = v10_new + d_v10*deltm1*t_minus_t2
         if (interpolate_meteo) then
            t2  = t2_new  + d_t2 *deltm1*t_minus_t2
            hum = hum_new + d_hum*deltm1*t_minus_t2
         end if
      else
         swr = swr_new + d_swr*deltm1*t_minus_t2
      end if
      if (fwf_method.eq.2 .or. fwf_method.eq.3) then
         precip = precip_new + d_precip*deltm1*t_minus_t2
      end if
      if (nudge_sst) then
         sst = sst_new + d_sst*deltm1*t_minus_t2
      end if
      if (nudge_sss) then
         sss = sss_new + d_sss*deltm1*t_minus_t2
      end if

      first = .false.

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
   character(len=19)         :: tbuf
   integer         :: junit,sunit,j1,s1,j2,s2
   integer         :: n,err,idum
   logical         :: first=.true.
   integer, save   :: lon_id=-1,lat_id=-1,time_id=-1,id=-1
   integer, save   :: time_var_id=-1
   character(len=256) :: dimname
   logical         :: have_southpole
!
!EOP
!-------------------------------------------------------------------------
#ifdef DEBUG
   integer, save   :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'open_meteo_file() # ',Ncall
#endif

      if (first) open(iunit,file=meteo_file,status='old',action='read',err=80)

      do
         if (.not. first) then
            err = nf90_close(ncid)
            if (err .NE. NF90_NOERR) go to 10
         end if
         read(iunit,*,err=85,end=90) fn
         LEVEL3 'Trying meteo from:'
         LEVEL4 trim(fn)
         err = nf90_open(fn,NF90_NOWRITE,ncid)
         if (err .ne. NF90_NOERR) go to 10

         err = nf90_inq_dimid(ncid,name_time,time_id)
         if (err .NE. NF90_NOERR) go to 10

         err = nf90_inquire_dimension(ncid,time_id,len=idum)
         if (err .ne. NF90_NOERR) go to 10
         if(idum .gt. size(met_times)) then
            if (allocated(met_times)) then
            deallocate(met_times,stat=err)
            if (err /= 0) stop      &
               'open_meteo_file(): Error de-allocating memory (met_times)'
            end if
            allocate(met_times(idum),stat=err)
            if (err /= 0) stop &
               'open_meteo_file(): Error allocating memory (met_times)'
         end if
         textr = idum
         LEVEL4 'time_id --> ',time_id,', len = ',textr
!        if (tmax .lt. 0) tmax=textr
         tmax=textr

         err = nf90_inq_varid(ncid,name_time,time_var_id)
         if (err .NE. NF90_NOERR) go to 10
         err =  nf90_get_att(ncid,time_var_id,'units',time_units)
         if (err .NE. NF90_NOERR) go to 10
         call string_to_julsecs(time_units,junit,sunit)
         err = nf90_get_var(ncid,time_var_id,met_times(1:textr))
         if (err .ne. NF90_NOERR) go to 10

         call add_secs(junit,sunit,nint(met_times(    1)),j1,s1)
         call write_time_string(j1,s1,tbuf)
         LEVEL4 'Datafile starts:   ',tbuf

         call add_secs(junit,sunit,nint(met_times(textr)),j2,s2)
         call write_time_string(j2,s2,tbuf)
         LEVEL4 'Datafile ends  :   ',tbuf

!        KK-TODO: if not first, julianday and secondsofday lag one timestep behind!!!
         if (in_interval(j1,s1,julianday,secondsofday,j2,s2)) exit

         if ( time_diff(julianday,secondsofday,junit,sunit) .lt. met_times(1) ) then
            if (first) then
               FATAL 'Model simulation starts before available meteo data'
               call write_time_string()
               FATAL 'Simulation starts: ',timestr
               call write_time_string(j1,s1,tbuf)
               FATAL 'Datafile starts:   ',tbuf
               stop 'open_meteo_file()'
            else
               LEVEL4 'WARNING: possible time gap - check start of Datafile!'
               call write_time_string()
               LEVEL4 'Simulation clock: ',timestr
               call write_time_string(j1,s1,tbuf)
               LEVEL4 'Datafile starts:   ',tbuf
               exit
            end if
         end if
         if (.not. first) then
!           KK-TODO: Or should we allow to cycle?
            FATAL 'Datafile does not contain new records'
            stop 'open_meteo_file()'
         end if

      end do

         if (first) then
            first = .false.
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
            end do
            if(lon_id .eq. -1) then
               FATAL 'could not find longitude coordinate in meteo file'
               stop 'open_meteo_file()'
            end if
            if(lat_id .eq. -1) then
               FATAL 'could not find latitude coordinate in meteo file'
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

!           first we check for CF compatible grid_mapping_name
            err = nf90_inq_varid(ncid,'rotated_pole',id)
            if (err .eq. NF90_NOERR) then
               LEVEL4 'Reading CF-compliant rotated grid specification'
               err = nf90_get_att(ncid,id, &
                                  'grid_north_pole_latitude',southpole(1))
               if (err .ne. NF90_NOERR) go to 10
               err = nf90_get_att(ncid,id, &
                                  'grid_north_pole_longitude',southpole(2))
#if 0
STDERR 'Inside rotated_pole'
STDERR 'grid_north_pole_latitude ',southpole(1)
STDERR 'grid_north_pole_longitude ',southpole(2)
#endif
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
#if 0
STDERR 'After transformation:'
STDERR 'grid_north_pole_latitude ',southpole(1)
STDERR 'grid_north_pole_longitude ',southpole(2)
#endif
               southpole(3) = _ZERO_
               have_southpole = .true.
               rotated_meteo_grid = .true.
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
                  rotated_meteo_grid = .true.
               end if
            end if
            if (rotated_meteo_grid) then
               LEVEL4 'south pole:'
!              changed indices - kb 2014-12-15
               LEVEL4 '      lon ',southpole(2)
               LEVEL4 '      lat ',southpole(1)
            end if
         end if

      airp_id = ncdf_meteo_inq_varid(ncid,name_airp)

      if (fwf_method .eq. 2) then
         evap_id   = ncdf_meteo_inq_varid(ncid,name_evap  ,required=.true.)
      end if
      if (fwf_method .eq. 2 .or. fwf_method .eq. 3) then
         precip_id = ncdf_meteo_inq_varid(ncid,name_precip,required=.true.)
      end if

      if (calc_met) then
         u10_id  = ncdf_meteo_inq_varid(ncid,name_u10 )
         v10_id  = ncdf_meteo_inq_varid(ncid,name_v10 )
         t2_id   = ncdf_meteo_inq_varid(ncid,name_t2  )
         tcc_id  = ncdf_meteo_inq_varid(ncid,name_tcc )

         hum_id = -1
         err = nf90_inq_varid(ncid,name_hum1,hum_id)
         if (err .NE. NF90_NOERR) then
            err = nf90_inq_varid(ncid,name_hum2,hum_id)
            if (err .NE. NF90_NOERR) then
               err = nf90_inq_varid(ncid,name_hum3,hum_id)
               if (err .NE. NF90_NOERR) then
                  hum_id = -1
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

         if (t2_id.eq.-1 .and. tcc_id.eq.-1 .and. hum_id.eq.-1) then
            LEVEL4 'WARNING: Incomplete input data for bulk formulas.'
            LEVEL4 '         Stresses will be calculated with constant drag coefficient.'
            constant_cd = .true.
         else if (airp_id.eq.-1 .or. t2_id.eq.-1 .or. tcc_id.eq.-1 .or. hum_id.eq.-1) then
            FATAL 'Incomplete input data for bulk formulas.'
            stop 'init_meteo_input_ncdf()'
         end if

      else
         tausx_id = ncdf_meteo_inq_varid(ncid,name_tausx)
         tausy_id = ncdf_meteo_inq_varid(ncid,name_tausy)
         swr_id   = ncdf_meteo_inq_varid(ncid,name_swr  )
         shf_id   = ncdf_meteo_inq_varid(ncid,name_shf  )
      end if

      if (nudge_sst) then
         sst_id = ncdf_meteo_inq_varid(ncid,name_sst,required=.true.)
      end if

      if (nudge_sss) then
         sss_id = ncdf_meteo_inq_varid(ncid,name_sss,required=.true.)
      end if

      offset = time_diff(jul0,secs0,junit,sunit)
      LEVEL3 'Using meteo from:'
      LEVEL4 trim(fn)
      LEVEL3 'Meteorological offset time ',offset


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
   subroutine read_data(indx)
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
! !INPUT PARAMETERS:
   integer,intent(in) :: indx
!
! !LOCAL VARIABLES:
   integer         :: i,j,err
   REALTYPE,dimension(edges(1),edges(2)) :: wrk!,wrk_dp
   REALTYPE        :: angle,uu,vv,sinconv,cosconv
!EOP
!-----------------------------------------------------------------------

   call write_time_string()
   LEVEL3 timestr,': reading meteo data ...',indx
   start(3) = indx

   if (airp_id .gt. 0) then
      err = nf90_get_var(ncid,airp_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (on_grid) then
         if (point_source) then
            airp_input = wrk(1,1)
         else
            airp_input(ill:ihl,jll:jhl) = wrk
         end if
      else
         !KBKwrk_dp = _ZERO_
         !call copy_var(grid_scan,wrk,wrk_dp)
         !call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,airp_input)
         call flip_var(wrk)
         call do_grid_interpol(az,wrk,gridmap,ti,ui,airp_input)
      end if
   end if

   if (evap_id .ge. 0) then
      err = nf90_get_var(ncid,evap_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (on_grid) then
         if (point_source) then
            evap_input = wrk(1,1)
         else
            evap_input(ill:ihl,jll:jhl) = wrk
         end if
      else
         !call copy_var(grid_scan,wrk,wrk_dp)
         !call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,evap_input)
         call flip_var(wrk)
         call do_grid_interpol(az,wrk,gridmap,ti,ui,evap_input)
      end if
      if (evap_factor .ne. _ONE_) then
         evap_input = evap_input * evap_factor
      end if
   end if

   if (precip_id .gt. 0) then
      err = nf90_get_var(ncid,precip_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (on_grid) then
         if (point_source) then
            precip_input = wrk(1,1)
         else
            precip_input(ill:ihl,jll:jhl) = wrk
         end if
      else
         !call copy_var(grid_scan,wrk,wrk_dp)
         !call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,precip_input)
         call flip_var(wrk)
         call do_grid_interpol(az,wrk,gridmap,ti,ui,precip_input)
      end if
      if (precip_factor .ne. _ONE_) then
         precip_input = precip_input * precip_factor
      end if
   end if

   if (calc_met) then

   if (u10_id .gt. 0) then
      err = nf90_get_var(ncid,u10_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (on_grid) then
         if (point_source) then
            u10_input = wrk(1,1)
         else
            u10_input(ill:ihl,jll:jhl) = wrk
         end if
      else
         !KBKwrk_dp = _ZERO_
         !call copy_var(grid_scan,wrk,wrk_dp)
         !call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,u10_input)
         call flip_var(wrk)
         call do_grid_interpol(az,wrk,gridmap,ti,ui,u10_input)
      end if
   end if

   if (v10_id .gt. 0) then
      err = nf90_get_var(ncid,v10_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (on_grid) then
         if (point_source) then
            v10_input = wrk(1,1)
         else
            v10_input(ill:ihl,jll:jhl) = wrk
         end if
      else
         !KBKwrk_dp = _ZERO_
         !call copy_var(grid_scan,wrk,wrk_dp)
         !call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,v10_input)
         call flip_var(wrk)
         call do_grid_interpol(az,wrk,gridmap,ti,ui,v10_input)
      end if
   end if

!     Rotation of wind due to the combined effect of possible rotation of
!     meteorological grid and possible hydrodynamic grid convergence
!     (cartesian and curvi-linear grids where conv <> 0.)
      do j=jll,jhl
         do i=ill,ihl
!KBK            angle=-convc(i,j)*deg2rad
!KBK            angle=beta(i,j)
            angle=beta(i,j)-convc(i,j)*deg2rad
            if(angle .ne. _ZERO_) then
               sinconv=sin(angle)
               cosconv=cos(angle)
               uu=u10_input(i,j)
               vv=v10_input(i,j)
               u10_input(i,j)= uu*cosconv+vv*sinconv
               v10_input(i,j)=-uu*sinconv+vv*cosconv
            end if
         end do
      end do

      if (wind_factor .ne. _ONE_) then
         u10_input = u10_input * wind_factor
         v10_input = v10_input * wind_factor
      end if

   if (t2_id .gt. 0) then
      err = nf90_get_var(ncid,t2_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (on_grid) then
         if (point_source) then
            t2_input = wrk(1,1)
         else
            t2_input(ill:ihl,jll:jhl) = wrk
         end if
      else
         !KBKwrk_dp = _ZERO_
         !call copy_var(grid_scan,wrk,wrk_dp)
         !call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,t2_input)
         call flip_var(wrk)
         call do_grid_interpol(az,wrk,gridmap,ti,ui,t2_input)
      end if
   end if

   if (hum_id .gt. 0) then
      err = nf90_get_var(ncid,hum_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (on_grid) then
         if (point_source) then
            hum_input = wrk(1,1)
         else
            hum_input(ill:ihl,jll:jhl) = wrk
         end if
      else
         !KBKwrk_dp = _ZERO_
         !call copy_var(grid_scan,wrk,wrk_dp)
         !call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,hum_input)
         call flip_var(wrk)
         call do_grid_interpol(az,wrk,gridmap,ti,ui,hum_input)
      end if
   end if

   if (tcc_id .gt. 0) then
      err = nf90_get_var(ncid,tcc_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (on_grid) then
         if (point_source) then
            tcc_input = wrk(1,1)
         else
            tcc_input(ill:ihl,jll:jhl) = wrk
         end if
      else
         !KBKwrk_dp = _ZERO_
         !call copy_var(grid_scan,wrk,wrk_dp)
         !call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,tcc_input)
         call flip_var(wrk)
         call do_grid_interpol(az,wrk,gridmap,ti,ui,tcc_input)
      end if
   end if

   else

      if (tausx_id .gt. 0) then
         err = nf90_get_var(ncid,tausx_id,wrk,start,edges)
         if (err .ne. NF90_NOERR) go to 10
         if (point_source) then
            tausx_input = wrk(1,1)
         else
            tausx_input(ill:ihl,jll:jhl) = wrk
         end if
      end if
      if (tausy_id .gt. 0) then
         err = nf90_get_var(ncid,tausy_id,wrk,start,edges)
         if (err .ne. NF90_NOERR) go to 10
         if (point_source) then
            tausy_input = wrk(1,1)
         else
            tausy_input(ill:ihl,jll:jhl) = wrk
         end if
      end if
      if (swr_id .gt. 0) then
         err = nf90_get_var(ncid,swr_id,wrk,start,edges)
         if (err .ne. NF90_NOERR) go to 10
         if (point_source) then
            swr_input = wrk(1,1)
         else
            swr_input(ill:ihl,jll:jhl) = wrk
         end if
      end if
      if (shf_id .gt. 0) then
         err = nf90_get_var(ncid,shf_id,wrk,start,edges)
         if (err .ne. NF90_NOERR) go to 10
         if (point_source) then
            shf_input = wrk(1,1)
         else
            shf_input(ill:ihl,jll:jhl) = wrk
         end if
      end if

   end if

   if (sst_id .gt. 0) then
      err = nf90_get_var(ncid,sst_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (on_grid) then
         if (point_source) then
            sst_input = wrk(1,1)
         else
            sst_input(ill:ihl,jll:jhl) = wrk
         end if
      else if (calc_met) then
         !call copy_var(grid_scan,wrk,wrk_dp)
         !call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,sst_input)
         call flip_var(wrk)
         call do_grid_interpol(az,wrk,gridmap,ti,ui,sst_input)
      end if
   end if

   if (sss_id .gt. 0) then
      err = nf90_get_var(ncid,sss_id,wrk,start,edges)
      if (err .ne. NF90_NOERR) go to 10
      if (on_grid) then
         if (point_source) then
            sss_input = wrk(1,1)
         else
            sss_input(ill:ihl,jll:jhl) = wrk
         end if
      else if (calc_met) then
         !call copy_var(grid_scan,wrk,wrk_dp)
         !call do_grid_interpol(az,wrk_dp,gridmap,ti,ui,sss_input)
         call flip_var(wrk)
         call do_grid_interpol(az,wrk,gridmap,ti,ui,sss_input)
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
! !IROUTINE: flip_var -
!
! !INTERFACE:
   subroutine flip_var(var)
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout) :: var(edges(1),edges(2))
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------

   select case (grid_scan)
      case (0)
         var = var(:,edges(2):1:-1)
   end select
   return
   end subroutine flip_var
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
   integer         :: iextr,jextr
!EOP
!-----------------------------------------------------------------------

   iextr = edges(1) ; jextr = edges(2)

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
!BOP
!
! !IROUTINE: ncdf_meteo_inq_varid -
!
! !INTERFACE:
   integer function ncdf_meteo_inq_varid(ncid,varname,required)
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   integer         , intent(in)        :: ncid
   character(len=*), intent(in)        :: varname
   logical         , intent(in), optional :: required
!
! !LOCAL VARIABLES:
   integer         :: varid,err
   logical         :: break_on_missing
!EOP
!-----------------------------------------------------------------------

   if ( present(required) ) then
      break_on_missing = required
   else
      break_on_missing = .false.
   end if

   LEVEL4 ' ... checking variable ',trim(varname)
   err = nf90_inq_varid(ncid,varname,varid)
   if (err .NE. NF90_NOERR) then
      if ( break_on_missing ) then
         FATAL 'ncdf_meteo_inq_varid: ',nf90_strerror(err)
         stop 'ncdf_meteo_inq_varid()'
      else
         LEVEL4 '       missing - continue with '//trim(varname)//'=0'
         varid = -1
      end if
   end if

   ncdf_meteo_inq_varid = varid
   end function ncdf_meteo_inq_varid
!EOC
!-----------------------------------------------------------------------

   end module ncdf_meteo

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
