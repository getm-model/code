#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: domain - sets up the calculation domain.
!
! !INTERFACE:
   module domain
!
! !DESCRIPTION:
!  This module provides all variables related to the bathymetry and
!  model grid. The public subroutine $init\_domain()$ is called once
!  and upon successful completion the bathymetry has been read and
!  optionally modified, the calculation masks have been setup and
!  all grid related variables have been initialised.\newline
!  The $domain$-module depends on another module doing the actual
!  reading of variables from files. This is provided through the
!  generic subroutine $read\_topo\_file$. This subroutine takes two
!  parameters - 1) a fileformat and 2) a filename. Adding a new
!  input file format is thus straight forward and can be done
!  without any changes to $domain$.
!  Public variables defined in this module is used through out the
!  code.
!
! !USES:
   use exceptions
   use halo_zones,     only: update_2d_halo,wait_halo
   use halo_zones,     only: H_TAG,U_TAG,V_TAG
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   integer                             :: bathy_format   = NETCDF

   integer                             :: grid_type      = 1
   integer                             :: vert_cord      = 1
   integer                             :: iextr_topo,jextr_topo
!  global index range
   integer                             :: ilg=-1,ihg=-1,jlg=-1,jhg=-1
!  local index range
   integer                             :: ill=-1,ihl=-1,jll=-1,jhl=-1

   logical                             :: have_lonlat    = .true.
   logical                             :: have_xy        = .true.

   REALTYPE                            :: rearth

   REALTYPE                            :: maxdepth       = -1.
   REALTYPE                            :: ddu            = -_ONE_
   REALTYPE                            :: ddl            = -_ONE_
   REALTYPE                            :: d_gamma        = 20.
   logical                             :: gamma_surf     = .true.
   REALTYPE, allocatable, dimension(:) :: ga

   integer                             :: NWB=-1,NNB=-1,NEB=-1,NSB=-1,NOB
   integer                             :: calc_points
   logical                             :: openbdy        = .false.

   REALTYPE                            :: Hland=-10.0
   REALTYPE                            :: min_depth,crit_depth

   REALTYPE                            :: longitude      = _ZERO_
   REALTYPE                            :: latitude       = _ZERO_
   logical                             :: f_plane        = .true.

#ifdef STATIC
#include "static_domain.h"
#else
#include "dynamic_declarations_domain.h"
#endif
   integer                             :: nsbv

   integer                             :: ioff=0,joff=0
   integer, dimension(:), allocatable  :: bdy_2d_type
   integer, dimension(:), allocatable  :: bdy_3d_type
   integer, dimension(:), allocatable  :: wi,wfj,wlj
   integer, dimension(:), allocatable  :: nj,nfi,nli
   integer, dimension(:), allocatable  :: ei,efj,elj
   integer, dimension(:), allocatable  :: sj,sfi,sli
   integer, allocatable                :: bdy_index(:),bdy_map(:,:)
   logical                             :: have_boundaries=.false.

   character(len=64)                   :: bdy_2d_desc(5)
   logical                             :: need_2d_bdy_elev = .false.
   logical                             :: need_2d_bdy_u    = .false.
   logical                             :: need_2d_bdy_v    = .false.

   REALTYPE                            :: cori= _ZERO_

!  method for specifying bottom roughness (0=const, 1=from topo.nc)
   integer                             :: z0_method=0
   REALTYPE                            :: z0_const=0.01d0

! !DEFINED PARAMETERS:
   integer,           parameter        :: INNER          = 1
   REALTYPE, private, parameter        :: pi             = 3.141592654
   REALTYPE, private, parameter        :: deg2rad        = pi/180.
   REALTYPE, private, parameter        :: omega          = 2.*pi/86400.
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   REALTYPE, parameter                  :: rearth_default = 6378815.
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_domain() - initialise the computational domain

!
! !INTERFACE:
   subroutine init_domain(input_dir)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  This routine is responsible for setting up the bathymetry and
!  the grid information.\newline
!  The following steps are done in $init\_domain()$:
!  \begin{itemize}
!     \item[1:] partition of the calculation domain - important for
!               parallel runs
!     \item[2:] reading bathymetry and grid information through the
!               generic subroutine $read\_topo\_file$
!     \item[3:] optionally set minimum depth in regions
!     \item[4:] optionally adjust the depth in regions
!     \item[5:] optionally adjust the depth in regions
!     \item[6:] calculate the mask for T-points
!     \item[7:] optionally adjust the mask in regions
!     \item[8:] read boundary information and adjust masks
!     \item[9:] calculate masks for U-, V- and X-points
!     \item[10:] calculate additional grid-information - like $latu$ and
!                $latv$
!     \item[11:] calculate metrics - i.e. all necessary grid-spacings
!     \item[12:] calculate Coriolis parameter - can be constant or
!                spatially varying
!  \end{itemize}
!
! !INPUT/OUTPUT PARAMETERS:
   character(len=*)                    :: input_dir
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
   integer                   :: rc
   integer                   :: np,sz
   integer                   :: i,j,n
   integer                   :: kdum
   character(len=PATH_MAX)   :: bathymetry               = 'topo.nc'
   integer                   :: vel_depth_method=0
   character(len=PATH_MAX)   :: bdyinfofile              = 'bdyinfo.dat'
   character(len=PATH_MAX)   :: min_depth_file           = 'minimum_depth.dat'
   character(len=PATH_MAX)   :: bathymetry_adjust_file   = 'bathymetry.adjust'
   character(len=PATH_MAX)   :: mask_adjust_file         = 'mask.adjust'
   integer                   :: il=-1,ih=-1,jl=-1,jh=-1
   integer                   :: iskipl,jskipl
   namelist /domain/ &
             vert_cord,maxdepth,                               &
             bathy_format,bathymetry,vel_depth_method,         &
             longitude,latitude,f_plane,openbdy,bdyinfofile,   &
             crit_depth,min_depth,kdum,ddu,ddl,                &
             d_gamma,gamma_surf,il,ih,jl,jh,z0_method,z0_const
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_domain()'
#endif

   bdy_2d_desc(ZERO_GRADIENT)           = "Zero gradient"
   bdy_2d_desc(SOMMERFELD)              = "Sommerfeld rad."
   bdy_2d_desc(CLAMPED)                 = "Clamped"
   bdy_2d_desc(FLATHER_ELEV)            = "Flather (elev)"
   bdy_2d_desc(FLATHER_VEL)             = "Flather (vel)"

   LEVEL1 'init_domain'

!  Read domain specific things from the namelist.
   read(NAMLST,domain)

   if (crit_depth .lt. 2.5*min_depth)  then
      stop 'crit_depth must be larger than 2.5 time min_depth'
   end if

   call open_topo_file(bathy_format,bathymetry,grid_type,iextr_topo,jextr_topo)

!   LEVEL2 'grid_type: ',grid_type
   select case (grid_type)
      case(1)
         LEVEL2 'using cartesian grid.'
      case(2)
         LEVEL2 'using spherical grid.'
      case(3)
         LEVEL2 'using plane curvilinear grid.'
      case(4)
         LEVEL2 'using spherical curvilinear grid.'
         call getm_error("init_domain()","grid_type=4 has not been properly implemented yet.")
      case default
         call getm_error("init_domain()","Invalid grid type. Choose grid_type=1-4.")
   end select

#ifdef STATIC
   if (iextr_topo.ne.iextr .or. jextr_topo.ne.jextr) then
      FATAL "init_domain: inconsistent grid size:"
      FATAL "   expected ",iextr     ," by ",jextr
      FATAL "   read     ",iextr_topo," by ",jextr_topo
      call getm_error("init_domain()","inconsistent grid size")
   end if
#else
   iextr = iextr_topo
   jextr = jextr_topo
   kmax=kdum
#endif
   LEVEL3 'iextr, jextr: ',iextr,jextr

!  prepare parallel run
   call part_domain()
   il=imin ; ih=imax ; jl=jmin ; jh=jmax

#ifndef STATIC
#include "dynamic_allocations_domain.h"
#endif

!  GLOBAL index range
   ilg = max(imin-HALO+ioff,1); ihg = min(imax+HALO+ioff,iextr)
   jlg = max(jmin-HALO+joff,1); jhg = min(jmax+HALO+joff,jextr)
   iskipl= ilg - (imin-HALO+ioff)
   jskipl= jlg - (jmin-HALO+joff)

!  LOCAL index range
!  (different from GLOBAL range only for parallel runs)
   ill = imin-HALO+iskipl; jll = jmin-HALO+jskipl;
   ihl = ihg-ilg+ill;      jhl = jhg-jlg+jll;

   call read_topo_file(bathy_format,grid_type)

   select case (vert_cord)
      case(_SIGMA_COORDS_)
         LEVEL2 'Using sigma coordinates'
      case(_Z_COORDS_)
         LEVEL2 'Using z-level coordinates'
      case(_GENERAL_COORDS_)
         LEVEL2 'Using general vertical coordinates'
      case (_HYBRID_COORDS_) ! hybrid vertical coordinates
         LEVEL2 'using hybrid vertical coordinates'
         STDERR 'domain: hybrid_coordinates not coded yet'
         stop
      case (_ADAPTIVE_COORDS_) ! adaptive vertical coordinates
         LEVEL2 'using adaptive vertical coordinates'
      case default
         call getm_error("init_domain()", &
                         "A non valid vertical coordinate system has been chosen");
   end select

!  Calculation masks
!  Do we want to set a minimum depth for certain regions
   call set_min_depth(trim(input_dir) // min_depth_file)

!  Do we want to do adjust the bathymetry
   call adjust_bathymetry(trim(input_dir) // bathymetry_adjust_file)

   call update_2d_halo(H,H,az,imin,jmin,imax,jmax,H_TAG,mirror=.false.)
   call wait_halo(H_TAG)

   az = 0
   where (H .gt. Hland+SMALL)
      az=1
   end where

   call update_2d_halo(H,H,az,imin,jmin,imax,jmax,H_TAG,mirror=.true.)
   call wait_halo(H_TAG)

   call uv_depths(vel_depth_method)

!  Reads boundary location information
   if (openbdy) then
      call bdy_spec(trim(input_dir) // bdyinfofile)
      call print_bdy('Global Boundary Information')
      call have_bdy()
      call print_bdy('Local Boundary Information')
   end if

#define BOUNDARY_POINT 2
!  western boundary - at present elev only
   do n=1,NWB
      az(wi(n),wfj(n):wlj(n)) = BOUNDARY_POINT
   end do
!  northern boundary - at present elev only
   do n=1,NNB
      az(nfi(n):nli(n),nj(n)) = BOUNDARY_POINT
   end do
!  easter boundary - at present elev only
   do n=1,NEB
      az(ei(n),efj(n):elj(n)) = BOUNDARY_POINT
   end do
!  southern boundary - at present elev only
   do n=1,NSB
      az(sfi(n):sli(n),sj(n)) = BOUNDARY_POINT
   end do
#undef BOUNDARY_POINT

!  Do we want to further adjust the mask
   call adjust_mask(trim(input_dir) // mask_adjust_file)

   mask = _ONE_*az
   call update_2d_halo(mask,mask,az,imin,jmin,imax,jmax,H_TAG,mirror=.false.)
   call wait_halo(H_TAG)
   az=mask

!  mask for U-points
   mask=0
   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO-1
         if (az(i,j) .eq. 1 .and. az(i+1,j) .eq. 1) then
            mask(i,j)=1
         end if
         if ((az(i,j) .eq. 1 .and. az(i+1,j) .eq. 2).or.    &
             (az(i,j) .eq. 2 .and. az(i+1,j) .eq. 1)) then
            mask(i,j)=2
         end if
         if (az(i,j) .eq. 2 .and. az(i+1,j) .eq. 2) then
            mask(i,j)=3
         end if
      end do
   end do
   call update_2d_halo(mask,mask,az,imin,jmin,imax,jmax,H_TAG,mirror=.false.)
   call wait_halo(H_TAG)
   au = mask
   
!  mask for V-points
   mask=_ZERO_
   do j=jmin-HALO,jmax+HALO-1
      do i=imin-HALO,imax+HALO
         if (az(i,j) .eq. 1 .and. az(i,j+1) .eq. 1) then
            mask(i,j)=1
         end if
         if ((az(i,j) .eq. 1 .and. az(i,j+1) .eq. 2).or.    &
             (az(i,j) .eq. 2 .and. az(i,j+1) .eq. 1)) then
            mask(i,j)=2
         end if
         if (az(i,j) .eq. 2 .and. az(i,j+1) .eq. 2) then
            mask(i,j)=3
         end if
      end do
   end do
   call update_2d_halo(mask,mask,az,imin,jmin,imax,jmax,H_TAG,mirror=.false.)
   call wait_halo(H_TAG)
   av = mask

!  mask for X-points
   mask=0
   do j=jmin-HALO,jmax+HALO-1
      do i=imin-HALO,imax+HALO-1
         if (az(i  ,j) .ge. 1 .and. az(i  ,j+1) .ge. 1 .and.    &
             az(i+1,j) .ge. 1 .and. az(i+1,j+1) .ge. 1) then
            mask(i,j)=1
         end if
      end do
   end do
   call update_2d_halo(mask,mask,az,imin,jmin,imax,jmax,H_TAG)
   call wait_halo(H_TAG)
   ax = mask

   call mirror_bdy_2d(H,H_TAG)

!  Interpolate (latx,lonx) and (xx,yx) to the u, v, and T-points.
   call x2uvc()

!  calculate the metric coeff.
   call metric()

   if ( .not. have_lonlat ) then
      LEVEL2 "Setting constant longitude (swr) - lon = ",longitude
      lonc = longitude
      LEVEL2 "Setting constant latitude (meteo) - lat = ",latitude
      latc = latitude
   end if

!  Compute Coriolis parameter
   if (f_plane) then
      LEVEL2 "Assuming constant Coriolis parameter - lat = ",latitude
      cori = 2.*omega*sin(deg2rad*latitude)
      coru = cori
      corv = cori
   else
      if (have_lonlat) then
         LEVEL2 "Computing spatially varying Coriolis parameter"

         coru = 2.*omega*sin(deg2rad*latu)
         corv = 2.*omega*sin(deg2rad*latv)
      else
         call getm_error("init_domain()",   &
              "f_plane=.false. only possible if (lat,lon) exist.")
      endif
   endif

   call update_2d_halo(coru,coru,au,imin,jmin,imax,jmax,U_TAG)
   call wait_halo(U_TAG)

   call update_2d_halo(corv,corv,av,imin,jmin,imax,jmax,V_TAG)
   call wait_halo(V_TAG)

   select case (z0_method)
      case(0)
         LEVEL2 'Using constant bottom roughness'
      case(1)
         LEVEL2 'Using space varying bottom roughness'
         LEVEL2 '..  will read z0 from the topo file ..'
         call update_2d_halo(z0,z0,az,imin,jmin,imax,jmax,H_TAG)
         call wait_halo(H_TAG)
      case default
         call getm_error("init_domain()", &
                         "A non valid z0 method has been chosen");
   end select

#ifdef DEBUG
   STDERR 'az'
   call print_mask(az)
   STDERR 'au'
   call print_mask(au)
   STDERR 'av'
   call print_mask(av)
#endif

   np = count(az(imin:imax,jmin:jmax) .gt. 0)
   sz = (imax-imin+1)*(jmax-jmin+1)
   LEVEL2 'Dimensions: ',imin,':',imax,',',jmin,':',jmax,',',0,':',kmax
   LEVEL2 '# waterpoints = ',np,' of ',sz

   calc_points = np

#ifdef DEBUG
   write(debug,*) 'Leaving init_domain()'
   write(debug,*)
#endif
   return
   end subroutine init_domain
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: x2uvc() - interpolate grid-points
!
! !INTERFACE:
   subroutine x2uvc()
   IMPLICIT NONE
!
! !DESCRIPTION:
!  This routine interpolates (latx,lonx), (xx,yx), and convx to the
!  u-points, v-points, and the central T-points. The data at the T-points
!  are only updated from values of the X-points if the logical flags {\tt updateXYC},
!  {\tt updateXYC}, and {\tt updateXYC} are {\tt .true.}. This is not necessary
!  if data at the T-points have been read from the topo input file.
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
! !LOCAL VARIABLES:
   integer                   :: i,j,n
   REALTYPE                  :: x
!EOP
!------------------------------------------------------------------------
!BOC

   select case (grid_type)

      case(1)

         if ( have_lonlat ) then

!           latu and latv are only needed for active momentum points
            do j=jll,jhl-1
               do i=ill,ihl-1

                  if ( au(i,j) .eq. 1 .or. au(i,j) .eq. 2 ) then

                     latu(i,j) = ( latc(i  ,j) + latc(i+1,j+1) &
                                 + latc(i+1,j) + latc(i  ,j+1) ) / 4
                  end if

                  if ( av(i,j) .eq. 1 .or. av(i,j) .eq. 2 ) then
                     latv(i,j) = ( latc(i  ,j) + latc(i+1,j+1) &
                                 + latc(i+1,j) + latc(i  ,j+1) ) / 4
                  end if

               end do
            end do

! this is just a check and can be deleted if nobody experiences problems
#if 1
            if ( joff+jhl .eq. jextr ) then
               do i=ill,ihl
                  if ( au(i,jhl) .eq. 1 .or. au(i,j) .eq. 2 ) then
                     latu(i,jhl) = 1000
                     LEVEL0 'x2uvc() - warning: latu is set to illegal value'
                     LEVEL0 'please report the problem on getm-users'
                     stop
                  end if
               end do
            end if

            if ( ioff+ihl .eq. iextr ) then
               do j=jll,jhl
                  if ( av(ihl,j) .eq. 1 .or. av(i,j) .eq. 2 ) then
                     latv(i,jhl) = 1000
                     LEVEL0 'x2uvc() - warning: latv is set to illegal value'
                     LEVEL0 'please report the problem on getm-users'
                     stop
                  end if
               end do
            end if
#endif

         end if


      case(2)

!        we need latx to calculate dxv - utilize equidistance
         latx(ill:ihl,jll-1) = latc(ill:ihl,jll) - dlat/2.
STDERR ill,jll,dlat/2.
STDERR latc(1,1),latx(1,0)
!stop
         n=1
         do j=jll,jhl
            latx(ill:ihl,j) = latx(ill:ihl,jll-1) + n*dlat
            n=n+1
         end do

         latu = latc
         latv(ill:ihl,jll:jhl) = latx(ill:ihl,jll:jhl)

      case(3)

         do j=jll,jhl
            do i=ill,ihl
               xu(i,j)   = ( xx(i,j) +   xx(i,j-1) ) / 2
               yu(i,j)   = ( yx(i,j) +   yx(i,j-1) ) / 2

               xv(i,j)   = ( xx(i,j) +   xx(i-1,j) ) / 2
               yv(i,j)   = ( yx(i,j) +   yx(i-1,j) ) / 2
            end do
         end do

         do j=jll,jhl
            do i=ill+1,ihl
               xc(i,j)   = ( xu(i,j) +  xu(i-1,j) ) / 2
            end do
         end do

         do j=jll+1,jhl
            do i=ill,ihl
               yc(i,j)   = ( yv(i,j) +   yv(i,j-1) ) / 2
            end do
         end do


         if ( have_lonlat ) then

            do j=jll,jhl
               do i=ill,ihl

                  latu(i,j)  = ( latx(i,j) + latx(i,j-1) ) / 2

                  latv(i,j)  = ( latx(i,j) + latx(i-1,j) ) / 2

                  lonc(i,j)  = ( lonx(i-1,j-1) + lonx(i-1,j) &
                               + lonx(i  ,j-1) + lonx(i,j  ) ) / 4

                  latc(i,j)  = ( latx(i-1,j-1) + latx(i-1,j) &
                               + latx(i  ,j-1) + latx(i,j  ) ) / 4

                  convc(i,j) = ( convx(i-1,j-1) + convx(i-1,j) &
                               + convx(i  ,j-1) + convx(i,j  ) ) / 4
               end do
            end do

         end if

      case default
         call getm_error("x2uvc()","A non valid grid type has been chosen.")
   end select

   return
   end subroutine x2uvc
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: metric() - calculate metric coefficients
!
! !INTERFACE:
   subroutine metric()
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Computes the grid increments and areas related to the metric
!  coefficients.
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
! !LOCAL VARIABLES:
   integer                   :: i,j
!EOP
!------------------------------------------------------------------------
!BOC
   rearth = rearth_default

   LEVEL2 'calculating metrics'

   select case (grid_type)

      case(1)  ! Cartesian

         ard1 = _ONE_/(dx*dy)

      case(2)  ! Spherical

!        note that all dy? are identical on constant

         do j=jll,jhl
            dxc(ill:ihl,j)=deg2rad*dlon*rearth*cos(deg2rad*latc(ill:ihl,j))
         end do
         dyc = deg2rad*dlat*rearth

         dxu = dxc
         dyu = dyc

         do j=jll,jhl
            dxv(ill:ihl,j)=deg2rad*dlon*rearth*cos(deg2rad*latx(ill:ihl,j))
         end do
         dyv = dyc

         dxx = dxv
         dyx = dyc

         LEVEL3 'dxc= [ ',minval(dxc,mask=(az .gt. 0)), &
                          maxval(dxc,mask=(az .gt. 0)), ' ]'
         LEVEL3 'dxu= [ ',minval(dxu,mask=(au .gt. 0)), &
                          maxval(dxu,mask=(au .gt. 0)), ' ]'
         LEVEL3 'dxv= [ ',minval(dxv,mask=(av .gt. 0)), &
                          maxval(dxv,mask=(av .gt. 0)), ' ]'
         LEVEL3 'dxx= [ ',minval(dxx,mask=(ax .gt. 0)), &
                          maxval(dxx,mask=(ax .gt. 0)), ' ]'
         LEVEL3 'dy[cuvx]=',minval(dyc,mask=(az .gt. 0))

      case(3) ! planar curvi-linear

         do j=jll+1,jhl
            do i=ill+1,ihl
               dxc(i,j)=sqrt((xu(i,j)-xu(i-1,j))**2+(yu(i,j)-yu(i-1,j))**2)
               dyc(i,j)=sqrt((xv(i,j)-xv(i,j-1))**2+(yv(i,j)-yv(i,j-1))**2)
            end do
         end do

         do j=jll+1,jhl
            do i=ill+1,ihl-1
               dxu(i,j)=sqrt((xc(i+1,j)-xc(i,j))**2+(yc(i+1,j)-yc(i,j))**2)
            end do
         end do

         do j=jll,jhl
            do i=ill,ihl
               dyu(i,j)=sqrt((xx(i,j)-xx(i,j-1))**2+(yx(i,j)-yx(i,j-1))**2)
            end do
         end do

         do j=jll,jhl
            do i=ill,ihl
               dxv(i,j)=sqrt((xx(i,j)-xx(i-1,j))**2+(yx(i,j)-yx(i-1,j))**2)
            end do
         end do

         do j=jll+1,jhl-1
            do i=ill+1,ihl
               dyv(i,j)=sqrt((xc(i,j+1)-xc(i,j))**2+(yc(i,j+1)-yc(i,j))**2)
            end do
         end do

         do j=jll,jhl
            do i=ill,ihl-1
               dxx(i,j)=sqrt((xv(i+1,j)-xv(i,j))**2+(yv(i+1,j)-yv(i,j))**2)
            end do
         end do

         do j=jll,jhl-1
            do i=ill,ihl
               dyx(i,j)=sqrt((xu(i,j+1)-xu(i,j))**2+(yu(i,j+1)-yu(i,j))**2)
            end do
         end do

         LEVEL3 'dxc= [ ',minval(dxc,mask=(az .gt. 0)), &
                          maxval(dxc,mask=(az .gt. 0)), ' ]'
         LEVEL3 'dyc= [ ',minval(dyc,mask=(az .gt. 0)), &
                          maxval(dyc,mask=(az .gt. 0)), ' ]'
         LEVEL3 'dxu= [ ',minval(dxu,mask=(au .gt. 0)), &
                          maxval(dxu,mask=(au .gt. 0)), ' ]'
         LEVEL3 'dyu= [ ',minval(dyu,mask=(au .gt. 0)), &
                          maxval(dyu,mask=(au .gt. 0)), ' ]'
         LEVEL3 'dxv= [ ',minval(dxv,mask=(av .gt. 0)), &
                          maxval(dxv,mask=(av .gt. 0)), ' ]'
         LEVEL3 'dyv= [ ',minval(dyv,mask=(av .gt. 0)), &
                          maxval(dyv,mask=(av .gt. 0)), ' ]'
         LEVEL3 'dxx= [ ',minval(dxx,mask=(ax .gt. 0)), &
                          maxval(dxx,mask=(ax .gt. 0)), ' ]'
         LEVEL3 'dyx= [ ',minval(dyx,mask=(ax .gt. 0)), &
                          maxval(dyx,mask=(ax .gt. 0)), ' ]'

   case(4) ! sperical curvi-linear

      do j=jmin,jmax
         do i=imin,imax
            dx = deg2rad*(lonu(i,j)-lonu(i-1,j))*rearth*cos(deg2rad*latc(i,j))
            dy = deg2rad*(latu(i,j)-latu(i-1,j))*rearth
            dxc(i,j)= sqrt(dx*dx+dy*dy)
            dx = deg2rad*(lonv(i,j)-lonv(i,j-1))*rearth*cos(deg2rad*latc(i,j))
            dy = deg2rad*(latv(i,j)-latv(i,j-1))*rearth
            dyc(i,j)= sqrt(dx*dx+dy*dy)
         end do
      end do

      do j=jmin,jmax
         do i=imin,imax
            dx = deg2rad*(lonc(i+1,j)-lonc(i,j))*rearth*cos(deg2rad*latu(i,j))
            dy = deg2rad*(latc(i+1,j)-latc(i,j))*rearth
            dxu(i,j)= sqrt(dx*dx+dy*dy)
            dx = deg2rad*(lonx(i,j)-lonx(i,j-1))*rearth*cos(deg2rad*latu(i,j))
            dy = deg2rad*(latx(i,j)-latx(i,j-1))*rearth
            dyu(i,j)= sqrt(dx*dx+dy*dy)
         end do
      end do

      do j=jmin,jmax
         do i=imin,imax
            dx = deg2rad*(lonx(i,j)-lonx(i-1,j))*rearth*cos(deg2rad*latv(i,j))
            dy = deg2rad*(latx(i,j)-latx(i-1,j))*rearth
            dxv(i,j)= sqrt(dx*dx+dy*dy)
            dx = deg2rad*(lonc(i,j+1)-lonc(i,j))*rearth*cos(deg2rad*latv(i,j))
            dy = deg2rad*(latc(i,j+1)-latc(i,j))*rearth
            dyv(i,j)= sqrt(dx*dx+dy*dy)
         end do
      end do

      do j=jmin,jmax
         do i=imin,imax
            dx = deg2rad*(lonv(i+1,j)-lonv(i,j))*rearth*cos(deg2rad*latx(i,j))
            dy = deg2rad*(latv(i+1,j)-latv(i,j))*rearth
            dxx(i,j)= sqrt(dx*dx+dy*dy)
            dx = deg2rad*(lonu(i,j+1)-lonu(i,j))*rearth*cos(deg2rad*latx(i,j))
            dy = deg2rad*(latu(i,j+1)-latu(i,j))*rearth
            dyx(i,j)= sqrt(dx*dx+dy*dy)
         end do
      end do

      case default
         call getm_error("metric()","A non valid grid type has been chosen.")
   end select

   if ( grid_type .ne. 1 ) then

      call update_2d_halo(dxc,dxc,az,imin,jmin,imax,jmax,H_TAG)
      call wait_halo(H_TAG)

      call update_2d_halo(dyc,dyc,az,imin,jmin,imax,jmax,H_TAG)
      call wait_halo(H_TAG)

      call update_2d_halo(dxu,dxu,au,imin,jmin,imax,jmax,U_TAG)
      call wait_halo(U_TAG)

      call update_2d_halo(dyu,dyu,au,imin,jmin,imax,jmax,U_TAG)
      call wait_halo(U_TAG)

      call update_2d_halo(dxv,dxv,av,imin,jmin,imax,jmax,V_TAG)
      call wait_halo(V_TAG)

      call update_2d_halo(dyv,dyv,av,imin,jmin,imax,jmax,V_TAG)
      call wait_halo(V_TAG)

      call update_2d_halo(dxx,dxx,ax,imin,jmin,imax,jmax,H_TAG)
      call wait_halo(H_TAG)

      call update_2d_halo(dyx,dyx,ax,imin,jmin,imax,jmax,H_TAG)
      call wait_halo(H_TAG)


!     compute differently centered areas of grid boxes
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO

            if( az(i,j) .gt. 0) then
               arcd1(i,j)=_ONE_/(dxc(i,j)*dyc(i,j))
            end if

            if( au(i,j) .gt. 0) then
               arud1(i,j)=_ONE_/(dxu(i,j)*dyu(i,j))
            end if

            if( av(i,j) .gt. 0) then
               arvd1(i,j)=_ONE_/(dxv(i,j)*dyv(i,j))
            end if

         end do
      end do

      call update_2d_halo(arcd1,arcd1,az,imin,jmin,imax,jmax,H_TAG)
      call wait_halo(H_TAG)

      call update_2d_halo(arud1,arud1,au,imin,jmin,imax,jmax,U_TAG)
      call wait_halo(U_TAG)

      call update_2d_halo(arvd1,arvd1,av,imin,jmin,imax,jmax,V_TAG)
      call wait_halo(V_TAG)

   end if

   return
   end subroutine metric
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_min_depth() - set the minimum depth in regions
!
! !INTERFACE:
   subroutine set_min_depth(fn)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Read region definitions and minimum depth for those regions. Adjust the
!  bathymetry (variable $H$) accordingly.
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
   integer                   :: unit = 25 ! kbk
   character(len=255)        :: line
   integer                   :: iostat
   integer                   :: i,j,k=0,n=-1
   integer                   :: il,jl,ih,jh
   integer                   :: i1,j1
   REALTYPE                  :: dmin
!EOP
!-----------------------------------------------------------------------
!BOC
!   open(unit,file=fn,action='read',iostat=iostat,status='old',err=90)
   open(unit,file=fn,action='read',iostat=iostat,status='old')

   do while (iostat == 0)
      read(unit,'(A)',iostat=iostat,end=91,err=92) line
!     skip comments and empty lines
      if (line(1:1) == '#' .or. line(1:1) == '!' .or. len(trim(line)) == 0 ) then
      else if ( n .eq. -1 ) then
         read(line,*) n
         if(n .gt. 1) then
            LEVEL2 'setting minimum depths according to:'
            LEVEL3 trim(fn)
         end if
      else
         read(line,*,iostat=iostat) il,jl,ih,jh,dmin
         if (iostat .ne. 0) goto 93
         k = k+1
         LEVEL3 il,jl,ih,jh,dmin
         do j=jl,jh
            do i=il,ih
               if(imin+ioff .le. i .and. i .le. imax+ioff .and. &
                  jmin+joff .le. j .and. j .le. jmax+joff ) then
                  i1 = i-ioff
                  j1 = j-joff
                  if(H(i1,j1) .gt. -9. .and. H(i1,j1) .lt. dmin) then
                     H(i1,j1) = dmin
                  end if
               end if
            end do
         end do
      end if
   end do

   close(unit)
   return

90 LEVEL2 'could not open ',trim(fn),' no minimum depth adjustments done'
91 LEVEL2 'done'
   return
92 call getm_error("set_min_depth()","End of file "//trim(fn)//".")
93 call getm_error("set_min_depth()","Error reading line: "//trim(line))
   end subroutine set_min_depth
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: adjust_bathymetry() - read mask adjustments from file.
!
! !INTERFACE:
   subroutine adjust_bathymetry(fn)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Read bathymetry adjustments from file.
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
   integer                   :: unit = 25 ! kbk
   character(len=255)        :: line
   integer                   :: iostat
   integer                   :: i,j,k=0,n=-1
   integer                   :: il,jl,ih,jh
   REALTYPE                  :: x
!EOP
!-----------------------------------------------------------------------
!BOC
!   open(unit,file=fn,action='read',iostat=iostat,status='old',err=90)
   open(unit,file=fn,action='read',iostat=iostat,status='old')

   do while (iostat == 0)
      read(unit,'(A)',iostat=iostat,end=91,err=92) line
!     skip comments and empty lines
      if (line(1:1) == '#' .or. line(1:1) == '!' .or. len(trim(line)) == 0 ) then
      else if ( n .eq. -1 ) then
         read(line,*) n
         if(n .gt. 1) then
            LEVEL2 'adjusting bathymetry according to:'
            LEVEL3 trim(fn)
         end if
      else
         read(line,*,iostat=iostat) il,jl,ih,jh,x
         if (iostat .ne. 0) goto 93
         k = k+1
         LEVEL3 il,jl,ih,jh,x
         do j=jl,jh
            do i=il,ih
               if(imin+ioff-HALO .le. i .and. i .le. imax+ioff+HALO .and. &
                  jmin+joff-HALO .le. j .and. j .le. jmax+joff+HALO ) then
                  H(i-ioff,j-joff) = x
               end if
            end do
         end do
      end if
   end do

   close(unit)
   return

90 LEVEL2 'could not open ',trim(fn),' no bathymetry adjustments done'
91 LEVEL2 'done'
   return
92 call getm_error("adjust_bathymetry()","End of file "//trim(fn)//".")
93 call getm_error("adjust_bathymetry()","Error reading line: "//trim(line))
   end subroutine adjust_bathymetry
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: adjust_mask() - read mask adjustments from file.
!
! !INTERFACE:
   subroutine adjust_mask(fn)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Read mask adjustments from file. The file format allows comments.
!  Comment characters are ! or \# - they MUST be in column 1.
!  Lines with white-spaces are skipped. Conversion errors
!  are caught and an error condition occurs.
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
   integer                   :: unit = 25 ! kbk
   character(len=255)        :: line
   integer                   :: iostat
   integer                   :: i,j,k=0,n=-1
   integer                   :: il,jl,ih,jh
!EOP
!-----------------------------------------------------------------------
!BOC
!   open(unit,file=fn,action='read',iostat=iostat,status='old',err=90)
   open(unit,file=fn,action='read',iostat=iostat,status='old')

   do while (iostat == 0)
      read(unit,'(A)',iostat=iostat,end=91,err=92) line
!     skip comments and empty lines
      if (line(1:1) == '#' .or. line(1:1) == '!' .or. len(trim(line)) == 0 ) then
      else if ( n .eq. -1 ) then
         read(line,*) n
         if(n .ge. 1) then
            LEVEL2 'adjusting mask according to:'
            LEVEL3 trim(fn)
         end if
      else
         read(line,*,iostat=iostat) il,jl,ih,jh
         if (iostat .ne. 0) goto 93
         k = k+1
         LEVEL3 il,jl,ih,jh
         do j=jl,jh
            do i=il,ih
               if(imin+ioff-HALO .le. i .and. i .le. imax+ioff+HALO .and. &
                  jmin+joff-HALO .le. j .and. j .le. jmax+joff+HALO ) then
                  az(i-ioff,j-joff) = 0
               end if
            end do
         end do
      end if
   end do

   close(unit)
   return

90 LEVEL2 'could not open ',trim(fn),' no mask adjustments done'
91 LEVEL2 'done'
   return
92 call getm_error("adjust_mask()","Error reading "//trim(fn))
93 call getm_error("adjust_mask()","Error reading line: "//trim(line))
   end subroutine adjust_mask
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: print_mask() - prints a mask in readable format
!
! !INTERFACE:
   subroutine print_mask(mask)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Prints a integer mask in a human readable form.
!
! !INPUT PARAMETERS:
   integer, intent(in), dimension(E2DFIELD) :: mask
!
! !REVISION HISTORY:
!
!  22Apr99   Karsten Bolding & Hans Burchard  Initial code.
!
! !LOCAL VARIABLES:
   integer                   :: i,j
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
#endif

#if 0
   do j=jmax+HALO,jmin-HALO,-1
!      write(0,'(5000(i1,1x))') (mask(i,j), i=imin,imax)
      write(0,'(5000(i1))') (mask(i,j), i=imin-HALO,imax+HALO,1)
   end do
#else
   do j=jmax,jmin,-1
!      write(0,'(5000(i1,1x))') (mask(i,j), i=imin,imax)
      write(0,'(5000(i1))') (mask(i,j), i=imin,imax,1)
   end do
#endif

   return
   end subroutine print_mask
!EOC

!-----------------------------------------------------------------------

   end module domain

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
