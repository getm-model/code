!$Id: domain.F90,v 1.33 2009-08-31 10:37:03 bjb Exp $
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
!
! !USES:
   use exceptions
   use topo_interface, only: check_grid,get_grid
   use halo_zones,     only: update_2d_halo,wait_halo
   use halo_zones,     only: H_TAG,U_TAG,V_TAG

   IMPLICIT NONE
!
!
! !PUBLIC DATA MEMBERS:
   integer                             :: bathy_format   = NETCDF

   integer                             :: grid_type      = 1
   integer                             :: vert_cord      = 1
   integer                             :: proj_type      = 1
!  global index range
   integer                             :: ilg=-1,ihg=-1,jlg=-1,jhg=-1
!  local index range
   integer                             :: ill=-1,ihl=-1,jll=-1,jhl=-1

   logical                             :: latlon_exists  = .false.
   logical                             :: xy_exists      = .false.
   logical                             :: updateXYC      = .true.
   logical                             :: updateLatLonC  = .true.
   logical                             :: updateConvC    = .true.

   logical                             :: proj_exists    = .true.
   REALTYPE                            :: proj_lon,proj_lat,proj_rot
   REALTYPE                            :: rearth

   REALTYPE                            :: maxdepth       = -1.
   REALTYPE                            :: ddu            = -_ONE_
   REALTYPE                            :: ddl            = -_ONE_
   REALTYPE                            :: d_gamma        = 20.
   logical                             :: gamma_surf     = .true.
   REALTYPE, allocatable, dimension(:) :: ga

   integer                             :: NWB,NNB,NEB,NSB,NOB
   integer                             :: calc_points
   logical                             :: openbdy        = .false.

   REALTYPE                            :: Hland
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

   character(len=64)                   :: bdy_2d_desc(5)
   logical                             :: need_2d_bdy_elev = .false.
   logical                             :: need_2d_bdy_u    = .false.
   logical                             :: need_2d_bdy_v    = .false.

   REALTYPE                            :: cori= _ZERO_

!  method for specifying bottom roughness (0=const, 1=from topo.nc)
   integer                             :: z0_method=0
   REALTYPE                            :: z0_const=0.001

! !DEFINED PARAMETERS:
   integer,           parameter        :: INNER          = 1
   REALTYPE, private, parameter        :: pi             = 3.141592654
   REALTYPE, private, parameter        :: deg2rad        = pi/180.
   REALTYPE, private, parameter        :: omega          = 2.*pi/86400.
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: domain.F90,v $
!  Revision 1.33  2009-08-31 10:37:03  bjb
!  Consistent treatment of topo in halo zones
!
!  Revision 1.32  2009-05-15 06:59:10  bjb
!  typo fix
!
!  Revision 1.31  2009-05-07 15:50:46  kb
!  added global and local horizontal index range
!
!  Revision 1.30  2009-05-07 10:10:15  kb
!  fixed tag in wait_halo() - Buchmann
!
!  Revision 1.29  2008-12-09 00:31:57  kb
!  added new 2D open boundaries
!
!  Revision 1.28  2007-10-16 06:22:56  kbk
!  curvi-linear now runs in parallel
!
!  Revision 1.27  2007-03-30 13:10:59  hb
!  Use of adaptive and hybrid vertical coordinates technically enabled
!
!  Revision 1.26  2007-02-08 06:43:27  kbk
!  update HALOS for z0 - and changed loop boundaries for zub0,zvb0
!
!  Revision 1.25  2007-02-07 16:32:22  kbk
!  added spatial varying bottom roughness
!
!  Revision 1.24  2007-02-07 16:27:06  kbk
!  changed and fixed loop boundaries
!
!  Revision 1.23  2006-08-25 09:05:20  kbk
!  metric coefficients also calculated in HALO zones
!
!  Revision 1.22  2006-06-03 11:43:16  kbk
!  added namelist fallback longitude - heatfluxes
!
!  Revision 1.21  2005-06-27 07:18:13  frv-bjb
!  Changed STOP statements to call getm_error(...)
!
!  Revision 1.20  2005/06/17 07:40:19  frv-bjb
!  Added check/bailout for zero dlat and dlon
!
!  Revision 1.19  2005/05/25 10:43:42  kbk
!  fixed ax calculation
!
!  Revision 1.18  2005/05/25 10:32:12  kbk
!  merged from stabe branch v1_2_1
!
!  Revision 1.17  2005/04/25 09:32:34  kbk
!  added NetCDF IO rewrite + de-stag of velocities - Umlauf
!
!  Revision 1.16  2004/11/04 11:07:00  kbk
!  fixed format statement in print_mask
!
!  Revision 1.15  2004/01/05 13:24:27  kbk
!  maxdepth from domain namelist - should be calculated later
!
!  Revision 1.14  2003/09/02 14:12:14  kbk
!  au and av also in HALO-zones
!
!  Revision 1.13  2003/08/28 10:36:30  kbk
!  also calculate ax in HALO-zones
!
!  Revision 1.12  2003/08/21 15:28:29  kbk
!  re-enabled update_2d_halo for lonc and latc + cleaning
!
!  Revision 1.11  2003/08/15 12:52:49  kbk
!  moved az mask calculation + removed print statements
!
!  Revision 1.10  2003/08/03 09:52:11  kbk
!  nicer print statements
!
!  Revision 1.9  2003/06/29 17:09:04  kbk
!  removed reference to barrier
!
!  Revision 1.8  2003/05/09 11:52:08  kbk
!  do not mirror coordinate info + use mask for inverse area calculation
!
!  Revision 1.7  2003/05/02 08:32:31  kbk
!  re-ordering mask calculation
!
!  Revision 1.6  2003/04/23 11:59:39  kbk
!  update_2d_halo on spherical variables + TABS to spaces
!
!  Revision 1.5  2003/04/07 14:34:42  kbk
!  parallel support, proper spherical grid init. support
!
!  Revision 1.1.1.1  2002/05/02 14:01:11  gotm
!  recovering after CVS crash
!
!  Revision 1.19  2001/10/23 14:15:55  bbh
!  Moved ga from coordinates.F90 to domain.F90
!
!  Revision 1.18  2001/10/22 12:10:26  bbh
!  Partly support for SPHERICAL grid is coded
!
!  Revision 1.17  2001/09/26 10:01:41  bbh
!  lat and lon maps now read in ncdf_topo.F90
!
!  Revision 1.16  2001/09/24 07:49:32  bbh
!  Include .h files for memory declaration/allocation
!
!  Revision 1.15  2001/09/21 11:52:47  bbh
!  Minimum depth for specific areas through - set_min_depth()
!
!  Revision 1.14  2001/09/14 12:04:15  bbh
!  Added xc,yc to hold coordinates + cleaning
!
!  Revision 1.13  2001/09/04 08:00:14  bbh
!  Fill coru and corv arrays
!
!  Revision 1.12  2001/09/04 07:36:32  bbh
!  We need ioff and joff in parallel runs
!
!  Revision 1.11  2001/09/03 15:14:22  bbh
!  Bug with Coriolis removed
!
!  Revision 1.10  2001/09/01 17:10:25  bbh
!  Vertical coordinate definition now specified via namelist
!
!  Revision 1.9  2001/08/27 11:55:02  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.8  2001/08/01 08:19:57  bbh
!  Fields for CURVILINEAR - now done
!
!  Revision 1.7  2001/07/26 14:31:43  bbh
!  Manual merge
!
!  Revision 1.6  2001/07/26 14:20:02  bbh
!  Added grid_type, vert_cord, lonmap and latmap
!
!  Revision 1.5  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.4  2001/05/14 12:38:58  bbh
!  Set minimum detph to 10. meters if not COAST_TEST - to be fixed later.
!
!  Revision 1.3  2001/05/06 18:51:55  bbh
!  Towards proper implementation of specified 2D bdy.
!
!  Revision 1.2  2001/04/24 08:24:58  bbh
!  Use runtype instead of macro
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:

!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_domain - initialise the computational domain

!
! !INTERFACE:
   subroutine init_domain(input_dir)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
   character(len=*)                    :: input_dir
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Initialize the calculation domain - that is obtain information of the
!  extent of the domain - imin,imax,jmin and jmax.
!  This subroutine operates in two principally different ways.
!  In the first case it has been compiled with \#define STATIC which means that
!  most arrays are statically allocated. In this case the information of
!  the extenstion of the domain is obtained from a file called dimensions.h,
!  where imin,imax,jmin and jmax are specified as integer parameters.
!  If STATIC has not been defined during compilation all arrays are allocatable.
!  In this case information of the extension of the domain is obtained from
!  the bathymetry file. After imin,imax,jmin and jmax are known arrays are
!  allocated.
!  Next step is to read in the batymetry. Then location of boundaries are read.
!  Hereafter calculation masks are setup - note that the masks are explicitely
!  set to 0 in the halo zones.
!  Finally the halo zones are filled. In the case of only one process the
!  information is simply copied where as in a parallel run the information is
!  communicated from the neighboring processes - this is all done in
!  update\_2d\_halo().
!
! !REVISION HISTORY:
!
!  See log for module
!
! !LOCAL VARIABLES:
   integer                   :: rc
   integer                   :: np,sz
   integer                   :: i,j,n
   integer                   :: kdum
   character(len=PATH_MAX)   :: bathymetry               = 'topo.nc'
   character(len=PATH_MAX)   :: bdyinfofile              = 'bdyinfo.dat'
   character(len=PATH_MAX)   :: min_depth_file           = 'minimum_depth.dat'
   character(len=PATH_MAX)   :: bathymetry_adjust_file   = 'bathymetry.adjust'
   character(len=PATH_MAX)   :: mask_adjust_file         = 'mask.adjust'
   integer                   :: il=-1,ih=-1,jl=-1,jh=-1
   REALTYPE                  :: mask(E2DFIELD)
   namelist /domain/ &
             vert_cord,maxdepth,bathy_format,bathymetry,       &
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
   bdy_2d_desc(SOMMERFELDT)             = "Sommerfeldt rad."
   bdy_2d_desc(CLAMPED)                 = "Clamped"
   bdy_2d_desc(FLATHER_ELEV)            = "Flather (elev)"
   bdy_2d_desc(FLATHER_VEL)             = "Flather (vel)"

   LEVEL1 'init_domain'

!  Read domain specific things from the namelist.
   read(NAMLST,domain)

!  check grid and dimensions
   call check_grid(bathymetry,bathy_format,iextr,jextr)

   select case (vert_cord)
      case(1)
         LEVEL2 'Using sigma coordinates'
      case(2)
         LEVEL2 'Using z-level coordinates'
      case(3)
         LEVEL2 'Using general vertical coordinates'
      case (4) ! hybrid vertical coordinates
         LEVEL2 'using hybrid vertical coordinates'
         STDERR 'domain: hybrid_coordinates not coded yet'
         stop
      case (5) ! adaptive vertical coordinates
         LEVEL2 'using adaptive vertical coordinates'
         STDERR 'domain: adaptive_coordinates not coded yet'
         stop
      case default
         call getm_error("init_domain()", &
                         "A non valid vertical coordinate system has been chosen");
   end select

   select case (z0_method)
      case(0)
         LEVEL2 'Using constant bottom roughness'
      case(1)
         LEVEL2 'Using space varying bottom roughness'
         LEVEL2 '..  will read z0 from the topo file ..'
      case default
         call getm_error("init_domain()", &
                         "A non valid z0 method has been chosen");
   end select

#ifndef STATIC
   kmax=kdum
#endif
! prepare parallel run
   call part_domain()
   il=imin ; ih=imax ; jl=jmin ; jh=jmax
#ifndef STATIC
#include "dynamic_allocations_domain.h"
#endif

!  GLOBAL horizontal index range for 2D and 3D variables
   ilg = max(imin+ioff,1); ihg = min(imax+ioff,iextr)
   jlg = max(jmin+joff,1); jhg = min(jmax+joff,jextr)

!  LOCAL horizontal index range for 2D and 3D variables
   ill = max(imin-ioff,1); ihl = min(imax,iextr-ioff)
   jll = max(jmin-joff,1); jhl = min(jmax,jextr-joff)

   H = -10.
   HU = -10.
   HV = -10.

   lonc = -1000.
   latc = -1000.

!  read grid and bathymetry
   call get_grid(bathy_format,H,Hland,iextr,jextr,ioff,joff, &
                imin,imax,jmin,jmax)

!  prepare parallel setup
   call update_2d_halo(lonc,lonc,az,imin,jmin,imax,jmax,H_TAG, &
                       mirror=.false.)
   call wait_halo(H_TAG)
   call update_2d_halo(latc,latc,az,imin,jmin,imax,jmax,H_TAG, &
                       mirror=.false.)
   call wait_halo(H_TAG)

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
      if(wfj(n) .eq. jmin) az(wi(n),jmin-1) = az(wi(n),jmin)
      if(wlj(n) .eq. jmax) az(wi(n),jmax+1) = az(wi(n),jmax)
   end do
!  northern boundary - at present elev only
   do n=1,NNB
      az(nfi(n):nli(n),nj(n)) = BOUNDARY_POINT
      if(nfi(n) .eq. imin) az(imin-1,nj(n)) = az(imin,nj(n))
      if(nli(n) .eq. imax) az(imax+1,nj(n)) = az(imax,nj(n))
   end do
!  easter boundary - at present elev only
   do n=1,NEB
      az(ei(n),efj(n):elj(n)) = BOUNDARY_POINT
      if(efj(n) .eq. jmin) az(ei(n),jmin-1) = az(ei(n),jmin)
      if(elj(n) .eq. jmax) az(ei(n),jmax+1) = az(ei(n),jmax)
   end do
!  southern boundary - at present elev only
   do n=1,NSB
      az(sfi(n):sli(n),sj(n)) = BOUNDARY_POINT
      if(sfi(n) .eq. imin) az(imin-1,sj(n)) = az(imin,sj(n))
      if(sli(n) .eq. imax) az(imax+1,sj(n)) = az(imax,sj(n))
   end do
#undef BOUNDARY_POINT

!  Do we want to further adjust the mask
   call adjust_mask(trim(input_dir) // mask_adjust_file)

   mask = _ONE_*az

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

!  Compute grid points and metric coefficients for different grid types
   select case (grid_type)

!  CARTESIAN
   case(1)

!     Generate xx and yx from grid spacing
      do j = jmin-1,jmax
         xx(imin-1,j) = ioff*dx + x0
         do i=imin,imax
            xx(i,j) = xx(i-1,j) + dx
         end do
      end do

      do i=imin-1,imax
         yx(i,jmin-1) = joff*dy + y0
         do j=jmin,jmax
            yx(i,j) = yx(i,j-1) + dy
         end do
      end do

!     Exchange halo information for parallel runs
      call update_2d_halo(xx,xx,ax,imin,jmin,imax,jmax,H_TAG, mirror=.false.)
      call wait_halo(H_TAG)

      call update_2d_halo(yx,yx,ax,imin,jmin,imax,jmax,H_TAG,mirror=.false.)
      call wait_halo(H_TAG)

!     xx and yx points exist now
      xy_exists = .true.

!     Interpolate (latx,lonx) and (xx,yx) to the u, v, and T-points.
      call x2uvc(updateXYC,updateLatLonC,updateConvc)

!     Compute metric coefficients
      ard1 = _ONE_/(dx*dy)

!  SPHERICAL
   case(2)

!     Generate lonx,latx from grid spacing and offset
!
!     First make test that we have good values for dlat and dlon
!     (if they are zer the program will fail later due to dx=0 etc)
      if (dlat .eq. 0.0) then
         call getm_error("init_domain()", &
                         "Delta lat (dlat) seems to be zero!")
      end if
      if (dlon .eq. 0.0) then
         call getm_error("init_domain()", &
                         "Delta lon (dlon) seems to be zero!")
      end if
!
!     Then actually generate lonx,latx
      do j = jmin-1,jmax+1
         lonx(imin-1,j) = ioff*dlon + lon0
         do i=imin,imax+1
            lonx(i,j) = lonx(i-1,j) + dlon
         end do
      end do

      do i=imin-1,imax+1
         latx(i,jmin-1) = joff*dlat + lat0
         do j=jmin,jmax+1
            latx(i,j) = latx(i,j-1) + dlat
         end do
      end do

!     Exchange halo information for parallel runs
      call update_2d_halo(lonx,lonx,ax,imin,jmin,imax,jmax,H_TAG,mirror=.false.)
      call wait_halo(H_TAG)

      call update_2d_halo(latx,latx,ax,imin,jmin,imax,jmax,H_TAG,mirror=.false.)
      call wait_halo(H_TAG)

!     lat and long exist now
      latlon_exists = .true.

!     Interpolate (latx,lonx) and (xx,yx) to the u, v, and T-points.
      call x2uvc(updateXYC,updateLatLonC,updateConvc)

!     Compute metric coefficients
      call metric(grid_type)

!  PLANE/SPHERICAL CURVILINEAR
   case(3,4)

!     Interpolate (latx,lonx) and (xx,yx) to the u, v, and T-points.
      call x2uvc(updateXYC,updateLatLonC,updateConvc)

!     Compute metric coefficients
      call metric(grid_type)

   case default
      call getm_error("init_domain()","A non valid grid type has been chosen.")
   end select


   if ( .not. latlon_exists ) then
      LEVEL2 "Setting constant longitude (swr)"
      lonc = longitude
   end if

!  Compute Coriolis parameter
   if (f_plane) then
      LEVEL2 "Assuming constant Coriolis parameter."
      cori = 2.*omega*sin(deg2rad*latitude)
      coru = cori
      corv = cori
   else
      if (latlon_exists) then
         LEVEL2 "Computing spatially varying Coriolis parameter from (lat,lon)."

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

   if (z0_method .eq. 1) then
      call update_2d_halo(z0,z0,az,imin,jmin,imax,jmax,H_TAG)
      call wait_halo(H_TAG)
   end if

#ifdef DEBUG
   STDERR 'az'
   call print_mask(az)
   STDERR 'au'
   call print_mask(au)
   STDERR 'av'
   call print_mask(av)
#endif

   np = count(az(1:imax,1:jmax) .gt. 0)
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
   subroutine x2uvc(updateXYC,updateLatLonC,updateConvc)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
  logical, intent(in)        :: updateXYC
  logical, intent(in)        :: updateLatLonC
  logical, intent(in)        :: updateConvC

!
!
! !DESCRIPTION:
!  This routine interpolates (latx,lonx), (xx,yx), and convx to the
!  u-points, v-points, and the central T-points. The data at the T-points
!  are only updated from values of the X-points if the logical flags {\tt updateXYC},
!  {\tt updateXYC}, and {\tt updateXYC} are {\tt .true.}. This is not necessary
!  if data at the T-points have been read from the topo input file.
!
! !REVISION HISTORY:
! Original author(s): Lars Umlauf
!
! !LOCAL VARIABLES:
   integer                   :: i,j
!
!EOP
!------------------------------------------------------------------------
!BOC
   do i=imin-HALO,imax+HALO
      do j=jmin-HALO+1,jmax+HALO
         xu(i,j)   = 0.5*(  xx(i,j) +   xx(i,j-1))
         yu(i,j)   = 0.5*(  yx(i,j) +   yx(i,j-1))
         latu(i,j) = 0.5*(latx(i,j) + latx(i,j-1))
         lonu(i,j) = 0.5*(lonx(i,j) + lonx(i,j-1))
      end do
   end do

   call update_2d_halo(xu,xu,au,imin,jmin,imax,jmax,U_TAG,mirror=.false.)
   call wait_halo(U_TAG)
   call update_2d_halo(yu,yu,au,imin,jmin,imax,jmax,U_TAG,mirror=.false.)
   call wait_halo(U_TAG)

   call update_2d_halo(latu,latu,au,imin,jmin,imax,jmax,U_TAG,mirror=.false.)
   call wait_halo(U_TAG)
   call update_2d_halo(lonu,lonu,au,imin,jmin,imax,jmax,U_TAG,mirror=.false.)
   call wait_halo(U_TAG)


   do i=imin-HALO+1,imax+HALO
      do j=jmin-HALO,jmax+HALO
         xv(i,j)   = 0.5*(  xx(i,j) +   xx(i-1,j))
         yv(i,j)   = 0.5*(  yx(i,j) +   yx(i-1,j))
         latv(i,j) = 0.5*(latx(i,j) + latx(i-1,j))
         lonv(i,j) = 0.5*(lonx(i,j) + lonx(i-1,j))
      end do
   end do

   call update_2d_halo(xv,xv,av,imin,jmin,imax,jmax,V_TAG,mirror=.false.)
   call wait_halo(V_TAG)
   call update_2d_halo(yv,yv,av,imin,jmin,imax,jmax,V_TAG,mirror=.false.)
   call wait_halo(V_TAG)

   call update_2d_halo(latv,latv,av,imin,jmin,imax,jmax,V_TAG,mirror=.false.)
   call wait_halo(V_TAG)
   call update_2d_halo(lonv,lonv,av,imin,jmin,imax,jmax,V_TAG,mirror=.false.)
   call wait_halo(V_TAG)


   if (updateXYC) then
      do i=imin,imax
         do j=jmin,jmax
            xc(i,j)   = 0.5*(  xu(i,j) +  xu(i-1,j))
         end do
      end do

      do i=imin,imax
         do j=jmin,jmax
            yc(i,j)   = 0.5*(  yv(i,j) +   yv(i,j-1))
         end do
      end do
   endif

   if (updateLatLonC) then
      do i=imin,imax
         do j=jmin,jmax
            latc(i,j) = 0.5*(latu(i,j) +latu(i-1,j))
         end do
      end do

      do i=imin,imax
         do j=jmin,jmax
            lonc(i,j) = 0.5*(lonv(i,j) + lonv(i,j-1))
         end do
      end do
   endif

   if (updateConvc) then
      do i=imin,imax
         do j=jmin,jmax
            convc(i,j)   = 0.25*(  convx(i-1,j-1) + convx(i-1,j)    &
                                 + convx(i  ,j-1) + convx(i,j  ) )
         end do
      end do
   endif


   call update_2d_halo(xc,xc,az,imin,jmin,imax,jmax,H_TAG,mirror=.false.)
   call wait_halo(H_TAG)

   call update_2d_halo(yc,yc,az,imin,jmin,imax,jmax,H_TAG,mirror=.false.)
   call wait_halo(H_TAG)

   call update_2d_halo(latc,latc,az,imin,jmin,imax,jmax,H_TAG,mirror=.false.)
   call wait_halo(H_TAG)

   call update_2d_halo(lonc,lonc,az,imin,jmin,imax,jmax,H_TAG,mirror=.false.)
   call wait_halo(H_TAG)

   call update_2d_halo(convc,convc,az,imin,jmin,imax,jmax,H_TAG,mirror=.false.)
   call wait_halo(H_TAG)


   return
   end subroutine x2uvc
!EOC



!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: compute the metric coefficients
!
! !INTERFACE:
   subroutine metric(grid_type)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: grid_type
!
! !DESCRIPTION:
! Computes the grid increments and areas related to the metric coefficients
! of non-Cartesian grids.
!
! !REVISION HISTORY:
! Original author(s): Lars Umlauf
!
! !LOCAL VARIABLES:
   integer                   :: i,j
!
!EOP
!------------------------------------------------------------------------
!BOC

   select case (grid_type)

   case(1)

   case(2)

      do j=jmin-HALO+1,jmax+HALO
         do i=imin-HALO+1,imax+HALO
            dxc(i,j)=deg2rad*(lonu(i,j)-lonu(i-1,j))*rearth &
                 *cos(deg2rad*latc(i,j))
            dyc(i,j)=deg2rad*(latv(i,j)-latv(i,j-1))*rearth
         end do
      end do

      do j=jmin-HALO+1,jmax+HALO
         do i=imin-HALO,imax+HALO-1
            dxu(i,j)=deg2rad*(lonc(i+1,j)-lonc(i,j))*rearth &
                 *cos(deg2rad*latc(i,j))
            dyu(i,j)=deg2rad*(latx(i,j)-latx(i,j-1))*rearth
         end do
      end do

      do j=jmin-HALO,jmax+HALO-1
         do i=imin-HALO+1,imax+HALO
            dxv(i,j)=deg2rad*(lonx(i,j)-lonx(i-1,j))*rearth &
                 *cos(deg2rad*latv(i,j))
            dyv(i,j)=deg2rad*(latc(i,j+1)-latc(i,j))*rearth
         end do
      end do

      do j=jmin-HALO,jmax+HALO-1
         do i=imin-HALO,imax+HALO-1
            dxx(i,j)=deg2rad*(lonv(i+1,j)-lonv(i,j))*rearth &
                 *cos(deg2rad*latx(i,j))
            dyx(i,j)=deg2rad*(latu(i,j+1)-latu(i,j))*rearth
         end do
      end do

   case(3)

      do i=imin,imax
         do j=jmin,jmax
            dxc(i,j)=sqrt((xu(i,j)-xu(i-1,j))**2+(yu(i,j)-yu(i-1,j))**2)
            dyc(i,j)=sqrt((xv(i,j)-xv(i,j-1))**2+(yv(i,j)-yv(i,j-1))**2)
         end do
      end do

      do j=jmin,jmax
         do i=imin,imax
            dxu(i,j)=sqrt((xc(i+1,j)-xc(i,j))**2+(yc(i+1,j)-yc(i,j))**2)
            dyu(i,j)=sqrt((xx(i,j)-xx(i,j-1))**2+(yx(i,j)-yx(i,j-1))**2)
         end do
      end do

      do j=jmin,jmax
         do i=imin,imax
            dxv(i,j)=sqrt((xx(i,j)-xx(i-1,j))**2+(yx(i,j)-yx(i-1,j))**2)
            dyv(i,j)=sqrt((xc(i,j+1)-xc(i,j))**2+(yc(i,j+1)-yc(i,j))**2)
         end do
      end do

      do j=jmin,jmax
         do i=imin,imax
            dxx(i,j)=sqrt((xv(i+1,j)-xv(i,j))**2+(yv(i+1,j)-yv(i,j))**2)
            dyx(i,j)=sqrt((xu(i,j+1)-xu(i,j))**2+(yu(i,j+1)-yu(i,j))**2)
         end do
      end do

   case(4)

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


!  compute differently centered areas of grid boxes
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

   return
   end subroutine metric
!EOC



!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_par_setup - reads domain partition
!
! !INTERFACE:
   subroutine read_par_setup(myid)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: myid
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Reads the partitioning of the domain in parallel run
!
! !REVISION HISTORY:
!
!  22Apr99   Karsten Bolding & Hans Burchard  Initial code.
!
! !LOCAL VARIABLES:
   integer                   :: id
!
!EOP
!------------------------------------------------------------------------
!BOC

#ifndef STATIC
!   open(PARSETUP,file=input_dir // 'par_setup')
   read(PARSETUP,*)

100  read(PARSETUP,*,ERR=110,END=111) id,imin,imax,jmin,jmax


   if(id .eq. myid ) then
      close(PARSETUP)
      ioff=imin-1 ; joff=jmin-1
      imax=imax-imin+1 ; imin=1
      jmax=jmax-jmin+1 ; jmin=1
      LEVEL2 'From read_par_setup ',id,ioff,imin,imax,joff,jmin,jmax
      return
   end if

   goto 100

110 call getm_error("read_par_setup()","reading domain partition information.")

111 call getm_error("read_par_setup()","End of file reached.")
#endif
   return
   end subroutine read_par_setup
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_min_depth -
!
! !INTERFACE:
   subroutine set_min_depth(fn)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Read mask adjustments from file.
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
   integer                   :: unit = 25 ! kbk
   integer                   :: i,j,k,n
   integer                   :: il,jl,ih,jh
   integer                   :: i1,j1
   REALTYPE                  :: dmin
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Should read in to a buffer at some time to allow for #
   open(unit,file=fn,action='read',status='old',err=90)
   read(unit,*,end=91,err=92) n
   if(n .ge. 1) then
      LEVEL2 'setting minimum depths according to:'
      LEVEL3 trim(fn)
   end if
   do k=1,n
      read(unit,*,end=91,err=92) il,jl,ih,jh,dmin
      LEVEL3 'setting min depth in ',il,jl,ih,jh,' to ',dmin
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
   end do
   close(unit)

   return

90 LEVEL2 'could not open ',trim(fn),' no minimum depth adjustments done'
   return
91 call getm_error("set_min_depth()","End of file "//trim(fn)//".")
92 call getm_error("set_min_depth()","Error reading "//trim(fn)//".")
   end subroutine set_min_depth
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: adjust_bathymetry - read mask adjustments from file.
!
! !INTERFACE:
   subroutine adjust_bathymetry(fn)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Read bathymetry adjustments from file.
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
   integer                   :: unit = 25 ! kbk
   integer                   :: i,j,k,n
   integer                   :: il,jl,ih,jh
   integer                   :: i1,j1
   REALTYPE                  :: x
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Should read in to a buffer at some time to allow for #
   open(unit,file=fn,action='read',status='old',err=90)
   read(unit,*,end=91,err=92) n
   if(n .gt. 1) then
      LEVEL2 'adjusting bathymetry according to:'
      LEVEL3 trim(fn)
   end if
   do k=1,n
      read(unit,*,end=91,err=92) il,jl,ih,jh,x
      LEVEL3 'setting bathymetry in ',il,jl,ih,jh,' to ',x
      do j=jl,jh
         do i=il,ih
            if(imin+ioff .le. i .and. i .le. imax+ioff .and. &
               jmin+joff .le. j .and. j .le. jmax+joff ) then
               i1 = i-ioff
               j1 = j-joff
               H(i1,j1) = x
            end if
         end do
      end do
   end do
   close(unit)

   return

90 LEVEL2 'could not open ',trim(fn),' no bathymetry adjustments done'
   return
91 call getm_error("adjust_bathymetry()","End of file "//trim(fn)//".")
92 call getm_error("adjust_bathymetry()","Error reading "//trim(fn)//".")
   end subroutine adjust_bathymetry
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: adjust_mask - read mask adjustments from file.
!
! !INTERFACE:
   subroutine adjust_mask(fn)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Read mask adjustments from file.
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
   integer                   :: unit = 25 ! kbk
   integer                   :: i,j,k,n
   integer                   :: il,jl,ih,jh
   integer                   :: i1,j1
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Should read in to a buffer at some time to allow for #
   open(unit,file=fn,action='read',status='old',err=90)
   read(unit,*,end=91,err=92) n
   if(n .gt. 1) then
      LEVEL2 'adjusting mask according to:'
      LEVEL3 trim(fn)
   end if
   do k=1,n
      read(unit,*,end=91,err=92) il,jl,ih,jh
      LEVEL3 'masking area ',il,jl,ih,jh
      do j=jl,jh
         do i=il,ih
            if(imin+ioff .le. i .and. i .le. imax+ioff .and. &
               jmin+joff .le. j .and. j .le. jmax+joff ) then
               i1 = i-ioff
               j1 = j-joff
               az(i1,j1) = 0
            end if
         end do
      end do
   end do
   close(unit)

   return

90 LEVEL2 'could not open ',trim(fn),' no mask adjustments done'
   return
91 call getm_error("adjust_mask()","End of file "//trim(fn)//".")
92 call getm_error("adjust_mask()","Error reading "//trim(fn)//".")
   end subroutine adjust_mask
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: print_mask - prints a mask in readable format
!
! !INTERFACE:
   subroutine print_mask(mask)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in), dimension(E2DFIELD) :: mask
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Prints a integer mask in a human readable form.
!
! !REVISION HISTORY:
!
!  22Apr99   Karsten Bolding & Hans Burchard  Initial code.
!
! !LOCAL VARIABLES:
   integer                   :: i,j
!
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
