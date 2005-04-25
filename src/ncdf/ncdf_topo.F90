!$Id: ncdf_topo.F90,v 1.6 2005-04-25 09:32:34 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdf_topo - check and read topography and grid
!
! !INTERFACE:
module ncdf_topo
!
! !DESCRIPTION:
! This module is responsible for checking and reading all grid-related quantities
! contained in a netCDF file. The grid types available so far are:
! \begin{itemize}
!   \item {grid_type=1}. Cartesian grid with constant grid spacing {\tt dx} and {\tt dy}. 
!         The first $X$-point of the grid may start at {\tt x0} and {\tt y0}. 
!         If these two offsets are not given, they are set to zero.
!   \item {grid_type=2}. Spherical grid with constant grid spacing {\tt dlon} and {\tt dlat}. 
!         The first grid point may start at {\tt lon0} and {\tt lat0}. If these two offsets
!         are not given they are set to zero. Earth's radius, {\tt rearth}, necessary to 
!         construct a spherical grid, may also be given in the input file. 
!         Otherwise it is set to a default value.
!   \item {grid_type=3}. Plane curvilinear grid where the $x$ and $y$ positions of the $X$-points 
!         have to be specified in the bathymetry file.
!   \item {grid_type=4}. Spherical  curvilinear grid where the latitudes and longitudes
!          of the $X$-points have to be specified in the bathymetry file.
! \end{itemize}
!
! GETM requires that the grid positions read from the netCDF file correspond to $X$-points,
!  from which the positions of all other points ($u$, $v$, $T$) can be interpolated 
! straightforwardly. 
! The only exception is the {\tt bathymetry}, which is has to be given on the central $T$-points.
! Since the bathymetery in GETM is defined on $T$-points this avoids unnecessary interpolation 
! of an already polished bathymetry.
!
! Other quantities than those mentioned above may also be specified in the bathymetry file. 
! These  are not mandatory and, if GETM doesn't find them, only a warning will be written.
! To plot data when working with spherical grids, for example, a decision has to be made 
! about the geographical mapping from (lat,lon) to the ($x$,$y$)-plane.  
! Information about both sets of grid points can be read in by GETM - 
! and will be written to the output files used for plotting.
! GETM checks the bathymetery file for corresponding variables called 
! {\tt (latx,lonx)} {\tt ($xx$,$xy$)}. If it is known, also the type of the projection, 
! {\tt proj_type}, and its specifications like the projection center and rotation (
! {\tt proj_lon, proj_lat, proj_rot}) can be given in the input file. The definition
! follows that by the Seagrid grid-generation tool. MATLAB scripts to generate
! simple grids and bathymeteries can be found it the {\tt matlab} subdirectory. 
! Ready-to-run bathymetries for some simple basins can be found among the test cases.
!
! !USES:
  use exceptions
  use domain, only                    : grid_type,proj_type
  use domain, only                    : proj_exists
  use domain, only                    : proj_lon,proj_lat,proj_rot
  use domain, only                    : rearth
  use domain, only                    : latlon_exists,xy_exists
  use domain, only                    : updateXYC,updateLatLonC,updateConvC
  use domain, only                    : dx,dy,x0,y0
  use domain, only                    : dlon,dlat,lon0,lat0
  use domain, only                    : xx,yx,xc,yc
  use domain, only                    : latx,lonx,latc,lonc
  use domain, only                    : convx,convc


  IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
  public                                ncdf_check_grid,ncdf_get_grid
!
! !PUBLIC DATA MEMBERS:
!

! !DEFINED PARAMETERS:
  integer,  parameter                  :: missing_id     =-999
  integer,  parameter                  :: missing_int    =-999
  REALTYPE, parameter                  :: missing_double =-999.
  REALTYPE, parameter                  :: rearth_default = 6378815
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf (adapted from an earlier version of
!                      Karsten Bolding and Hans Burchard)
!
!  $Log: ncdf_topo.F90,v $
!  Revision 1.6  2005-04-25 09:32:34  kbk
!  added NetCDF IO rewrite + de-stag of velocities - Umlauf
!
!
!EOP
!
! !LOCAL VARIABLES:
  private                                 ncdf_read_2d
  integer, private                     :: ncbathy
  integer, private                     :: grid_type_id, proj_type_id
  integer, private                     :: proj_lon_id,proj_lat_id,proj_rot_id
  integer, private                     :: rearth_id
  integer, private                     :: bathymetry_id
  integer, private                     :: dx_id,dy_id,x0_id,y0_id
  integer, private                     :: dlon_id, dlat_id,lon0_id,lat0_id
  integer, private                     :: convx_id,convc_id
  integer, private                     :: latx_id,lonx_id,latc_id,lonc_id
  integer, private                     :: xx_id,yx_id,xc_id,yc_id
  integer, private, dimension(2)       :: dimidsT(2)
  integer, private, dimension(2)       :: dimidsX(2)

  logical, private                     :: latlonx_exists = .true.
  logical, private                     :: latlonc_exists = .true.
  logical, private                     :: xyx_exists     = .true.
  logical, private                     :: xyc_exists     = .true.
  logical, private                     :: convx_exists   = .true.
  logical, private                     :: convc_exists   = .true.


!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncdf_check_grid
!
! !INTERFACE:
  subroutine ncdf_check_grid(filename,iextr,jextr)
!
! !USES:
    IMPLICIT NONE
!
! !DESCRIPTION:
! This routine checks if all mandatory quantities are in the netCDF bathymtery file
! called {\tt filename} and if additonal information about the grid can be found.
! If the run used {\tt STATIC} allocation, the consistency of {\tt iextr} and
! {\tt jextr} with the size of the {\tt bathymetry} in the netCDF file is checked.
! If the run uses {\tt DYNAMICAL} allocation, {\tt iextr} and {\tt jextr} obtained
! from the size of the {\tt bathymetry}. 
!
! If non-mandatory variables are missing, a warning is output and some logical flags
! are set to account for the missing variable.
! 
!
! !INPUT PARAMETERS:
    character(len=*), intent(in)        :: filename
#ifdef STATIC 
    integer, intent(in)                 :: iextr
    integer, intent(in)                 :: jextr
#else
    integer, intent(out)                :: iextr
    integer, intent(out)                :: jextr
#endif
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!
! !LOCAL VARIABLES:
    integer                             :: status
    integer                             :: ndims
    integer                             :: dimlen
!
!-------------------------------------------------------------------------
#include"netcdf.inc"

!  Look for things in the bathymetry file that should be there
!  for all grid types.


!   Open file
    status = nf_open(filename,nf_nowrite,ncbathy)
    if (status .ne. NF_NOERR) then
       call netcdf_error(status,"ncdf_check_grid()",   &
                         "Error opening "//trim(filename)//".")
    endif
  
!   What kind of grid is it?
    status = nf_inq_varid(ncbathy,"grid_type",grid_type_id)
    if (status .ne. NF_NOERR) then
       call netcdf_error(status,"ncdf_check_grid()",   &
                         "Could not find 'grid_type' in "//trim(filename)//".")
    endif

    status = nf_get_var_int(ncbathy,grid_type_id,grid_type)
    if (status .ne. NF_NOERR) then
       call netcdf_error(status,"ncdf_check_grid()",   &
                         "Could not read 'grid_type' in "//trim(filename)//".")
    endif


!   Look for 'bathymetry'
    status = nf_inq_varid(ncbathy,"bathymetry",bathymetry_id)
    if (status .ne. NF_NOERR) then
       call netcdf_error(status,"ncdf_check_grid()",   &
                         "Could not find 'bathymetry' in "//trim(filename)//".")
    endif
    
!   Is 'bathymetry' a matrix?
    status = nf_inq_varndims(ncbathy,bathymetry_id,ndims)
    if (status .ne. NF_NOERR) then
       call netcdf_error(status,"ncdf_check_grid()",    &
                         "Could not get 'ndims' of 'bathymetry' in "//trim(filename)//".")
    endif

    if (ndims.ne.2) then 
       call getm_error("ncdf_check_grid()","'bathymetry' must have 2 dimensions.")
    endif
    
!   Is the size of 'bathymetry' consistent?
    status = nf_inq_vardimid(ncbathy,bathymetry_id,dimidsT)
    if (status .ne. NF_NOERR) then
       call netcdf_error(status,"ncdf_check_grid()",   &
                         "Could not get 'dimidsT' of 'bathymetry' in "//trim(filename)//".")
    endif
    
    status = nf_inq_dimlen(ncbathy,dimidsT(1),dimlen)
    if (status .ne. NF_NOERR) then
       call netcdf_error(status,"ncdf_check_grid()",   &
                         "Could not get 'dimlen' of 'bathymetry' in "//trim(filename)//".")
    endif
    
#ifdef STATIC
    if (dimlen.ne.iextr) then 
       call getm_error("ncdf_check_grid()",   &
                       "Length of first dimension in 'bathymetry' inconsistent.")
    endif
#else
!   Get i-dimension for dynamic allocation
    iextr = dimlen
#endif
    
    status = nf_inq_dimlen(ncbathy,dimidsT(2),dimlen)
    if (status .ne. NF_NOERR) then
       call netcdf_error(status,"ncdf_check_grid()",   &
                         "Could not get 'dimlen' of 'bathymetry' in "//trim(filename)//".")
    endif
    
#ifdef STATIC
    if (dimlen.ne.jextr) then 
       call getm_error("ncdf_check_grid()",   &
                       "Length of second dimension in 'bathymetry' inconsistent.")
    endif
#else
!   Get j-dimension for dynamic allocation
    jextr = dimlen
#endif
    

!   Set grid rotation to zero
!   (may be updated below)
    convx = _ZERO_
    convc = _ZERO_

!   Look for x and y positions of 'X'-points unless grid is Cartesian.
    if (grid_type.ne.1) then 

       status = nf_inq_varid(ncbathy,"xx",xx_id)
       if (status .ne. NF_NOERR) then
          if (grid_type.eq.2) then
             call netcdf_warning(status,"ncdf_check_grid()",   &
                                 "Could not find 'xx' in "//trim(filename)//". I'll try 'xc'.")
          else
             call netcdf_warning(status,"ncdf_check_grid()",   &
                                 "Could not find 'xx' in "//trim(filename)//". Proceeding.")
          endif
          xyx_exists     = .false.
          xx_id          = missing_id
          xx             = missing_double
       endif

       status = nf_inq_varid(ncbathy,"yx",yx_id)
       if (status .ne. NF_NOERR) then
          if (grid_type.eq.2) then
             call netcdf_warning(status,"ncdf_check_grid()",   &
                                 "Could not find 'yx' in "//trim(filename)//". I'll try 'yc'.")
          else
             call netcdf_warning(status,"ncdf_check_grid()",   &
                                 "Could not find 'yx' in "//trim(filename)//". Proceeding.")
          endif
          xyx_exists     = .false.
          yx_id          = missing_id
          yx             = missing_double
       endif

    endif

!   Look for lat,lon,conv at 'X'-points unless grid is equally spaced spherical.
    if (grid_type.ne.2) then

       status = nf_inq_varid(ncbathy,"lonx",lonx_id)
       if (status .ne. NF_NOERR) then
          if (grid_type.eq.1) then
             call netcdf_warning(status,"ncdf_check_grid()",   &
                              "Could not find 'lonx' in "//trim(filename)//". I'll try 'lonc'.")
          else
             call netcdf_warning(status,"ncdf_check_grid()",   &
                              "Could not find 'lonx' in "//trim(filename)//". Proceeding.")
          endif
          latlonx_exists = .false.
          lonx_id       = missing_id
          lonx          = missing_double
       endif

       status = nf_inq_varid(ncbathy,"latx",latx_id)
       if (status .ne. NF_NOERR) then
          if (grid_type.eq.1) then
             call netcdf_warning(status,"ncdf_check_grid()",   &
                               "Could not find 'latx' in "//trim(filename)//". I'll try 'latc'.")
          else
             call netcdf_warning(status,"ncdf_check_grid()",   &
                               "Could not find 'latx' in "//trim(filename)//". Proceeding.")
          endif
          latlonx_exists = .false.
          latx_id       = missing_id
          latx          = missing_double
       endif
          
       status = nf_inq_varid(ncbathy,"convx",convx_id)
       if (status .ne. NF_NOERR) then
          if (grid_type.eq.1) then
             call netcdf_warning(status,"ncdf_check_grid()",   &
                                "Could not find 'convx' in "//trim(filename)//". I'll try 'convc'.")
          else
             call netcdf_warning(status,"ncdf_check_grid()",   &
                                "Could not find 'convx' in "//trim(filename)//". Proceeding.")
          endif
          convx_exists  = .false.
          convx_id      = missing_id
       endif

    endif


!   Has a geographic projection type been specified?
    status = nf_inq_varid(ncbathy,"proj_type",proj_type_id)
    if (status .ne. NF_NOERR) then
       proj_exists    = .false.
       proj_type_id   = missing_id
       proj_type      = 99           ! which means 'unkown'
       call netcdf_warning(status,"ncdf_check_grid()",   &
                           "Could not find 'proj_type' in "//trim(filename)//". Ingnored.")
    endif

!   What parameters are used for the geographic projection?
    if (proj_exists) then 

       status = nf_inq_varid(ncbathy,"proj_lon",proj_lon_id)
       if (status .ne. NF_NOERR) then
          call netcdf_error(status,"ncdf_check_grid()",   &
                           "Could not find 'proj_lon' in "//trim(filename)//".")
       endif

       status = nf_inq_varid(ncbathy,"proj_lat",proj_lat_id)
       if (status .ne. NF_NOERR) then
          call netcdf_error(status,"ncdf_check_grid()",   &
                            "Could not find 'proj_lat' in "//trim(filename)//".")
       endif

       status = nf_inq_varid(ncbathy,"proj_rot",proj_rot_id)
       if (status .ne. NF_NOERR) then
          call netcdf_error(status,"ncdf_check_grid()",   &
                            "Could not find 'proj_rot' in "//trim(filename)//".")
       endif

    endif


!   Look for radius of earth (needed for geographic projection and/or spherical grid)
    status = nf_inq_varid(ncbathy,"rearth",rearth_id)
    if (status .ne. NF_NOERR) then
       if (proj_exists) then
          ! rearth is either needed for the projection ...
          call netcdf_error(status,"ncdf_check_grid()",   &
                            "Could not find 'rearth' in "//trim(filename)//".")
       else
          if ((grid_type.eq.2).or.(grid_type.eq.4)) then
             ! ... or for constructing spherical grids.
             call netcdf_warning(status,"ncdf_check_grid()",   &
                              "Could not find 'rearth' in "//trim(filename)//". Set to default.")
             rearth_id      = missing_id
             rearth         = rearth_default
          endif
       endif
    endif


!   Is all we need for that particular grid_type in the topo-file?
    select case (grid_type)
    case(1)

#if  ( defined(SPHERICAL) || defined(CURVILINEAR) ) 
       call getm_error("ncdf_check_grid()",  &
                       "Cannot use Cartesian grid with SPHERICAL or CURVILINEAR #defined.")
#endif

!      Look for grid spacing
       status = nf_inq_varid(ncbathy,"dx",dx_id)
       if (status .ne. NF_NOERR) then
          call netcdf_error(status,"ncdf_check_grid()",   & 
                            "Could not find 'dx' in "//trim(filename)//".")
       endif
       
       status = nf_inq_varid(ncbathy,"dy",dy_id)
       if (status .ne. NF_NOERR) then
          call netcdf_error(status,"ncdf_check_grid()",   &
                            "Could not find 'dy' in "//trim(filename)//".")
       endif

!      Look for offset
       status = nf_inq_varid(ncbathy,"x0",x0_id)
       if (status .ne. NF_NOERR) then
          call netcdf_warning(status,"ncdf_check_grid()",   &
                              "Could not find 'x0' in "//trim(filename)//". Set to zero." )
          x0_id = missing_id
          x0    = _ZERO_
       endif
       
       status = nf_inq_varid(ncbathy,"y0",y0_id)
       if (status .ne. NF_NOERR) then
          call netcdf_warning(status,"ncdf_check_grid()",    &
                              "Could not find 'y0' in "//trim(filename)//". Set to zero.")
          y0_id = missing_id
          y0    = _ZERO_
       endif
       

!      Try if you can find something more on the T-points
 
       if (.not.latlonx_exists) then

          status = nf_inq_varid(ncbathy,"lonc",lonc_id)
          if (status .ne. NF_NOERR) then
             call netcdf_warning(status,"ncdf_check_grid()",   &
                                 "Could not find 'lonc' in "//trim(filename)//". Ignored.")
             latlonc_exists = .false.
             lonc_id       = missing_id
             lonc          = missing_double
          endif

          status = nf_inq_varid(ncbathy,"latc",latc_id)
          if (status .ne. NF_NOERR) then
             call netcdf_warning(status,"ncdf_check_grid()",   &
                                "Could not find 'latc' in "//trim(filename)//". Ignored.")
             latlonc_exists = .false.
             latc_id       = missing_id
             latc          = missing_double
          endif

       endif


       if (.not.convx_exists) then

          status = nf_inq_varid(ncbathy,"convc",convc_id)
          if (status .ne. NF_NOERR) then
             call netcdf_warning(status,"ncdf_check_grid()",   &
                                "Could not find 'convc' in "//trim(filename)//". Ignored.")
             convc_exists   = .false.
             convc_id       = missing_id
          endif

       endif


       LEVEL3 'Using Cartesian grid.'

 case(2)

#ifndef SPHERICAL
       call getm_error("ncdf_check_grid()",   &
                       "Cannot use spherical grid with SPHERICAL not #defined.")
#endif

!      Look for grid spacing
       status = nf_inq_varid(ncbathy,"dlon",dlon_id)
       if (status .ne. NF_NOERR) then
          call netcdf_error(status,"ncdf_check_grid()",   &
                            "Could not find 'dlon' in "//trim(filename))
       endif

       status = nf_inq_varid(ncbathy,"dlat",dlat_id)
       if (status .ne. NF_NOERR) then
          call netcdf_error(status,"ncdf_check_grid()",   &
                            "Could not find 'dlat' in "//trim(filename))
       endif

!      Look for offset
       status = nf_inq_varid(ncbathy,"lon0",lon0_id)
       if (status .ne. NF_NOERR) then
          call netcdf_warning(status,"ncdf_check_grid()",   &
                              "Could not find 'lon0' in "//trim(filename)//". Set to zero." )
          lon0_id = missing_id
          lon0    = _ZERO_
       endif

       status = nf_inq_varid(ncbathy,"lat0",lat0_id)
       if (status .ne. NF_NOERR) then
          call netcdf_warning(status,"ncdf_check_grid()",    &
               "Could not find 'lat0' in "//trim(filename)//". Set to zero.")
          lat0_id = missing_id 
          lat0    = _ZERO_
       endif
       

!      Try if you can find something on T-points
       if (.not.xyx_exists) then
          
          status = nf_inq_varid(ncbathy,"xc",xc_id)
          if (status .ne. NF_NOERR) then
             call netcdf_warning(status,"ncdf_check_grid()",   &
                                 "Could not find 'xc' in "//trim(filename)//". Ignored.")
             xyc_exists    = .false.
             xc_id         = missing_id
             xc            = missing_double
          endif
          
          status = nf_inq_varid(ncbathy,"yc",yc_id)
          if (status .ne. NF_NOERR) then
             call netcdf_warning(status,"ncdf_check_grid()",   &
                                "Could not find 'yc' in "//trim(filename)//". Ignored.")
             xyc_exists    = .false.
             yc_id         = missing_id
             yc            = missing_double
          endif

       endif

       LEVEL3 'Using spherical grid.'

    case(3)
#ifndef CURVILINEAR
       call getm_error("ncdf_check_grid()",   &
                       "Cannot use curvlinear grid with CURVILINEAR not #defined")
#endif

!      You need at least 'xx', 'xy' and 'convx' for a curvilinear grid
       if (.not.xyx_exists) then
          call getm_error("ncdf_check_grid()",   &
                          "Curvilinear grid needs 'xx' and 'yx' in "//trim(filename)//".")
       endif

       if (.not.convx_exists) then
          call getm_error("ncdf_check_grid()",   &
                          "Curvilinear grid needs 'convx' in "//trim(filename)//".")
       endif

       LEVEL3 'Using plane curvilinear grid.'

    case(4)
#ifndef (CURVILINEAR && SPHERICAL)
       call getm_error("ncdf_check_grid()",                      & 
                     & "Cannot use spherical curvlinear grid with&
                     & CURVILINEAR and SPHERICAL not #defined")
#endif

!      You need at least 'latx', 'lonx', and 'convx' for a curvilinear spherical grid
       if (.not.latlonx_exists) then
          call getm_error("ncdf_check_grid()",    &
                          "Curvilinar spherical grid needs 'lonx' and 'latx' in "//trim(filename)//".")
       endif

       if (.not.convx_exists) then
          call getm_error("ncdf_check_grid()",   &
                          "Spherical curvilinear grid needs 'convx' in "//trim(filename)//".")
       endif


       LEVEL3 'Using spherical curvilinear grid.'

    case default
       call getm_error("ncdf_check_grid()","Invalid grid type. Choose grid_type=1-4.")
    end select
    
    return
  end subroutine ncdf_check_grid
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ncdf_get_grid - reads grid and bathymetry from netCDF file
!
! !INTERFACE:
   subroutine ncdf_get_grid(H,Hland,iextr,jextr,ioff,joff,imin,imax,jmin,jmax)
!
! !USES:

   IMPLICIT NONE
!
! !DESCRIPTION:
! This subroutine reads grid-related quantities and the bathymetry from the 
! netCDF file. It relies on a previous call to {\tt ncdf_check_grid} setting
! all netCDF variable ID's and doing some consistency checks.
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: iextr,jextr,ioff,joff
   integer, intent(in)                 :: imin,imax,jmin,jmax
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: H(E2DFIELD)
   REALTYPE, intent(out)               :: Hland
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: status
   integer                   :: il,ih,jl,jh,iloc,jloc,i,j
   integer                   :: ilocl,iloch,jlocl,jloch
!
!-------------------------------------------------------------------------
#include"netcdf.inc"


!  'Land' value for masking
   Hland = -10.
   H     = Hland

   where ( H .gt. 20000.)
      H = Hland
   end where

!  GLOBAL index range for variable to be read
   il    = max(imin+ioff,1);   ih    = min(imax+ioff,iextr)
   jl    = max(jmin+joff,1);   jh    = min(jmax+joff,jextr)

!  LOCAL index range for variable to be read 
!  (different from GLOBAL range only for parallel runs)  
   ilocl = max(imin-ioff,1);   jlocl = max(jmin-joff,1)
   iloch = ih-il+ilocl;        jloch = jh-jl+jlocl; 

!  Read bathymetry
   call ncdf_read_2d(ncbathy,bathymetry_id,H(ilocl:iloch,jlocl:jloch),il,ih,jl,jh)


!  Read information about (x,y), if available
   if (grid_type.ne.1) then

      if (xyx_exists) then
         call ncdf_read_2d(ncbathy, xx_id, xx(ilocl-1:iloch,jlocl-1:jloch), il,ih+1,jl,jh+1 )
         call ncdf_read_2d(ncbathy, yx_id, yx(ilocl-1:iloch,jlocl-1:jloch), il,ih+1,jl,jh+1 )

!        x-y information exists now
         xy_exists  = .true.

      elseif (xyc_exists.and.(grid_type.eq.2)) then
         
         call ncdf_read_2d(ncbathy,xc_id,xc(ilocl:iloch,jlocl:jloch),il,ih,jl,jh)
         call ncdf_read_2d(ncbathy,yc_id,yc(ilocl:iloch,jlocl:jloch),il,ih,jl,jh)

!        generate positions of X-points
         call c2x(imin,imax,jmin,jmax,xc,xx)
         call c2x(imin,imax,jmin,jmax,yc,yx)


!        x-y information exists now
         xy_exists  = .true.

!        don't update (y,x) at T-points, because it has been read in 
         updateXYC  = .false. 
         
      endif

   endif


!  Read information about (lat,lon) and conv, if available
   if (grid_type.ne.2) then

      if (latlonx_exists) then

         call ncdf_read_2d(ncbathy, lonx_id, lonx(ilocl-1:iloch,jlocl-1:jloch), il,ih+1,jl,jh+1 )
         call ncdf_read_2d(ncbathy, latx_id, latx(ilocl-1:iloch,jlocl-1:jloch), il,ih+1,jl,jh+1 )

!        lat-lon information exists now
         latlon_exists  = .true.

      elseif (latlonc_exists.and.(grid_type.eq.1)) then

         call ncdf_read_2d(ncbathy,lonc_id,lonc(ilocl:iloch,jlocl:jloch),il,ih,jl,jh)
         call ncdf_read_2d(ncbathy,latc_id,latc(ilocl:iloch,jlocl:jloch),il,ih,jl,jh)

!        generate positions of X-points
         call c2x(imin,imax,jmin,jmax,lonc,lonx)
         call c2x(imin,imax,jmin,jmax,latc,latx)


!        lat-lon information exists now
         latlon_exists  = .true.

!        don't update (lat,lon) at T-points, because it has been read in 
         updateLatLonC  = .false. 
         
      endif

      if (convx_exists) then

         call ncdf_read_2d(ncbathy,convx_id,convx(ilocl-1:iloch,jlocl-1:jloch),il,ih+1,jl,jh+1 )

      elseif (convc_exists.and.(grid_type.eq.1)) then

         call ncdf_read_2d(ncbathy,convc_id,convc(ilocl:iloch,jlocl:jloch),il,ih,jl,jh)

!        generate convergence at X-points
         call c2x(imin,imax,jmin,jmax,convc,convx)

!        don't update (lat,lon) at T-points, because it has been read in 
         updateConvC  = .false. 

      endif

   endif


!  Get information about the projection, if available
   if (proj_exists) then

      status = nf_get_var_double(ncbathy,proj_lon_id,proj_lon)
      if (status .ne. NF_NOERR) then
         call netcdf_error(status,"ncdf_get_grid()","Could not read 'proj_lon'.")
      endif
   
      status = nf_get_var_double(ncbathy,proj_lat_id,proj_lat)
      if (status .ne. NF_NOERR) then
         call netcdf_error(status,"ncdf_get_grid()","Could not read 'proj_lat'.")
      endif
   
      status = nf_get_var_double(ncbathy,proj_rot_id,proj_rot)
      if (status .ne. NF_NOERR) then
         call netcdf_error(status,"ncdf_get_grid()","Could not read 'proj_rot'.")
      endif

      status = nf_get_var_double(ncbathy,rearth_id,rearth)
      if (status .ne. NF_NOERR) then
         call netcdf_error(status,"ncdf_get_grid()","Could not read 'rearth'.")
      endif

   endif

!  Get quantities for the specific grid_types
   select case (grid_type)
   case(1)

!     Get grid spacing
      status = nf_get_var_double(ncbathy,dx_id,dx)
      if (status .ne. NF_NOERR) then
         call netcdf_error(status,"ncdf_get_grid()","Could not read 'dx'.")
      endif

      status = nf_get_var_double(ncbathy,dy_id,dy)
      if (status .ne. NF_NOERR) then
         call netcdf_error(status,"ncdf_get_grid()","Could not read 'dy'.")
      endif

!     Get global offsets for x and y      
      if (x0_id.ne.missing_id) then
         status = nf_get_var_double(ncbathy,x0_id,x0)
         if (status .ne. NF_NOERR) then
            call netcdf_error(status,"ncdf_get_grid()","Could not read 'x0'.")
         endif
      endif

      if (y0_id.ne.missing_id) then
         status = nf_get_var_double(ncbathy,y0_id,y0)
         if (status .ne. NF_NOERR) then
            call netcdf_error(status,"ncdf_get_grid()","Could not read 'y0'.")
         endif
      endif
      

    case(2)

!     Get grid spacing
      status = nf_get_var_double(ncbathy,dlon_id,dlon)
      if (status .ne. NF_NOERR) then
         call netcdf_error(status,"ncdf_get_grid()","Could not read 'dlon'.")
      endif

      status = nf_get_var_double(ncbathy,dlat_id,dlat)
      if (status .ne. NF_NOERR) then
         call netcdf_error(status,"ncdf_get_grid()","Could not read 'dlat'.")
      endif


!     Get global offsets for lon and lat      
      if (lon0_id.ne.missing_id) then
         status = nf_get_var_double(ncbathy,lon0_id,lon0)
         if (status .ne. NF_NOERR) then
            call netcdf_error(status,"ncdf_get_grid()","Could not read 'lon0'.")
         endif
      endif

      if (lat0_id.ne.missing_id) then
         status = nf_get_var_double(ncbathy,lat0_id,lat0)
         if (status .ne. NF_NOERR) then
            call netcdf_error(status,"ncdf_get_grid()","Could not read 'lat0'.")
         endif
      endif

      if (rearth_id.ne.missing_id) then
         status = nf_get_var_double(ncbathy,rearth_id,rearth)
         if (status .ne. NF_NOERR) then
            call netcdf_error(status,"ncdf_get_grid()","Could not read 'rearth'.")
         endif
      endif

    case(3)
!       everything already read-in
    case(4)
!       everything already read-in
    case default
       call getm_error("ncdf_get_grid()","Invalid grid type.")
    end select

   return
   end subroutine ncdf_get_grid
!EOC



!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncdf_read_2d
!
! !INTERFACE:
   subroutine ncdf_read_2d(ncid,varid,field,il,ih,jl,jh)
!
! !USES:
   IMPLICIT NONE
!
! !DESCRIPTION:
!  A two-dimensional netCDF variable with specified global range 
! {\tt il < i < ih} and {\tt jl < j < jh} is read into {\tt field}.
! It is checked if the sizes of the fields correspond exactly. 
! When calling this funtions, remember that  FORTRAN netCDF variables 
! start with index 1.
!
! !INPUT PARAMETERS:
   integer,          intent(in)        :: ncid
   integer,          intent(in)        :: varid
   integer,          intent(in)        :: il,ih,jl,jh
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: field(:,:)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!EOP
!
! !LOCAL VARIABLES:
   integer                             :: status
   integer, dimension(2)               :: start
   integer, dimension(2)               :: edges
   integer, dimension(2)               :: ubounds
   character(len=20)                   :: varname
!
!-------------------------------------------------------------------------
#include"netcdf.inc"

   start(1) = il
   start(2) = jl
   edges(1) = ih-il+1
   edges(2) = jh-jl+1

   ubounds =  ubound(field)

   if ((ubounds(1) .ne. edges(1)) .or. ubounds(2) .ne. edges(2) ) then
      call getm_error("ncdf_read_2d()", "Array bounds inconsistent.")
      stop
   endif

   status = nf_inq_varname(ncid,varid,varname)
   if (status .ne. NF_NOERR) then 
      call netcdf_error(status,"read_2d()","Error inquiring name of variable.")
   endif

   status = nf_get_vara_double(ncid,varid,start,edges,field)
   if (status .ne. NF_NOERR) then
      call netcdf_error(status,"read_2d()","Error reading "//trim(varname)//".")
   endif

   return
   end subroutine ncdf_read_2d
!EOC

!-------------------------------------------------------------------------


 end module ncdf_topo

