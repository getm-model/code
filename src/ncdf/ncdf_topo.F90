#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdf_topo() - read bathymetry and grid info (NetCDF)
!
! !INTERFACE:
   module ncdf_topo
!
! !DESCRIPTION:
!  This module reads the bathymetry and grid information required by the
!  module $domain$. The file format is NetCDF and data are read from
!  the file specified as an paramater $ncdf\_read\_topo\_file()$. For a
!  full description of the required variables see the documention for
!  domain. The specific readings are guided by $grid\_type$.

! !USES:
  use netcdf
  use exceptions
  use domain, only                    : have_lonlat,have_xy
  use domain, only                    : iextr,jextr,ioff,joff
  use domain, only                    : imin,imax,jmin,jmax
  use domain, only                    : ilg,ihg,jlg,jhg
  use domain, only                    : ill,ihl,jll,jhl
  use domain, only                    : H, Hland
  use domain, only                    : xcord,ycord
  use domain, only                    : xxcord,yxcord
  use domain, only                    : dx,dy
  use domain, only                    : xc,yc
  use domain, only                    : xx,yx
  use domain, only                    : dlon,dlat
  use domain, only                    : latc,lonc
  use domain, only                    : latx,lonx
  use domain, only                    : convx,convc
  use domain, only                    : bottfric_method,z0
  IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
   public ncdf_open_topo_file,ncdf_read_topo_file
!
! !DEFINED PARAMETERS:
   REALTYPE, parameter                 :: missing_double =-999.
   REALTYPE, parameter                 :: rearth_default = 6378815
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf (adapted from an earlier version of
!                      Karsten Bolding and Hans Burchard)
!
! !LOCAL VARIABLES:
  private                                 ncdf_read_2d
   integer,private                     :: ncid
   integer,private                     :: bathymetry_id
   integer,private,dimension(2)        :: dimidsT
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncdf_open_topo_file() - opens topo file
!
! !INTERFACE:
  subroutine ncdf_open_topo_file(filename,grid_type,iextr,jextr)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: filename
!
! !OUTPUT PARAMETERS:
   integer,intent(out)                 :: grid_type
   integer,intent(out)                 :: iextr,jextr
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!
! !DESCRIPTION:
!  This routine checks for and opens a NetCDF file with GETM bathymetry and
!  grid information. The first variable read and checked is $grid\_type$.
!
!  The following steps are done in $ncdf\_read\_topo\_file()$:
!  \begin{itemize}
!     \item[1:] check and open NetCDF file specified by 'filename'
!     \item[2:] read $grid\_type$
!     \item[3:] inquire $bathymetry\_id$
!     \item[4:] some test related to $bathymetry\_id$
!     \item[5:] set local and global index ranges for reading
!  \end{itemize}
!
! !LOCAL VARIABLES:
   integer                             :: status
   integer                             :: ndims
   integer                             :: id
!EOP
!-------------------------------------------------------------------------

!  Look for things in the bathymetry file that should be there
!  for all grid types.

   LEVEL2 'using NetCDF version: ',trim(NF90_INQ_LIBVERS())

!  Open file
   status = nf90_open(filename,nf90_nowrite,ncid)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"ncdf_open_topo_file()",   &
                        "Error opening "//trim(filename)//".")
   endif

!  What kind of grid is it?
   status = nf90_inq_varid(ncid,"grid_type",id)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"ncdf_open_topo_file()",   &
                        "Could not find 'grid_type' in "//trim(filename)//".")
   endif

   status = nf90_get_var(ncid,id,grid_type)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"ncdf_open_topo_file()",   &
                        "Could not read 'grid_type' in "//trim(filename)//".")
   endif


!  Look for 'bathymetry'
   status = nf90_inq_varid(ncid,"bathymetry",bathymetry_id)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"ncdf_open_topo_file()",   &
                        "Could not find 'bathymetry' in "//trim(filename)//".")
   endif

!  Is 'bathymetry' a matrix?
   status = nf90_inquire_variable(ncid,bathymetry_id,ndims=ndims,dimids=dimidsT)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"ncdf_open_topo_file()",    &
                        "Could not get 'ndims' of 'bathymetry' in "//trim(filename)//".")
   endif

   if (ndims .ne. 2) then
      call getm_error("ncdf_open_topo_file()","'bathymetry' must have 2 dimensions.")
   endif

   status = nf90_inquire_dimension(ncid,dimidsT(1),len=iextr)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"ncdf_open_topo_file()",   &
                        "Could not get 'dimlen' of 'bathymetry' in "//trim(filename)//".")
   endif

   status = nf90_inquire_dimension(ncid,dimidsT(2),len=jextr)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"ncdf_open_topo_file()",   &
                        "Could not get 'dimlen' of 'bathymetry' in "//trim(filename)//".")
   endif

  end subroutine ncdf_open_topo_file
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncdf_read_topo_file() - read required variables
!
! !INTERFACE:
  subroutine ncdf_read_topo_file(grid_type)
!
! !USES:
   IMPLICIT NONE
!
! !DESCRIPTION:
!  This routine relies on a previous call to {\tt ncdf\_open_topo_file}
!  and reads the bathymetry ({\tt H}), coordinate information (depending
!  on {\tt grid_type}) and optionally a spatially variable $z_0$
!  from a NetCDF file.
!
! !INPUT PARAMETERS:
   integer,intent(in)                  :: grid_type
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
! !LOCAL VARIABLES:
   integer                             :: status
   integer                             :: id
   integer                             :: xaxis_id=-1
   integer                             :: yaxis_id=-1
   character*(NF90_MAX_NAME)           :: xaxis_name,yaxis_name
   integer                             :: i,j,n
   integer                             :: iskipl,jskipl
   integer, dimension(1)               :: start
   integer, dimension(1)               :: count
   logical                             :: have_dx=.true.,have_dy=.true.
   logical                             :: have_dlon=.true.,have_dlat=.true.
   logical                             :: have_lon=.false.
   logical                             :: have_lat=.false.
   logical                             :: have_xc=.false.
   logical                             :: have_yc=.false.
   REALTYPE                            :: a(2)
   integer                             :: rc
!EOP
!-------------------------------------------------------------------------

!  Does the bathymetry have proper axis defined?
!  We will obtain the names of the two dimensions and
!  then inquire if there are variables with the same names - if
!  that is the case they are by NetCDF definition coordinate
!  axis.

   if ( grid_type .le. 2 ) then
      status = nf90_inquire_dimension(ncid,dimidsT(1),name=xaxis_name)
      if (status .ne. NF90_NOERR) then
         call netcdf_error(status,"ncdf_check_grid()",   &
                           "Could not get name associated with dimidsT(1) in topo file.")
      endif

      status = nf90_inq_varid(ncid,xaxis_name,xaxis_id)
      if (status .ne. NF90_NOERR) then
         call netcdf_error(status,"ncdf_check_grid()",   &
                           "Could not get first coordinate name in topo file.")
      endif

      status = nf90_inquire_dimension(ncid,dimidsT(2),name=yaxis_name)
      if (status .ne. NF90_NOERR) then
         call netcdf_error(status,"ncdf_check_grid()",   &
                           "Could not get name associated with dimidsT(2) in topo file.")
      endif
      status = nf90_inq_varid(ncid,yaxis_name,yaxis_id)
      if (status .ne. NF90_NOERR) then
         call netcdf_error(status,"ncdf_check_grid()",   &
                           "Could not get second coordinate name in topo file.")
      endif
      LEVEL3 'axes names:    ',trim(xaxis_name),', ',trim(yaxis_name)
   end if

!  Read bathymetry
   call ncdf_read_2d(ncid,bathymetry_id,H(ill:ihl,jll:jhl),ilg,ihg,jlg,jhg)

   select case (grid_type)

      case(1)
!     cartesian - we check for: dx, dy
!     cartesian - we require:   xc, yc
!     cartesian - we check for: lonc, latc, convc
!     cartesian - later we calculate: latu, latv

#if ( defined(SPHERICAL) || defined(CURVILINEAR) )
         call getm_error("ncdf_check_grid()",  &
                         "Cannot use Cartesian grid with SPHERICAL or CURVILINEAR #defined.")
#endif
         LEVEL2 'checking for dx and dy'

         status = nf90_inq_varid(ncid,'dx',id)
         if (status .ne. NF90_NOERR) then
            have_dx=.false.
         else
            status = nf90_get_var(ncid,id,dx)
         end if

         status = nf90_inq_varid(ncid,'dy',id)
         if (status .ne. NF90_NOERR) then
            have_dy=.false.
         else
            status = nf90_get_var(ncid,id,dy)
         end if

         if (have_dx .and. have_dy) then
            LEVEL3 'done'
         else
            LEVEL3 'dx and dy will be calculated from axes information'
         end if

         LEVEL2 'reading/constructing coordinate variables: ', &
                 trim(xaxis_name),', ',trim(yaxis_name)
!        x
         status = nf90_inq_varid(ncid,trim(xaxis_name),id)
         if (status .ne. NF90_NOERR) then
            call netcdf_error(status,"ncdf_check_grid()",   &
                              "Can not find x-axis in topo file")
         end if
         count(1) = 1
         start(1) = 1
         status = nf90_get_var(ncid,id,a(1:1),start = start,count = count)
         start(1) = iextr
         status = nf90_get_var(ncid,id,a(2:2),start = start,count = count)

         if ( .not. have_dx ) dx = (a(2)-a(1))/(iextr-1)
         do i=imin-HALO,imax+HALO
            xcord(i) = a(1) + (i+ioff-1)*dx
         end do
#if 0
!        potentially do checks on consistency of spacing and axis
         start(1) = jlg ; count(1) = ihg-ilg+1
         status = nf90_get_var(ncid,id,xcord(ill:ihl), &
                               start = start, count = count)
         if (status .ne. NF90_NOERR) then
            call netcdf_error(status,"ncdf_check_grid()",   &
                              "Error reading y-axis in topo file")
         end if
#endif
!        y
         status = nf90_inq_varid(ncid,trim(yaxis_name),id)
         if (status .ne. NF90_NOERR) then
            call netcdf_error(status,"ncdf_check_grid()",   &
                              "Can not find y-axis in topo file")
         end if
         count(1) = 1
         start(1) = 1
         status = nf90_get_var(ncid,id,a(1:1),start = start,count = count)
         start(1) = jextr
         status = nf90_get_var(ncid,id,a(2:2),start = start,count = count)
         if ( .not. have_dy ) dy = (a(2)-a(1))/(jextr-1)
         do j=jmin-HALO,jmax+HALO
            ycord(j) = a(1) + (j+joff-1)*dy
         end do
#if 0
!        potentially do checks on consistency of spacing and axis
         start(1) = jlg ; count(1) = jhg-jlg+1
         status = nf90_get_var(ncid,id,ycord(jll:jhl), &
                               start = start, count = count)
         if (status .ne. NF90_NOERR) then
            call netcdf_error(status,"ncdf_check_grid()",   &
                              "Error reading y-axis in topo file")
         end if
#endif
         LEVEL3 'dx= ',dx,', dy= ',dy
#if 0
STDERR xcord(imin:imax)
STDERR ycord(jmin:jmax)
stop
#endif

!        checking for additional fields
         LEVEL2 'checking for additional variable(s): lonc, latc, convc '

         status = nf90_inq_varid(ncid,"lonc",id)
         if (status .ne. NF90_NOERR) then
            LEVEL3 'lonc is not in the file'
         else
            call ncdf_read_2d(ncid,id,lonc(ill:ihl,jll:jhl),ilg,ihg,jlg,jhg)
            LEVEL3 'lonc - OK'
            have_lon = .true.
         end if

         status = nf90_inq_varid(ncid,"latc",id)
         if (status .ne. NF90_NOERR) then
            LEVEL3 'latc is not in the file'
         else
            call ncdf_read_2d(ncid,id,latc(ill:ihl,jll:jhl),ilg,ihg,jlg,jhg)
            LEVEL3 'latc - OK'
            have_lat = .true.
         end if

         if ( have_lon .and. have_lat ) then
            have_lonlat = .true.
         else
            have_lonlat = .false.
         end if

         status = nf90_inq_varid(ncid,"convc",id)
         if (status .ne. NF90_NOERR) then
            LEVEL3 'convc is not in the file'
         else
            call ncdf_read_2d(ncid,id,convc(ill:ihl,jll:jhl),ilg,ihg,jlg,jhg)
            LEVEL3 'convc - OK'
         end if

      case(2)
!     spherical - we check for: dlon, dlat
!     spherical - we require:   lonc, latc
!     spherical - we check for: xc, yc
!     spherical - later we calculate: latu, latv
!

#if !( defined(SPHERICAL) && !defined(CURVILINEAR) )
         call getm_error("ncdf_check_grid()",   &
                          "Cannot use spherical grid with SPHERICAL not #defined or CURVILINEAR #defined.")
#endif
         LEVEL2 'checking for dlon and dlat'

         status = nf90_inq_varid(ncid,'dlon',id)
         if (status .ne. NF90_NOERR) then
            have_dlon=.false.
         else
            status = nf90_get_var(ncid,id,dlon)
         end if

         status = nf90_inq_varid(ncid,'dlat',id)
         if (status .ne. NF90_NOERR) then
            have_dlat=.false.
         else
            status = nf90_get_var(ncid,id,dlat)
         end if

         if (have_dlon .and. have_dlat) then
            LEVEL3 'done'
         else
            LEVEL3 'dlon and dlat will be calculated from axes information'
         end if

         LEVEL2 'reading/constructing coordinate variables: ', &
                 trim(xaxis_name),', ',trim(yaxis_name)
!        lon
         status = nf90_inq_varid(ncid,trim(xaxis_name),id)
         if (status .ne. NF90_NOERR) then
            call netcdf_error(status,"ncdf_check_grid()",   &
                              "Can not find 'x-axis' in topo file")
         end if
         count(1) = 1
         start(1) = 1
         status = nf90_get_var(ncid,id,a(1:1),start = start,count = count)
         start(1) = iextr
         status = nf90_get_var(ncid,id,a(2:2),start = start,count = count)

         if ( .not. have_dlon ) dlon = (a(2)-a(1))/(iextr-1)
         do i=imin-HALO,imax+HALO
            xcord(i) = a(1) + (i+ioff-1)*dlon
         end do
#if 0
!        potentially do checks on consistency of spacing and axis
         start(1) = ilg ; count(1) = ihg-ilg+1
         status = nf90_get_var(ncid,id,xcord(ill:ihl), &
                               start = start, count = count)
         if (status .ne. NF90_NOERR) then
            call netcdf_error(status,"ncdf_check_grid()",   &
                              "Error reading x-axis in topo file")
         end if
#endif
         do j=jmin-HALO,jmax+HALO
            lonc(:,j) = xcord(:)
         end do
!        lat
         status = nf90_inq_varid(ncid,trim(yaxis_name),id)
         if (status .ne. NF90_NOERR) then
            call netcdf_error(status,"ncdf_check_grid()",   &
                              "Can not find 'y-axis' in topo file")
         end if
         count(1) = 1
         start(1) = 1
         status = nf90_get_var(ncid,id,a(1:1),start = start,count = count)
         start(1) = jextr
         status = nf90_get_var(ncid,id,a(2:2),start = start,count = count)

         if ( .not. have_dlat ) dlat = (a(2)-a(1))/(jextr-1)
         do j=jmin-HALO,jmax+HALO
            ycord(j) = a(1) + (j+joff-1)*dlat
         end do
#if 0
!        potentially do checks on consistency of spacing and axis
         start(1) = jlg ; count(1) = jhg-jlg+1
         status = nf90_get_var(ncid,id,ycord(jll:jhl), &
                               start = start, count = count)
         if (status .ne. NF90_NOERR) then
            call netcdf_error(status,"ncdf_check_grid()",   &
                              "Error reading y-axis in topo file")
         end if
#endif
         do i=imin-HALO,imax+HALO
            latc(i,:) = ycord(:)
         end do

         LEVEL3 'dlon= ',dlon,', dlat= ',dlat

!        checking for additional fields
         LEVEL2 'checking for additional variable(s): xc, yc '

         status = nf90_inq_varid(ncid,"xc",id)
         if (status .ne. NF90_NOERR) then
            LEVEL3 'xc is not in the file'
            have_xc = .true.
         else
            call ncdf_read_2d(ncid,id,xc(ill:ihl,jll:jhl),ilg,ihg,jlg,jhg)
            LEVEL3 'xc - OK'
            have_xc = .true.
         end if

         status = nf90_inq_varid(ncid,"yc",id)
         if (status .ne. NF90_NOERR) then
            LEVEL3 'yc is not in the file'
         else
            call ncdf_read_2d(ncid,id,yc(ill:ihl,jll:jhl),ilg,ihg,jlg,jhg)
            LEVEL3 'yc - OK'
            have_yc = .true.
         end if

         if ( have_xc .and. have_yc ) then
            have_xy = .true.
         else
            have_xy = .false.
         end if

      case(3)
!     curvi-linear - we require:   xx, yx
!     curvi-linear - we check for: lonx, latx, convx
!     curvi-linear - later we calculate: lonc, latc, latu, latv

#if !( defined(CURVILINEAR) && !defined(SPHERICAL) )
         call getm_error("ncdf_check_grid()",   &
                         "Cannot use curvlinear grid with CURVILINEAR not #defined or SPHERICAL #defined")
#endif
         LEVEL3 'reading coordinate variables: xx, yx'
         status = nf90_inq_varid(ncid,"xx",id)
         if (status .ne. NF90_NOERR) then
            call netcdf_error(status,"ncdf_check_grid()",   &
                              "Can not find 'xx' in topo file")
         end if
         call ncdf_read_2d(ncid,id,xx(-1+ill:ihl,-1+jll:jhl),ilg,ihg+1,jlg,jhg+1 )

         status = nf90_inq_varid(ncid,"yx",id)
         if (status .ne. NF90_NOERR) then
            call netcdf_error(status,"ncdf_check_grid()",   &
                              "Can not find 'yx' in topo file")
         end if
         call ncdf_read_2d(ncid,id,yx(-1+ill:ihl,-1+jll:jhl),ilg,ihg+1,jlg,jhg+1 )
         LEVEL3 'done'

!        checking for additional fields
         LEVEL2 'checking for additional variable(s): lonx, latx, convx '

         status = nf90_inq_varid(ncid,"lonx",id)
         if (status .ne. NF90_NOERR) then
            LEVEL3 'lonx is not in the file'
         else
            call ncdf_read_2d(ncid,id,lonx(-1+ill:ihl,-1+jll:jhl),ilg,ihg+1,jlg,jhg+1 )
            LEVEL3 'lonx - OK'
            have_lon = .true.
         end if

         status = nf90_inq_varid(ncid,"latx",id)
         if (status .ne. NF90_NOERR) then
            LEVEL3 'latx is not in the file'
         else
            call ncdf_read_2d(ncid,id,latx(-1+ill:ihl,-1+jll:jhl),ilg,ihg+1,jlg,jhg+1 )
            LEVEL3 'latx - OK'
            have_lat = .true.
         end if

         if ( have_lon .and. have_lat ) then
            have_lonlat = .true.
         else
            have_lonlat = .false.
         end if

         status = nf90_inq_varid(ncid,"convx",id)
         if (status .ne. NF90_NOERR) then
            LEVEL3 'convx is not in the file'
         else
            call ncdf_read_2d(ncid,id,convx(-1+ill:ihl,-1+jll:jhl),ilg,ihg+1,jlg,jhg+1 )
            LEVEL3 'convx - OK'
         end if

         LEVEL3 'done'

      case(4)
!     curvi-linear (spherical) - we require:   lonx, latx, convx
!     curvi-linear (spherical) - we check for: xx, yx
!     curvi-linear (spherical) - later we calculate: lonu, latu, lonv, latv and xc, yc
#if !( defined(SPHERICAL) && defined(CURVILINEAR) )
         call getm_error("ncdf_check_grid()",                      &
                       & "Cannot use spherical curvlinear grid with&
                       &  CURVILINEAR or SPHERICAL not #defined")
#endif
         LEVEL3 'reading coordinate variables: lonx, latx'
         status = nf90_inq_varid(ncid,"lonx",id)
         if (status .ne. NF90_NOERR) then
            call netcdf_error(status,"ncdf_check_grid()",   &
                              "Can not find 'lonx' in topo file")
         end if
         call ncdf_read_2d(ncid,id,lonx(-1+ill:ihl,-1+jll:jhl),ilg,ihg+1,jlg,jhg+1)

         status = nf90_inq_varid(ncid,"latx",id)
         if (status .ne. NF90_NOERR) then
            call netcdf_error(status,"ncdf_check_grid()",   &
                              "Can not find 'latx' in topo file")
         end if
         call ncdf_read_2d(ncid,id,latx(-1+ill:ihl,-1+jll:jhl),ilg,ihg+1,jlg,jhg+1)
         LEVEL4 'done'

         LEVEL4 'convx:'
         status = nf90_inq_varid(ncid,"convx",id)
         if (status .ne. NF90_NOERR) then
            call netcdf_error(status,"ncdf_check_grid()",   &
                              "Can not find 'convx' in topo file")
         end if
         call ncdf_read_2d(ncid,id,convx(-1+ill:ihl,-1+jll:jhl),ilg,ihg+1,jlg,jhg+1)
         LEVEL3 'done'

      case default
   end select

   select case (grid_type)
      case(3,4)
!        pseudo coordinates for T- and X-points
         xxcord(imin-HALO-1) = imin-HALO-1+ioff
         do i=imin-HALO,imax+HALO
            xcord(i)  = i+ioff - _HALF_
            xxcord(i) = i+ioff
         end do
         yxcord(jmin-HALO-1) = jmin-HALO-1+joff
         do j=jmin-HALO,jmax+HALO
            ycord(j)  = j+joff - _HALF_
            yxcord(j) = j+joff
         end do
   end select

!  read bottom roughness
   if (bottfric_method .eq. 3) then
      status = nf90_inq_varid(ncid,"z0",id)
      if (status .ne. NF90_NOERR) then
         call netcdf_error(status,"ncdf_check_grid()",   &
                          "Could not find 'z0' in topo file")
      end if

      call ncdf_read_2d(ncid,id,z0(ill:ihl,jll:jhl),ilg,ihg,jlg,jhg)
   end if

    return
  end subroutine ncdf_read_topo_file
!EOC

#if 0
!!!!!!!!!

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: coords_and_grid_spacing
!
! !INTERFACE:
   subroutine coords_and_grid_spacing(ncid,varid,iextr,cordname,x0,dx)
!
! !USES:
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Computes x and dx given that the netcdf file contains the axis
!  (T-point) information.
!  It is assumed that the coordinate values are equidistantly spaced.
!  The equidistance is tested and warnings given if non-equidistant
!  values are noted.
!
!  The routine also works for y, lon, and lat.
!
! !INPUT PARAMETERS:
   integer,      intent(in)             :: ncid
   character(len=*), intent(in)         :: spacing_name
   character(len=*), intent(in)         :: cord_name
   integer,      intent(in)             ::
   character(len=*), intent(in)         :: cordname
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: x0, dx
!
! !REVISION HISTORY:
!  Original author(s): Bjarne Buchmann
!
! !LOCAL VARIABLES:
   integer                   :: status
   integer                   :: indx(1)
   integer                   :: i
   REALTYPE                  :: startval,endval
   REALTYPE                  :: expectval,readval,dval
!EOP
!-------------------------------------------------------------------------

#ifdef DEBUG
   write(debug,*) "ncdf_get_grid_dxy(): working on coordinate " // cordname
#endif
!
! x0 and dx are computed from the T-points.
! Remember (x0, y0, lon0,lat0) are relative to X-points.
! The procedure only works for equidistant grids.
! Necessary to comply with .../domain/domain.F90.
! Iextr (or jextr) is used to find the last entry in the coordinate axis.
!
   indx(1) = 1
   status = nf90_get_var(ncid,varid,startval,indx)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"ncdf_get_grid_dxy()",    &
                        "Could not read first value of " // cordname)
   endif
   indx(1) = iextr
   status = nf90_get_var(ncid,varid,endval,indx)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"ncdf_get_grid_dxy()",    &
                        "Could not read last value of "  // cordname)
   endif

#ifdef DEBUG
   write(debug,*) "ncdf_get_grid_dxy(): Range of " // cordname // " ",   &
        startval,endval
#endif

!
! Compute grid spacing based on first and last value:
!
!
   dx = (endval-startval)/(iextr-1)

!
! Compute x0 as dx/2 before first read value.
! x0 should be an X-point (not a T-point).
!
   x0 = startval - 0.50*dx

!
! Test that the read values are approximately equidistantly spaced.
! This implementation could be faster if we read the entire array, but the
! present implementation is low on memory/allocation. Also, it is really
! fast as it is 1D and executed once, so there should be no performance
! gain by reading the entire array.
!
   do i=1,iextr
! Note that startval and x0 no longer match at this point,
      expectval = startval + dx * (i-1)
      indx(1) = i
      status = nf90_get_var(ncid,varid,indx,readval)
      if (status .ne. NF90_NOERR) then
         call netcdf_error(status,"ncdf_get_grid_dxy()",   &
                           "Could not read one value of " // cordname)
      endif
      dval = abs(expectval-readval)
! Compare with a fairly lax criterion:
      if (dval .gt. 0.1 * dx) then
         LEVEL1 "Warning: Non-equidistant grid detected for " // cordname
         LEVEL1 "    Read value no. ",i
         LEVEL1 "    Expected value ",expectval
         LEVEL1 "    Actually read  ",readval
! Dont bother checking the rest.
! All values are set, so we might as well return from here.
         return
      end if
   end do

   return
   end subroutine coords_and_grid_spacing
!EOC

!-------------------------------------------------------------------------
!!!!!!!!!
#endif

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncdf_read_2d() - generic reading routine
!
! !INTERFACE:
   subroutine ncdf_read_2d(ncid,varid,field,il,ih,jl,jh)
!
! !USES:
   IMPLICIT NONE
!
! !DESCRIPTION:
!  A two-dimensional netCDF variable with specified global range
!  {\tt il < i < ih} and {\tt jl < j < jh} is read into {\tt field}.
!  It is checked if the sizes of the fields correspond exactly.
!  When calling this funtions, remember that  FORTRAN netCDF variables
!  start with index 1.
!
! !INPUT PARAMETERS:
   integer,          intent(in)        :: ncid
   integer,          intent(in)        :: varid
   integer,          intent(in)        :: il,ih,jl,jh
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(inout)             :: field(:,:)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
! !LOCAL VARIABLES:
   integer                             :: status
   integer, dimension(2)               :: start
   integer, dimension(2)               :: count
   integer, dimension(2)               :: ubounds
   character(len=20)                   :: varname
!EOP
!-------------------------------------------------------------------------

   start(1) = il
   start(2) = jl
   count(1) = ih-il+1
   count(2) = jh-jl+1

   ubounds =  ubound(field)

   if ((ubounds(1) .ne. count(1)) .or. ubounds(2) .ne. count(2) ) then
      call getm_error("ncdf_read_2d()", "Array bounds inconsistent.")
      stop
   endif

#if 0
   status = nf90_inquire_variable(ncid,varid,name=varname)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"read_2d()","Error inquiring name of variable.")
   endif
#endif

   status = nf90_get_var(ncid,varid,field,start = start,count = count)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"read_2d()","Error reading "//trim(varname)//".")
   endif

   return
   end subroutine ncdf_read_2d
!EOC

 end module ncdf_topo

!-----------------------------------------------------------------------
! Copyright (C) 2009 - Karsten Bolding and Hans Burchard
!-----------------------------------------------------------------------
