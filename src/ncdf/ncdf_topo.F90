!$Id: ncdf_topo.F90,v 1.24 2009-12-22 08:44:38 kb Exp $
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
  use domain, only                    : grid_type
  use domain, only                    : xcord,ycord
  use domain, only                    : xxcord,yxcord
  use domain, only                    : dx,dy
  use domain, only                    : xc,yc
  use domain, only                    : xx,yx
  use domain, only                    : dlon,dlat
  use domain, only                    : latc,lonc
  use domain, only                    : latx,lonx
  use domain, only                    : convx,convc
  use domain, only                    : z0_method,z0
  IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
   public                                ncdf_read_topo_file
!
! !DEFINED PARAMETERS:
   REALTYPE, parameter                 :: missing_double =-999.
   REALTYPE, parameter                 :: rearth_default = 6378815
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf (adapted from an earlier version of
!                      Karsten Bolding and Hans Burchard)
!
!  $Log: ncdf_topo.F90,v $
!  Revision 1.24  2009-12-22 08:44:38  kb
!  added conditional compilation checks - Klingbeil
!
!  Revision 1.23  2009-12-10 14:22:52  kb
!  fixed typos - Hofmeister
!
!  Revision 1.22  2009-10-13 13:15:11  kb
!  added psedo-coordinates when grid-type 3 or 4
!
!  Revision 1.21  2009-10-08 16:08:00  kb
!  axes defined in entire domain - cartesian, spherical
!
!  Revision 1.20  2009-10-05 11:40:03  kb
!  fixed behaviour for grid_type=3 and lonx, latx, convx not in topo.nc
!
!  Revision 1.19  2009-09-30 05:32:48  kb
!  fixed calculation of dx, dy when dlon, dlat not present
!
!  Revision 1.18  2009-09-23 10:09:20  kb
!  rewrite of grid-initialisation, optional grid info saved to file, -DSAVE_HALO, updated documentation
!
!  Revision 1.17  2009-09-23 10:04:40  kb
!  reverted to v1.15 - to allow for major update
!
!  Revision 1.15  2007-05-26 15:20:37  kbk
!  print NetCDF version info
!
!  Revision 1.14  2007-02-07 16:32:22  kbk
!  added spatial varying bottom roughness
!
!  Revision 1.13  2006-11-24 09:10:56  frv-bjb
!  Higher accuracy in x0,dx computations
!
!  Revision 1.12  2006-01-29 20:32:34  hb
!  Small LaTeX corrections to source code documentation
!
!  Revision 1.11  2005-11-17 13:50:22  kbk
!  fixes to compile with gfortran
!
!  Revision 1.10  2005/06/17 07:57:46  frv-bjb
!  Bug fix: fail on dlat/lat0 versions
!
!  Revision 1.9  2005/06/14 13:36:01  frv-bjb
!  temporary KBK stop statement deleted
!
!  Revision 1.8  2005/06/10 16:16:41  kbk
!  documentation updated
!
!  Revision 1.7  2005/06/10 16:01:22  kbk
!  test and use real axis before using axis offset+increment method
!
!  Revision 1.6  2005/04/25 09:32:34  kbk
!  added NetCDF IO rewrite + de-stag of velocities - Umlauf
!
! !LOCAL VARIABLES:
  private                                 ncdf_read_2d
!EOP

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncdf_read_topo_file() - read required variables
!
! !INTERFACE:
  subroutine ncdf_read_topo_file(filename)
!
! !USES:
   IMPLICIT NONE
!
! !DESCRIPTION:
!  This routine checks for and opens a NetCDF file with GETM bathymetry and 
!  grid information. The first variable read and checked is $grid\_type$.
!  Subsequent operations depends on the value of $grid\_type$.
!
!  The following steps are done in $ncdf\_read\_topo\_file()$:
!  \begin{itemize}
!     \item[1:] check and open NetCDF file specified by 'filename'
!     \item[2:] read $grid\_type$
!     \item[3:] inquire $bathymetry\_id$
!     \item[4:] some test related to $bathymetry\_id$
!     \item[5:] set local and global index ranges for reading
!     \item[6:] read bathymetry into $H$
!     \item[7:] depending on $grid\_type$ read axes and grid information -
!               also check for optional variables
!     \item[8:] finally - check for and read spatially $z_0$
!  \end{itemize}
!
! !INPUT PARAMETERS:
    character(len=*), intent(in)        :: filename
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
! !LOCAL VARIABLES:
   integer                             :: ncid
   integer                             :: status
   integer                             :: ndims
   integer                             :: dimlen
   integer                             :: id
   integer                             :: bathymetry_id
   integer                             :: xaxis_id=-1
   integer                             :: yaxis_id=-1
   integer, dimension(2)               :: dimidsT(2)
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
!EOP
!-------------------------------------------------------------------------

!  Look for things in the bathymetry file that should be there
!  for all grid types.

   LEVEL2 'using NetCDF version: ',trim(NF90_INQ_LIBVERS())

!  Open file
   status = nf90_open(filename,nf90_nowrite,ncid)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"ncdf_check_grid()",   &
                        "Error opening "//trim(filename)//".")
   endif

!  What kind of grid is it?
   status = nf90_inq_varid(ncid,"grid_type",id)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"ncdf_check_grid()",   &
                        "Could not find 'grid_type' in "//trim(filename)//".")
   endif

   status = nf90_get_var(ncid,id,grid_type)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"ncdf_check_grid()",   &
                        "Could not read 'grid_type' in "//trim(filename)//".")
   endif

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
      case default
         call getm_error("ncdf_check_grid()","Invalid grid type. Choose grid_type=1-4.")
   end select

!  Look for 'bathymetry'
   LEVEL2 'reading bathymetry: '
   status = nf90_inq_varid(ncid,"bathymetry",bathymetry_id)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"ncdf_check_grid()",   &
                        "Could not find 'bathymetry' in "//trim(filename)//".")
   endif

!  Is 'bathymetry' a matrix?
   status = nf90_inquire_variable(ncid,bathymetry_id,ndims=ndims,dimids=dimidsT)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"ncdf_check_grid()",    &
                        "Could not get 'ndims' of 'bathymetry' in "//trim(filename)//".")
   endif

   if (ndims .ne. 2) then
      call getm_error("ncdf_check_grid()","'bathymetry' must have 2 dimensions.")
   endif

!  Is the size of 'bathymetry' consistent?
   status = nf90_inquire_dimension(ncid,dimidsT(1), len = dimlen)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"ncdf_check_grid()",   &
                        "Could not get 'dimlen' of 'bathymetry' in "//trim(filename)//".")
   endif

#ifdef STATIC
   if (dimlen .ne. iextr) then
      call getm_error("ncdf_check_grid()",   &
                      "Length of first dimension in 'bathymetry' inconsistent.")
   endif
#else
!  Get i-dimension for dynamic allocation
   iextr = dimlen
#endif

   status = nf90_inquire_dimension(ncid,dimidsT(2), len = dimlen)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"ncdf_check_grid()",   &
                        "Could not get 'dimlen' of 'bathymetry' in "//trim(filename)//".")
   endif

#ifdef STATIC
    if (dimlen .ne. jextr) then
      call getm_error("ncdf_check_grid()",   &
                      "Length of second dimension in 'bathymetry' inconsistent.")
   endif
#else
!  Get j-dimension for dynamic allocation
   jextr = dimlen
#endif
   LEVEL3 'iextr, jextr: ',iextr,jextr

!  Does the bathymetry have proper axis defined?
!  We will obtain the names of the two dimensions and
!  then inquire if there are variables with the same names - if
!  that is the case they are by NetCDF definition coordinate
!  axis.

   status = nf90_inquire_dimension(ncid,dimidsT(1),name=xaxis_name)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"ncdf_check_grid()",   &
                        "Could not get name associated with dimidsT(1) in "//trim(filename)//".")
   endif

   if ( grid_type .le. 2 ) then
      status = nf90_inq_varid(ncid,xaxis_name,xaxis_id)
      if (status .ne. NF90_NOERR) then
         call netcdf_error(status,"ncdf_check_grid()",   &
                           "Could not get first coordinate name in "//trim(filename)//".")
      endif

      status = nf90_inquire_dimension(ncid,dimidsT(2),name=yaxis_name)
      if (status .ne. NF90_NOERR) then
         call netcdf_error(status,"ncdf_check_grid()",   &
                           "Could not get name associated with dimidsT(2) in "//trim(filename)//".")
      endif
      status = nf90_inq_varid(ncid,yaxis_name,yaxis_id)
      if (status .ne. NF90_NOERR) then
         call netcdf_error(status,"ncdf_check_grid()",   &
                           "Could not get second coordinate name in "//trim(filename)//".")
      endif
      LEVEL3 'axes names:    ',trim(xaxis_name),', ',trim(yaxis_name)
   end if

   ilg = max(imin-HALO+ioff,1); ihg = min(imax+HALO+ioff,iextr)
   jlg = max(jmin-HALO+joff,1); jhg = min(jmax+HALO+joff,jextr)
   iskipl= ilg - (imin-HALO+ioff)
   jskipl= jlg - (jmin-HALO+joff)

!  LOCAL index range for variable to be read
!  (different from GLOBAL range only for parallel runs)
   ill = imin-HALO+iskipl; jll = jmin-HALO+jskipl;
   ihl = ihg-ilg+ill;      jhl = jhg-jlg+jll;

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
                              "Can not find x-axis in "//trim(filename))
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
                              "Error reading y-axis in "//trim(filename))
         end if
#endif
!        y
         status = nf90_inq_varid(ncid,trim(yaxis_name),id)
         if (status .ne. NF90_NOERR) then
            call netcdf_error(status,"ncdf_check_grid()",   &
                              "Can not find y-axis in "//trim(filename))
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
                              "Error reading y-axis in "//trim(filename))
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
                              "Can not find 'x-axis' in "//trim(filename))
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
                              "Error reading x-axis in "//trim(filename))
         end if
#endif
         do j=jmin-HALO,jmax+HALO
            lonc(:,j) = xcord(:)
         end do
!        lat
         status = nf90_inq_varid(ncid,trim(yaxis_name),id)
         if (status .ne. NF90_NOERR) then
            call netcdf_error(status,"ncdf_check_grid()",   &
                              "Can not find 'y-axis' in "//trim(filename))
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
                              "Error reading y-axis in "//trim(filename))
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
                              "Can not find 'xx' in "//trim(filename))
         end if
         call ncdf_read_2d(ncid,id,xx(-1+ill:ihl,-1+jll:jhl),ilg,ihg+1,jlg,jhg+1 )

         status = nf90_inq_varid(ncid,"yx",id)
         if (status .ne. NF90_NOERR) then
            call netcdf_error(status,"ncdf_check_grid()",   &
                              "Can not find 'yx' in "//trim(filename))
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
                              "Can not find 'lonx' in "//trim(filename))
         end if
         call ncdf_read_2d(ncid,id,lonx(-1+ill:ihl,-1+jll:jhl),ilg,ihg+1,jlg,jhg+1)

         status = nf90_inq_varid(ncid,"latx",id)
         if (status .ne. NF90_NOERR) then
            call netcdf_error(status,"ncdf_check_grid()",   &
                              "Can not find 'latx' in "//trim(filename))
         end if
         call ncdf_read_2d(ncid,id,latx(-1+ill:ihl,-1+jll:jhl),ilg,ihg+1,jlg,jhg+1)
         LEVEL4 'done'

         LEVEL4 'convx:'
         status = nf90_inq_varid(ncid,"convx",id)
         if (status .ne. NF90_NOERR) then
            call netcdf_error(status,"ncdf_check_grid()",   &
                              "Can not find 'convx' in "//trim(filename))
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
   if (z0_method .eq. 1) then
      status = nf90_inq_varid(ncid,"z0",id)
      if (status .ne. NF90_NOERR) then
         call netcdf_error(status,"ncdf_check_grid()",   &
                          "Could not find 'z0' in "//trim(filename)//".")
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
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: status
   integer                   :: indx(1)
   integer                   :: i
   REALTYPE                  :: startval,endval
   REALTYPE                  :: expectval,readval,dval
!
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
   status = nf_get_var1_double(ncid,varid,indx,startval)
   if (status .ne. NF_NOERR) then
      call netcdf_error(status,"ncdf_get_grid_dxy()",    &
                        "Could not read first value of " // cordname)
   endif
   indx(1) = iextr
   status = nf_get_var1_double(ncid,varid,indx,endval)
   if (status .ne. NF_NOERR) then
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
      status = nf_get_var1_double(ncid,varid,indx,readval)
      if (status .ne. NF_NOERR) then
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
