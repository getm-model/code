!$Id: grid_interpol.F90,v 1.5 2003-06-30 05:45:26 kbk Exp $
#include "cppdefs.h"
#ifndef HALO
#define HALO 0
#endif

! This needs to be made smarter..
! If USE_VALID_LON_LAT_ONLY is set the inerpolation is only done for
! valid lat,lon specifications - however, sometimes it is only desirable
! to calculate using the mask information e.g. salinity and temperature
! climatologies.
! lat and lon should be initialised to something < -1000.
! Comment out the following line and re-compile then the mask-method is
! used
#define USE_VALID_LON_LAT_ONLY

!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  grid_interpol - various interpolation routines
!
! !INTERFACE:
   module grid_interpol
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
!
! !PRIVATE DATA MEMBERS:
   REALTYPE, parameter       :: pi=3.1415926535897932384626433832795029
   REALTYPE, parameter       :: deg2rad=pi/180.,rad2deg=180./pi
   REALTYPE, parameter       :: earth_radius=6370.9490e3

   integer :: il,ih,jl,jh
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: grid_interpol.F90,v $
!  Revision 1.5  2003-06-30 05:45:26  kbk
!  do interpolation in rotated space
!
!  Revision 1.4  2003/06/18 08:49:53  kbk
!  southpole has length 3: lat, lon, rotation
!
!  Revision 1.3  2003/04/23 12:02:43  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 15:25:06  kbk
!  parallel support
!
!  Revision 1.1.1.1  2002/05/02 14:01:21  gotm
!  recovering after CVS crash
!
!  Revision 1.3  2001/09/30 09:06:00  bbh
!  Cleaned up
!
!  Revision 1.2  2001/07/26 13:51:35  bbh
!  Grid interpolation (including rotated grid) now works
!
!  Revision 1.1  2001/05/25 18:52:21  bbh
!  Added grid_interpol
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_grid_interpol - initialise grid interpolation.
!
! !INTERFACE:
   subroutine init_grid_interpol(imin,imax,jmin,jmax,mask,      &
                         olon,olat,met_lon,met_lat,southpole,   &
                         gridmap,t,u,met_mask)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  To be written.
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: imin,imax,jmin,jmax
   integer, intent(in)                 :: mask(-HALO+1:,-HALO+1:)
   REALTYPE, intent(in)                :: olon(-HALO+1:,-HALO+1:)
   REALTYPE, intent(in)                :: olat(-HALO+1:,-HALO+1:)
   REALTYPE, intent(in)                :: met_lon(:),met_lat(:)
   REALTYPE, intent(in)                :: southpole(3)
   integer, optional, intent(in)       :: met_mask(:,:)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out) :: t(-HALO+1:,-HALO+1:),u(-HALO+1:,-HALO+1:)
   integer, intent(out) :: gridmap(-HALO+1:,-HALO+1:,1:)
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   integer                   :: rc
   integer                   :: i,j
   REALTYPE                  :: x(4),y(4)
   REALTYPE                  :: z
   REALTYPE                  :: xr,yr,zr
   REALTYPE, allocatable     :: beta(:,:)
!EOP
!-------------------------------------------------------------------------
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_grid_interpol() # ',Ncall
#endif
   LEVEL1 'init_grid_interpol'

   il = imin; ih = imax
   jl = jmin; jh = jmax

#ifdef USE_VALID_LON_LAT_ONLY
   LEVEL2 'interpolates for all valid lat-lon'
#else
   LEVEL2 'interpolates only when mask > 0'
#endif

   if(southpole(3) .ne. _ZERO_ ) then
      FATAL 'southpole(3) (rotation) is not coded yet'
      stop 'init_grid_interpol'
   endif

   if(southpole(1) .ne. _ZERO_ .or. southpole(2) .ne. -90.) then
      LEVEL2 'source field domain (rotated coordinates):'
   else
      LEVEL2 'source field domain (geo-graphical coordinates):'
   end if
   if(met_lon(1) .lt. met_lon(size(met_lon))) then
      LEVEL3 'lon: ',met_lon(1),met_lon(size(met_lon))
   else
      LEVEL3 'lon: ',met_lon(size(met_lon)),met_lon(1)
   end if
   if(met_lat(1) .lt. met_lat(size(met_lat))) then
      LEVEL3 'lat: ',met_lat(1),met_lat(size(met_lat))
   else
      LEVEL3 'lon: ',met_lat(size(met_lat)),met_lat(1)
   end if

   LEVEL2 'target field domain:'
   LEVEL3 'lon: ',olon(1,1),olon(imax,jmax)
   LEVEL3 'lat: ',olat(1,1),olat(imax,jmax)

   allocate(beta(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_gridinterpol: Error allocating memory (beta)'

   if (present(met_mask)) then
      call interpol_coefficients(mask,southpole, &
                   olon,olat,met_lon,met_lat,beta,gridmap,t,u,met_mask)
   else
      call interpol_coefficients(mask,southpole, &
                   olon,olat,met_lon,met_lat,beta,gridmap,t,u)
   end if

#ifdef DEBUG
   write(debug,*) 'leaving init_grid_interpol()'
   write(debug,*)
#endif
   return
   end subroutine init_grid_interpol
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_grid_interpol - do grid interpolation.
!
! !INTERFACE:
   subroutine do_grid_interpol(mask,ifield,gridmap,t,u,ofield)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  To be written.
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: mask(-HALO+1:,-HALO+1:)
   REALTYPE, intent(in)                :: ifield(:,:)
   integer, intent(in)                 :: gridmap(-HALO+1:,-HALO+1:,1:)
   REALTYPE, intent(in)                :: t(-HALO+1:,-HALO+1:)
   REALTYPE, intent(in)                :: u(-HALO+1:,-HALO+1:)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: ofield(-HALO+1:,-HALO+1:)
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   integer                   :: i,j
   integer                   :: i1,i2,j1,j2
   REALTYPE                  :: d11,d21,d22,d12
!
!EOP
!-------------------------------------------------------------------------
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_grid_interpol() # ',Ncall
#endif

   do j=jl,jh
      do i=il,ih
         i1 = gridmap(i,j,1)
         j1 = gridmap(i,j,2)
         if(i1 .gt. -999 .and. j1 .gt. -999) then
            if(i1 .ge. size(ifield,1) .or. j1 .ge. size(ifield,2)) then
               ofield(i,j) = ifield(i1,j1)
            else
               i2 = i1+1
               j2 = j1+1
               d11 = (_ONE_-t(i,j))*(_ONE_-u(i,j))
               d21 = t(i,j)*(_ONE_-u(i,j))
               d22 = t(i,j)*u(i,j)
               d12 = (_ONE_-t(i,j))*u(i,j)
               ofield(i,j) = d11*ifield(i1,j1)+d21*ifield(i2,j1)+      &
                             d22*ifield(i2,j2)+d12*ifield(i1,j2)
            end if
         else
            ofield(i,j) = _ZERO_
         end if
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving do_grid_interpol()'
   write(debug,*)
#endif
   return
   end subroutine do_grid_interpol
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_3d_grid_interpol - do grid interpolation.
!
! !INTERFACE:
   subroutine do_3d_grid_interpol(ifield,gridmap,t,u,ofield)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  To be written.
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: ifield(:,:,:)
   integer, intent(in)                 :: gridmap(-HALO+1:,-HALO+1:,1:)
   REALTYPE, intent(in)                :: t(-HALO+1:,-HALO+1:)
   REALTYPE, intent(in)                :: u(-HALO+1:,-HALO+1:)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: ofield(-HALO+1:,-HALO+1:,0:)
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   integer                   :: i,j
   integer                   :: i1,i2,j1,j2
   REALTYPE                  :: d11,d21,d22,d12
!
!EOP
!-------------------------------------------------------------------------
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_3d_grid_interpol() # ',Ncall
#endif

   do j=jl,jh
      do i=il,ih
         i1 = gridmap(i,j,1)
         i2 = i1+1
         j1 = gridmap(i,j,2)
         j2 = j1+1
         d11 = (_ONE_-t(i,j))*(_ONE_-u(i,j))
         d21 = t(i,j)*(_ONE_-u(i,j))
         d22 = t(i,j)*u(i,j)
         d12 = (_ONE_-t(i,j))*u(i,j)
         ofield(i,j,1:) = d11*ifield(i1,j1,:)+d21*ifield(i2,j1,:)+     &
                          d22*ifield(i2,j2,:)+d12*ifield(i1,j2,:)
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving do_3d_grid_interpol()'
   write(debug,*)
#endif
   return
   end subroutine do_3d_grid_interpol
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: to_rotated_lat_lon - from geographical to  rotated lat-lon
!
! !INTERFACE:
   subroutine to_rotated_lat_lon(southpole,alon,alat,rlon,rlat,beta)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  To be written.
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: southpole(3)
   REALTYPE, intent(in)                :: alon
   REALTYPE, intent(in)                :: alat
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: rlon
   REALTYPE, intent(out)               :: rlat
   REALTYPE, intent(out)               :: beta
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   REALTYPE                  :: sinphis,cosphis
   REALTYPE                  :: alpha,cosalpha,sinalpha
   REALTYPE                  :: phi,sinphi,cosphi
   REALTYPE                  :: SA,CA,SB,CB
!
!EOP
!-------------------------------------------------------------------------
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'to_rotated_lat_lon() # ',Ncall
#endif

!   LEVEL1 'to_rotated_lat_lon'

   sinphis=sin(deg2rad*southpole(1))
   cosphis=cos(deg2rad*southpole(1))

   alpha = deg2rad*(alon-southpole(2))
   cosalpha = cos(alpha)
   sinalpha = sin(alpha)

   phi = deg2rad*alat
   sinphi = sin(phi)
   cosphi = cos(phi)

   rlat = asin(-sinphis*sinphi-cosphis*cosphi*cosalpha)*rad2deg

   SA = sinalpha*cosphi
   CA = cosphis*sinphi-sinphis*cosphi*cosalpha
   rlon = atan2(SA,CA)*rad2deg

   SB =  sinalpha*cosphis
   CB = -sinphis*cosphi+cosphis*sinphi*cosalpha
   beta = atan2(SB,CB)

#ifdef DEBUG
   write(debug,*) 'Leaving to_rotated_lat_lon()'
   write(debug,*)
#endif
   return
   end subroutine to_rotated_lat_lon
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: from_rotated_lat_lon - from rotated lat-lon to geographical
!
! !INTERFACE:
   subroutine from_rotated_lat_lon(southpole,rlon,rlat,alon,alat,beta)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  To be written.
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: southpole(3)
   REALTYPE, intent(in)                :: rlon
   REALTYPE, intent(in)                :: rlat
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: alon
   REALTYPE, intent(out)               :: alat
   REALTYPE, intent(out)               :: beta
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   REALTYPE                  :: sinphis,cosphis
   REALTYPE                  :: lambda,coslambda,sinlambda
   REALTYPE                  :: phi,sinphi,cosphi
   REALTYPE                  :: SA,CA,SB,CB
!
!EOP
!-------------------------------------------------------------------------
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'to_rotated_lat_lon() # ',Ncall
#endif

   sinphis=sin(deg2rad*southpole(1))
   cosphis=cos(deg2rad*southpole(1))

   lambda = deg2rad*rlon
   coslambda = cos(lambda)
   sinlambda = sin(lambda)

   phi = deg2rad*rlat
   sinphi = sin(phi)
   cosphi = cos(phi)

   alat = asin(-sinphis*sinphi+cosphis*cosphi*coslambda)*rad2deg

   SA = -sinlambda*cosphi
   CA = cosphis*sinphi+sinphis*cosphi*coslambda
   alon = southpole(2)+atan2(-SA,-CA)*rad2deg

   SB =  -sinlambda*cosphis
   CB = -sinphis*cosphi-cosphis*sinphi*coslambda
   beta = atan2(SB,CB)

#ifdef DEBUG
   write(debug,*) 'Leaving to_rotated_lat_lon()'
   write(debug,*)
#endif
   return
   end subroutine from_rotated_lat_lon
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: spherical_dist - calculate distances on a sphere
!
! !INTERFACE:
   REALTYPE function spherical_dist(radius,lon1,lat1,lon2,lat2)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Calculates the distance - in meters - between the two poinst specified
!  by (lon1,lat1) and (lon2,lat2). Radius is the radius of the sphere -
!  usually the radius of the earth.
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: radius,lon1,lat1,lon2,lat2
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
   REALTYPE                  :: a,b,c
!
!EOP
!-------------------------------------------------------------------------
   a = sin(deg2rad*0.5*(lat2-lat1))
   b = sin(deg2rad*0.5*(lon2-lon1))
   c = a*a + cos(deg2rad*lat1)*cos(deg2rad*lat2)*b*b
   spherical_dist = radius*2.0*atan2(sqrt(c),sqrt(1-c))
   end function spherical_dist
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: interpol_coefficients - set up interpolation coeffcients
!
! !INTERFACE:
   subroutine interpol_coefficients(mask,sp,olon,olat,met_lon,met_lat, &
                                    beta,gridmap,t,u,met_mask)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  To be written.
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: mask(-HALO+1:,-HALO+1:)
   REALTYPE, intent(in)                :: sp(3)
   REALTYPE, intent(in)                :: olon(-HALO+1:,-HALO+1:)
   REALTYPE, intent(in)                :: olat(-HALO+1:,-HALO+1:)
   REALTYPE, intent(in)                :: met_lon(:),met_lat(:)
   integer, optional, intent(in)       :: met_mask(:,:)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: beta(-HALO+1:,-HALO+1:)
   integer, intent(out)                :: gridmap(-HALO+1:,-HALO+1:,1:)
   REALTYPE, intent(out)               :: t(-HALO+1:,-HALO+1:)
   REALTYPE, intent(out)               :: u(-HALO+1:,-HALO+1:)
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   logical                   :: rotated_grid
   integer                   :: i,j,im,jm
   REALTYPE                  :: alon,alat
   REALTYPE                  :: x,y,lon1,lat1,lon2,lat2
   integer                   :: ngood
   logical                   :: outside=.false.
   integer                   :: max_i,max_j
   logical                   :: increasing_lat,increasing_lon
!EOP
!-------------------------------------------------------------------------
!  first find the lower left (im,jm) in the m-grid which coresponds to
!  olon(i,j), olat(i,j)

   max_i = size(met_lon)
   max_j = size(met_lat)

   increasing_lon = met_lon(1) .lt. met_lon(max_i)
   increasing_lat = met_lat(1) .lt. met_lat(max_j)

   if(sp(1) .ne. _ZERO_ .or. sp(2) .ne. -90.) then
      rotated_grid = .true.
   else
      rotated_grid = .false.
   end if

   do j=jl,jh
      do i=il,ih
#ifdef USE_VALID_LON_LAT_ONLY
         if(olon(i,j) .gt. -1000. .and. olat(i,j) .gt. -1000.) then
#else
         if(mask(i,j) .ge. 1) then
#endif
            if (rotated_grid) then
               call to_rotated_lat_lon(sp,olon(i,j),olat(i,j), &
                                       alon,alat,beta(i,j))
            else
               alon = olon(i,j)
               alat = olat(i,j)
               beta(i,j) = _ZERO_
            end if

            if (increasing_lon) then
               if(met_lon(1) .le. alon .and. alon .le. met_lon(max_i)) then
                  do im=1,max_i
                     if(met_lon(im) .gt. alon) EXIT
                  end do
                  gridmap(i,j,1) = im-1
              else
                  outside = .true.
	       end if
            else
            endif

            if (increasing_lat) then
               if(met_lat(1) .le. alat .and. alat .le. met_lat(max_j)) then
                  do jm=1,max_j
                     if(met_lat(jm) .gt. alat) EXIT
                  end do
                  gridmap(i,j,2) = jm-1
	       else
                  outside = .true.
	       end if
            else
            endif
#if 0
!KBK - test this
	    if(im .gt. 2) then
               gridmap(i,j,1) = im-1
            else
               gridmap(i,j,1) = im
            end if
	    if(jm .gt. 2) then
               gridmap(i,j,2) = jm-1
            else
               gridmap(i,j,2) = jm
            end if
! ======
#endif
         end if
      end do
   end do

   if(outside) then
      STDERR 'WARNING - interpol_coefficients: Some points out side the area'
   end if

!  then calculated the t and u coefficients - via distances - the point of
!  interest is (x,y)
   do j=jl,jh
      do i=il,ih
#ifdef USE_VALID_LON_LAT_ONLY
         if(olon(i,j) .gt. -1000. .and. olat(i,j) .gt. -1000.) then
#else
         if(mask(i,j) .ge. 1) then
#endif
            if (rotated_grid) then
               call to_rotated_lat_lon(sp,olon(i,j),olat(i,j), &
	                               x,y,beta(i,j))
            else
	       x = olon(i,j)
	       y = olat(i,j)
            end if
            im = gridmap(i,j,1)
            jm = gridmap(i,j,2)
            if(im .ge. max_i .or. jm .ge. max_j) then
               t(i,j) = _ZERO_
               u(i,j) = _ZERO_
            else
               if(present(met_mask)) then
                  ngood = met_mask(im  ,jm  )+met_mask(im+1,jm  )+ &
                          met_mask(im+1,jm+1)+met_mask(im  ,jm+1)
               else
                  ngood = 4
               end if
               t(i,j) = _ZERO_
               u(i,j) = _ZERO_
ngood = 4
               select case (ngood)
                  case (0)
                  case (1,2,3)
                     t(i,j) = _ZERO_
!                    if(met_mask(im,jm) .eq. 0 .or. met_mask(im,jm+1) .eq. 0 ) &
!                         t(i,j) = _ONE_
                     u(i,j) = _ZERO_
!                    if(met_mask(im,jm) .eq. 0 .or. met_mask(im+1,jm) .eq. 0 ) &
!                          u(i,j) = _ONE_
                  case (4)
                     lon1 = met_lon(im)
                     lat1 = met_lat(jm)
                     lon2 = met_lon(im+1)
                     lat2 = met_lat(jm+1)
                     t(i,j) = spherical_dist(earth_radius,x,lat1,lon1,lat1)/ &
                              spherical_dist(earth_radius,lon1,lat1,lon2,lat1)
                     u(i,j) = spherical_dist(earth_radius,lon1,y,lon1,lat1)/ &
                              spherical_dist(earth_radius,lon1,lat1,lon1,lat2)
                  case default
               end select
            end if
         else if (mask(i,j) .gt. 0) then
            FATAL 'Could not find coefficients for all water points'
            FATAL 'ocean(i,j) = ',i,j
            stop 'interpol_coefficients()'
         end if
      end do
   end do

   end subroutine interpol_coefficients
!EOC

!-----------------------------------------------------------------------

   end module grid_interpol

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Karsten Bolding and Hans Burchard               !
!-----------------------------------------------------------------------
