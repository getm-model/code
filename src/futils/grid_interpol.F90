!$Id: grid_interpol.F90,v 1.4 2003-06-18 08:49:53 kbk Exp $
#include "cppdefs.h"
#ifndef HALO
#define HALO 0
#endif
#if 1
#define USE_VALID_LON_LAT_ONLY
#endif
#ifdef USE_VALID_LON_LAT_ONLY
#define USE_GRIDMAP_UNDEF
#endif
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

   REALTYPE, allocatable, dimension(:,:) :: rlon,rlat,beta
   integer :: il,ih,jl,jh
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: grid_interpol.F90,v $
!  Revision 1.4  2003-06-18 08:49:53  kbk
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
   REALTYPE, intent(in)                :: met_lon(:,:),met_lat(:,:)
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
   logical                   :: rotated_grid
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

   if(southpole(3) .ne. _ZERO_ ) then
      FATAL 'southpole(3) (rotation) is not coded yet'
      stop 'init_grid_interpol'
   endif

   if(southpole(1) .ne. _ZERO_ .or. southpole(2) .ne. -90.) then
      rotated_grid = .true.
   else
      rotated_grid = .false.
   end if
   LEVEL2 'southpole=',southpole
   LEVEL2 ' --> rotated_grid=',rotated_grid
   if( rotated_grid ) then
#if 0
      allocate(rot_lon(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_gridinterpol: Error allocating memory (rot_lon)'

      allocate(rot_lat(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_gridinterpol: Error allocating memory (rot_lat)'

      allocate(beta(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_gridinterpol: Error allocating memory (beta)'

      call rotated_lat_lon(mask,southpole,met_lon,met_lat,rot_lon,rot_lat,beta)
      call interpol_coefficients(mask,olon,olat,rot_lon,rot_lat,gridmap,t,u,met_mask)
#endif
STDERR 'Check this Karsten - rotated grid'
stop 'init_grid_interpol.F90'
   else
      call interpol_coefficients(mask,olon,olat,met_lon,met_lat,gridmap,t,u,met_mask)
   endif

#ifdef DEBUG
   write(debug,*) 'Leaving init_grid_interpol()'
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
#ifdef USE_GRIDMAP_UNDEF
         i1 = gridmap(i,j,1)
         j1 = gridmap(i,j,2)
         if(i1 .gt. -999 .and. j1 .gt. -999) then
#else
         if(mask(i,j) .ge. 1) then
            i1 = gridmap(i,j,1)
            j1 = gridmap(i,j,2)
#endif
            d11 = (_ONE_-t(i,j))*(_ONE_-u(i,j))
            if(i1 .ge. size(ifield,1) .or. j1 .ge. size(ifield,2)) then
               ofield(i,j) = ifield(i1,j1)
            else
               i2 = i1+1
               j2 = j1+1
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

#if 0
i = 50 ; j = 50
i1 = gridmap(i,j,1)
i2 = i1+1
j1 = gridmap(i,j,2)
j2 = j1+1
d11 = (_ONE_-t(i,j))*(_ONE_-u(i,j))
d21 = t(i,j)*(_ONE_-u(i,j))
d22 = t(i,j)*u(i,j)
d12 = (_ONE_-t(i,j))*u(i,j)
STDERR gridmap(i,j,1),gridmap(i,j,2)
STDERR d11,d21,d22,d12
STDERR ifield(i1,j1)
STDERR ifield(i2,j1)
STDERR ifield(i2,j2)
STDERR ifield(i1,j2)
STDERR ofield(i,j)
stop 'kbk'
#endif

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
! !IROUTINE: rotated_lat_lon - get the rotated lat-lon
!
! !INTERFACE:
   subroutine rotated_lat_lon(southpole,met_lon,met_lat,rot_lon,rot_lat,beta)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  To be written.
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: southpole(2)
   REALTYPE, intent(in)                :: met_lon(-HALO+1:,-HALO+1:)
   REALTYPE, intent(in)                :: met_lat(-HALO+1:,-HALO+1:)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: rot_lon(-HALO+1:,-HALO+1:)
   REALTYPE, intent(out)               :: rot_lat(-HALO+1:,-HALO+1:)
   REALTYPE, intent(out)               :: beta(-HALO+1:,-HALO+1:)
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   integer                   :: i,j
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
   write(debug,*) 'rotated_lat_lon() # ',Ncall
#endif

   LEVEL1 'rotated_lat_lon'

   sinphis=sin(deg2rad*southpole(1))
   cosphis=cos(deg2rad*southpole(1))

   do j=jl,jh
      do i=il,ih
         alpha = deg2rad*(met_lon(i,j)-southpole(2))
         cosalpha = cos(alpha)
         sinalpha = sin(alpha)

         phi = deg2rad*met_lat(i,j)
         sinphi = sin(phi)
         cosphi = cos(phi)

         rot_lat(i,j) = asin(-sinphis*sinphi-cosphis*cosphi*cosalpha)*rad2deg

         SA = sinalpha*cosphi
         CA = cosphis*sinphi-sinphis*cosphi*cosalpha
         rot_lon(i,j) = atan2(SA,CA)*rad2deg

         SB =  sinalpha*cosphis
         CB = -sinphis*cosphi+cosphis*sinphi*cosalpha
         beta(i,j) = atan2(SB,CB)
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving rotated_lat_lon()'
   write(debug,*)
#endif
   return
   end subroutine rotated_lat_lon
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
   subroutine interpol_coefficients(mask,olon,olat,met_lon,met_lat,gridmap,t,u,met_mask)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  To be written.
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: mask(-HALO+1:,-HALO+1:)
   REALTYPE, intent(in)                :: olon(-HALO+1:,-HALO+1:)
   REALTYPE, intent(in)                :: olat(-HALO+1:,-HALO+1:)
   REALTYPE, intent(in)                :: met_lon(:,:),met_lat(:,:)
   integer, optional, intent(in)       :: met_mask(:,:)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   integer, intent(out)                :: gridmap(-HALO+1:,-HALO+1:,1:)
   REALTYPE, intent(out)               :: t(-HALO+1:,-HALO+1:)
   REALTYPE, intent(out)               :: u(-HALO+1:,-HALO+1:)
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   integer                   :: i,j,im,jm
   REALTYPE                  :: x,y,lon1,lat1,lon2,lat2
   integer                   :: ngood
   logical                   :: outside=.false.
   integer                   :: max_i,max_j
!EOP
!-------------------------------------------------------------------------
!  first find the lower left (im,jm) in the m-grid which coresponds to
!  olon(i,j), olat(i,j)
   do j=jl,jh
      do i=il,ih
#ifdef USE_VALID_LON_LAT_ONLY
         if(olon(i,j) .gt. -1000. .and. olat(i,j) .gt. -1000.) then
#else
         if(mask(i,j) .ge. 1) then
#endif
            do jm=1,size(met_lat,2)
               if(met_lat(1,jm) .gt. olat(i,j)) EXIT
            end do
            if (jm .gt. size(met_lat,2)) then
               outside = .true.
            end if
            do im=1,size(met_lon,1)
               if(met_lon(im,1) .gt. olon(i,j)) EXIT
            end do
            if (im .gt. size(met_lon,1)) then
               outside = .true.
            end if
            gridmap(i,j,1) = im-1
            gridmap(i,j,2) = jm-1
         end if
      end do
   end do

   if(outside) then
      STDERR 'WARNING - interpol_coefficients: Some points out side the area'
   end if

!  then calculated the t and u coefficients - via distances - the point of
!  interest is (x,y)
max_i = size(met_lon,1)
max_j = size(met_lat,2) 
   do j=jl,jh
      do i=il,ih
#ifdef USE_VALID_LON_LAT_ONLY
         if(olon(i,j) .gt. -1000. .and. olat(i,j) .gt. -1000.) then
#else
         if(mask(i,j) .ge. 1) then
#endif
            x = olon(i,j)
            y = olat(i,j)
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
                     lon1 = met_lon(im,jm)
                     lat1 = met_lat(im,jm)
                     lon2 = met_lon(im+1,jm+1)
                     lat2 = met_lat(im+1,jm+1)
                     t(i,j) = spherical_dist(earth_radius,x,lat1,lon1,lat1)/ &
                              spherical_dist(earth_radius,lon1,lat1,lon2,lat1)
                     u(i,j) = spherical_dist(earth_radius,lon1,y,lon1,lat1)/ &
                              spherical_dist(earth_radius,lon1,lat1,lon1,lat2)
                  case default
               end select
            end if
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
