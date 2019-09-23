#include"cppdefs.h"
   program test_rotation
!
! !DESCRIPTION:
!  program to test the rotated grid features
!  To execute:
!  make test_rotation
!
! !USES:
   use grid_interpol, only: to_rotated_lat_lon,from_rotated_lat_lon
   use exceptions
   implicit none
! !LOCAL VARIABLES
   integer  :: i,j,n
   REALTYPE :: southpole(3)
   REALTYPE :: lon,lat,rlon,rlat,beta
   REALTYPE :: rlonv(6),rlatv(6)
!EOP
!-----------------------------------------------------------------------
!BOC
! rotated_pole:grid_north_pole_latitude  = 39.25 ;
! rotated_pole:grid_north_pole_longitude = -162. ;

!#define JUJ
#ifdef JUJ
   southpole(1) = 0.  ! lat
   southpole(2) = 80. ! lon
   southpole(3) = 0.

   STDERR 'sp1=  ',southpole(1)
   STDERR 'sp2=  ',southpole(2)
   STDERR 'sp3=  ',southpole(3)

   rlonv(1) = -30.127
   rlatv(1) = -10.691
   rlonv(2) = -20.087
   rlatv(2) = -13.651

   do n=1,2
      STDERR 'rlon= ',rlonv(n)
      STDERR 'rlat= ',rlatv(n)

      call from_rotated_lat_lon(southpole,rlonv(n),rlatv(n),lon,lat,beta)

      STDERR 'lon=  ',lon
      STDERR 'lat=  ',lat
      STDERR 'beta= ',beta
   end do

   lon = 80.
   lat = 0.
   call to_rotated_lat_lon(southpole,lon,lat,rlon,rlat,beta)

   STDERR 'rlon= ',rlon
   STDERR 'rlat= ',rlat
   STDERR 'beta= ',beta
#endif
#define ADOLF
#ifdef ADOLF
   STDERR 'ADOLF'
   southpole(2) = -162.+180. ! lon
   southpole(1) = -39.25     ! lat
   southpole(3) = 0.
   STDERR 'sp1=  ',southpole(1)
   STDERR 'sp2=  ',southpole(2)
   STDERR 'sp3=  ',southpole(3)

   rlonv(1) = -28.375 ; rlatv(1) = -23.375
   rlonv(2) =  18.155 ; rlatv(2) = -23.375
   rlonv(3) =  18.155 ; rlatv(3) =  21.835
   rlonv(4) = -28.375 ; rlatv(4) =  21.835
   rlonv(5) =  -3.625 ; rlatv(5) =   7.645
   rlonv(6) = -19.465 ; rlatv(6) =  -4.235


   do n=5,6
      STDERR 'rlon= ',rlonv(n)
      STDERR 'rlat= ',rlatv(n)

      call from_rotated_lat_lon(southpole,rlonv(n),rlatv(n),lon,lat,beta)

      STDERR 'lon=  ',lon
      STDERR 'lat=  ',lat
      STDERR 'beta= ',beta
   end do

   lon = -162.+180.
   lat = -39.25
   call to_rotated_lat_lon(southpole,lon,lat,rlon,rlat,beta)

   STDERR 'rlon= ',rlon
   STDERR 'rlat= ',rlat
   STDERR 'beta= ',beta
#endif
#define KURT
#ifdef KURT
!  http://se.mathworks.com/matlabcentral/fileexchange/43435-rotated-grid-transform
   STDERR 'KURT'
   southpole(2) = 18.
   southpole(1) = -39.3
   southpole(3) = 0.
   STDERR 'sp1=  ',southpole(1)
   STDERR 'sp2=  ',southpole(2)
   STDERR 'sp3=  ',southpole(3)

   rlonv(1) = 12.0; rlatv(1) = 55.
   rlonv(2) = 12.0; rlatv(2) = 54.
   rlonv(3) = 12.0; rlatv(3) = 53.


   do n=1,3
      STDERR 'rlon= ',rlonv(n)
      STDERR 'rlat= ',rlatv(n)

!      call from_rotated_lat_lon(southpole,rlonv(n),rlatv(n),lon,lat,beta)
      call to_rotated_lat_lon(southpole,rlonv(n),rlatv(n),lon,lat,beta)

      STDERR 'lon=  ',lon
      STDERR 'lat=  ',lat
      STDERR 'beta= ',beta
   end do

   call from_rotated_lat_lon(southpole,lon,lat,rlon,rlat,beta)

   STDERR 'rlon= ',rlon
   STDERR 'rlat= ',rlat
   STDERR 'beta= ',beta
#endif

   end program test_rotation
!EOC
