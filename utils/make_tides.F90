program make_tides

   use time, only: String2JulSecs
   implicit none

   include  'netcdf.inc'
   integer :: yy1=-1,mm1=-1,dd1=-1
   integer :: yy2=-1,mm2=-1,dd2=-1
   character(len=128) :: fn
   character(len=19) :: t1,t2
   real(4), dimension(:), allocatable :: lons,lats
   real(4), dimension(:), allocatable :: z
   real(8) :: hh=0.d0
   character(len=20) :: model='legi'
   integer :: i,j,n,nmax
   integer :: rc
   integer :: ncid,time_id,jul_id,secs_id,elev_id
   integer :: start(2), edges(2)
   integer :: jul0,secs0
   integer :: jul1,secs1
   integer :: jul2,secs2
   integer :: jd,ss,ia(1)

   integer :: jnasa=13222 ! jnasa=0 --> 1958/01/01
   character(len=19) :: t0 = '1958/01/01 00:00:00'

   call String2JulSecs(t0,jul0,secs0)
!
   write(0,*) 'Reading time info from tides.dates'

   open(unit=10,file="tides.dates",status='old',err=120)
   read(10,*) yy1,mm1,dd1
   read(10,*) yy2,mm2,dd2
   close(10)

   write(t1,'(i4.4,a1,i2.2,a1,i2.2xa8)')  yy1,'-',mm1,'-',dd1,'00:00:00'
   write(t2,'(i4.4,a1,i2.2,a1,i2.2xa8)')  yy2,'-',mm2,'-',dd2,'00:00:00'

   call String2JulSecs(t1,jul1,secs1)
   call String2JulSecs(t2,jul2,secs2)


   jnasa = jul1 - jul0
   nmax = (jul2-jul1)*24 + 1
   write(0,*) t0,jul0,secs0
   write(0,*) t1,jul1,secs1
   write(0,*) t2,jul2,secs2
   write(0,*) jnasa,nmax
   jd = jul1
   ss = secs1

   write(0,*) 'Start: ',t1
   write(0,*) 'Stop:  ',t2

   read(5,*) n
   allocate(lons(n),stat=rc)
   if (rc /= 0) stop 'make_tides: Error allocating memory (lons)'
   allocate(lats(n),stat=rc)
   if (rc /= 0) stop 'make_tides: Error allocating memory (lats)'
   allocate(z(n),stat=rc)
   if (rc /= 0) stop 'make_tides: Error allocating memory (z)'
   do i=1,n
      read(5,*) lats(i),lons(i)
   end do

#undef NETCDF
#define NETCDF
#ifdef NETCDF
   fn = 'bdy_2d.nc'

   write(fn,100)  'tides.',yy1,'-',mm1,'-',dd1,'_',yy2,'-',mm2,'-',dd2,'.nc'

100 format (a6,i4.4,a1,i2.2,a1,i2.2,a1,i4.4,a1,i2.2,a1,i2.2,a3)

   call create_nc_file(fn,n,t1,ncid,time_id,jul_id,secs_id,elev_id)

   do j=1,nmax
      z(1) = (j-1)*3600.
      start(1) = j
      edges(1) = 1
      rc = nf_put_vara_real(ncid, time_id, start, edges, z(1))
      ia(1) = jd
      rc = nf_put_vara_int(ncid, jul_id, start, edges, ia)
      ia(1) = ss
      rc = nf_put_vara_int(ncid, secs_id, start, edges, ia)
      ss = ss+3600
      do i=1,n
         call otide(model,z(i),lats(i),lons(i),jnasa,hh,rc)
	 if (z(i) .lt. 10000.) then
            z(i) =  0.01*z(i)
	 else
            z(i) = 0.00
	 end if
      end do
      start(1) = 1
      start(2) = j
      edges(1) = n
      edges(2) = 1
      rc = nf_put_vara_real(ncid, elev_id, start, edges, z)

      hh=hh+1.d0
      if( hh .gt. 23 ) then
         jnasa = jnasa + 1
	 hh = 0.d0
	 jd = jd+1
	 ss = 0
      end if
   end do
   rc = nf_close(ncid)
#else
   do j=1,nmax
      print *, (j-1)*3600.,' ,'
   end do
   do j=1,nmax
      do i=1,n
         call otide(model,z(i),lats(i),lons(i),jnasa,hh,rc)
	 if (z(i) .lt. 10000.) then
            print *, z(i)/100.,' ,'
	 else
            print *, -99.00,' ,'
	 end if
      end do
!  print *,jnasa, hh,  z(1)
      hh=hh+1.d0
      if( hh .gt. 23 ) then
         jnasa = jnasa + 1
	 hh = 0.d0
      end if
   end do
#endif

   stop

120 write(0,*) 'Unable to open tides.dates'
   stop 'make_tides'
end
