#include "cppdefs.h"
   program main

   use halo_zones, only: init_halo_zones,wait_halo
   use halo_zones, only: update_2d_halo,update_3d_halo,set_active_communicator
   use halo_zones, only: myid, nprocs, comm_hd
   IMPLICIT NONE

#ifndef SYLT_TEST
   integer, parameter 	:: iextr=175,jextr=158,kmax=5
#endif

#ifdef STATIC
#ifdef SYLT_TEST
!#include "/home/kbk/getm-cases/v1.x/sylt/sylt.dim"
#include "/home/kbk/BBH/frv/setups/dk_06nm/dk_06nm.dim
#else
   integer, parameter	:: imin=1,imax=6,jmin=1,jmax=6
!   integer, parameter	:: imin=1,imax=iextr,jmin=1,jmax=jextr
   integer, parameter	:: iimin=imin,iimax=imax,jjmin=jmin,jjmax=jmax
#endif
   REALTYPE		:: A(E2DFIELD)
   REALTYPE		:: B(I3DFIELD)
   integer		:: mask(E2DFIELD)
#else
#ifdef SYLT_TEST
   integer, parameter	:: iextr=135,jextr=160,kmax=10
#endif
   integer 		:: imin,imax,jmin,jmax
   integer 		:: iimin,iimax,jjmin,jjmax
   REALTYPE, DIMENSION(:,:), ALLOCATABLE :: A
   REALTYPE, DIMENSION(:,:,:), ALLOCATABLE :: B
   integer, DIMENSION(:,:), ALLOCATABLE :: mask
#endif

   integer		:: domain=19
   integer		:: west=20,north=21,east=22,south=23
   character(len=20)	:: fname,fw,fn,fe,fs
   integer  		:: i,j,k,n,ierr

   logical		:: parallel=.true.
!   logical		:: parallel=.false.
   integer		:: id

   if(parallel) then
      call init_halo_zones(parallel,iextr,jextr,imin,imax,jmin,jmax,kmax)
   else
      call init_halo_zones(parallel)
   end if

#ifndef STATIC
   if( .not. parallel) then
      imin = 1; imax = iextr
      jmin = 1; jmax = jextr
   end if
   iimin = imin; iimax = imax
   jjmin = jmin; jjmax = jmax
   allocate(A(E2DFIELD))
   allocate(B(I3DFIELD))
   allocate(mask(E2DFIELD))
   STDERR iextr,jextr,imin,imax,jmin,jmax,kmax
#endif

   mask = 1

   write(*,'(a,1x,5(i4))') 'part_domain:',myid,imin,imax,jmin,jmax

   if (parallel) then
      call set_active_communicator(comm_hd)
   end if

   if(myid .ne. -1) then
      write(fname,'(i2.2,a7)') myid,'-domain'
      write(fw,'(i2.2,a5)') myid,'-west'
      write(fn,'(i2.2,a6)') myid,'-north'
      write(fe,'(i2.2,a5)') myid,'-east'
      write(fs,'(i2.2,a6)') myid,'-south'
      id = myid+1
   else
      write(fname,'(a9)') 'yy-domain'
      write(fw,'(a7)') 'yy-west'
      write(fn,'(a8)') 'yy-north'
      write(fe,'(a7)') 'yy-east'
      write(fs,'(a8)') 'yy-south'
      id = 1
   end if

   open(unit=domain,file=fname,status='replace')
   open(unit=west,file=fw,status='replace')
   open(unit=north,file=fn,status='replace')
   open(unit=east,file=fe,status='replace')
   open(unit=south,file=fs,status='replace')

   A = id*10

   A(imin,jmin+1:jmax-1) = id*10+1 ! west
   A(imin+1:imax-1,jmax) = id*10+2 ! north
   A(imax,jmin+1:jmax-1) = id*10+3 ! east
   A(imin+1:imax-1,jmin) = id*10+4 ! south

   B = id*100

   do k=1,kmax
      B(imin:imax,jmin:jmax,k) = id*100+k
      B(imin,jmin+1:jmax-1,k) = id*100+10+k ! west
      B(imin+1:imax-1,jmax,k) = id*100+20+k ! north
      B(imax,jmin+1:jmax-1,k) = id*100+30+k ! east
      B(imin+1:imax-1,jmin,k) = id*100+40+k ! south
   end do

#if 0
   write(domain,*) 'id = ',myid
   write(domain,*) 'Before'
   do j=jmax+1,jmin-1,-1
      write(domain,*) (INT(A(i,j)), i=imin-1,imax+1)
   end do

   write(west,*)  'west:  myid = ',myid, ' before'
   write(east,*)  'east:  myid = ',myid, ' before'
   write(north,*) 'north: myid = ',myid, ' before'
   write(south,*) 'south: myid = ',myid, ' before'
   do k=kmax,0,-1
      write(west,*)  (INT(B(imin,:,k)))
      write(east,*)  (INT(B(imax,:,k)))
      write(north,*) (INT(B(:,jmax,k)))
      write(south,*) (INT(B(:,jmin,k)))
   end do
#endif

   do n=1,1

!     Communicate - p. 67
      call update_2d_halo(A,A,mask,imin,jmin,imax,jmax,10)
      call wait_halo(10)

      call update_3d_halo(B,B,mask,iimin,jjmin,iimax,jjmax,kmax,30)
      call wait_halo(30)

#if 1
      write(domain,*) 'Loop ',n
      do j=jmax+1,jmin-1,-1
         write(domain,*) (INT(A(i,j)), i=imin-1,imax+1)
      end do

      write(west,*)  'west:  myid = ',myid, '  loop = ',n
      write(east,*)  'east:  myid = ',myid, '  loop = ',n
      write(north,*) 'north: myid = ',myid, '  loop = ',n
      write(south,*) 'south: myid = ',myid, '  loop = ',n
      do k=kmax,0,-1
         write(west,*)  (INT(B(imin-1,:,k)))
         write(east,*)  (INT(B(imax+1,:,k)))
         write(north,*) (INT(B(:,jmax+1,k)))
         write(south,*) (INT(B(:,jmin-1,k)))
      end do
#endif

   end do

   close(unit=domain)
   close(unit=west)
   close(unit=north)
   close(unit=east)
   close(unit=south)

#ifdef STATIC
   if(myid .le. 0) STDERR 'STATIC, PARALLEL = ',parallel
#else
   if(myid .le. 0) STDERR 'DYNAMIC, PARALLEL = ',parallel
#endif

   if(parallel) then
      call MPI_FINALIZE(ierr)
   end if

   end
