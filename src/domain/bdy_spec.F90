#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: bdy_spec() - defines open boundaries
!
! !INTERFACE:
   subroutine bdy_spec(fn)
!
! !DESCRIPTION:
!  Read in the open boundary information from 'fn'.
!
! !USES:
   use exceptions
   use domain, only: NWB,NNB,NEB,NSB,NOB
   use domain, only: wi,wfj,wlj,nj,nfi,nli,ei,efj,elj,sj,sfi,sli
   use domain, only: bdy_index,bdy_map,nsbv
   use domain, only: bdy_2d_type,bdy_3d_type
   use domain, only: need_2d_bdy_elev,need_2d_bdy_u,need_2d_bdy_v
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   character(len=255)        :: line
   integer                   :: iostat
   integer                   :: i,j,k,l
   integer                   :: n,rc
   integer                   :: type_2d(4,10),type_3d(4,10)
!EOP
!-----------------------------------------------------------------------
!BOC

#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'bdy_spec() # ',Ncall
#endif

   LEVEL2 'Reading boundary location information from:'
   LEVEL3 trim(fn)


!   open(BDYINFO,file=fn,status='unknown',ERR=90)

!   open(unit,file=fn,action='read',iostat=iostat,status='old',err=90)
   open(BDYINFO,file=fn,action='read',iostat=iostat,status='old')

!  Western boundary info
   do while (NWB .eq. -1 .and. iostat == 0)
      read(BDYINFO,'(A)',iostat=iostat,end=91,err=92) line
!     skip comments and empty lines
      if (line(1:1) == '#' .or. line(1:1) == '!' .or. &
          len(trim(line)) == 0 ) then
      else
         read(line,*) NWB
      end if
   end do
   if (NWB .lt. 0) then
      call getm_error("bdy_spec()","NWB < 0")
   end if

   if (NWB .ge. 1) then
      allocate(wi(NWB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (wi)'
      allocate(wfj(NWB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (wfj)'
      allocate(wlj(NWB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (wlj)'

      n=0
      do while (n .ne. NWB .and. iostat == 0)
         read(BDYINFO,'(A)',iostat=iostat,end=91,err=92) line
!        skip comments and empty lines
         if (line(1:1) == '#' .or. line(1:1) == '!' .or. &
             len(trim(line)) == 0 ) then
         else
            n=n+1
            read(line,*,END=91,ERR=92) wi(n),wfj(n),wlj(n), &
                                       type_2d(1,n),type_3d(1,n)
            nsbv = nsbv + (wlj(n)-wfj(n)+1)
         end if
      end do
      do n = 1,NWB
         if (type_2d(1,n) .eq. CLAMPED) need_2d_bdy_elev = .true.
         if (type_2d(1,n) .eq. FLATHER_ELEV) then
            need_2d_bdy_elev = .true.
            need_2d_bdy_u    = .true.
         end if
      end do
   end if

!  Northern boundary info
   do while (NNB .eq. -1 .and. iostat == 0)
      read(BDYINFO,'(A)',iostat=iostat,end=91,err=92) line
!     skip comments and empty lines
      if (line(1:1) == '#' .or. line(1:1) == '!' .or. &
          len(trim(line)) == 0 ) then
      else
         read(line,*) NNB
      end if
   end do
   if (NNB .lt. 0) then
      call getm_error("bdy_spec()","NNB < 0")
   end if

   if (NNB .ge. 1 ) then
      allocate(nj(NNB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (nj)'
      allocate(nfi(NNB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (nfi)'
      allocate(nli(NNB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (nli)'

      n=0
      do while (n .ne. NNB .and. iostat == 0)
         read(BDYINFO,'(A)',iostat=iostat,end=91,err=92) line
!        skip comments and empty lines
         if (line(1:1) == '#' .or. line(1:1) == '!' .or. &
             len(trim(line)) == 0 ) then
         else
            n=n+1
            read(line,*,END=91,ERR=92) nj(n),nfi(n),nli(n), &
                                       type_2d(2,n),type_3d(2,n)
            nsbv = nsbv + (nli(n)-nfi(n)+1)
         end if
      end do
      do n = 1,NNB
         if (type_2d(2,n) .eq. CLAMPED) need_2d_bdy_elev = .true.
         if (type_2d(2,n) .eq. FLATHER_ELEV) then
            need_2d_bdy_elev = .true.
            need_2d_bdy_v    = .true.
         end if
      end do
   end if

!  Eastern boundary info
   do while (NEB .eq. -1 .and. iostat == 0)
      read(BDYINFO,'(A)',iostat=iostat,end=91,err=92) line
!     skip comments and empty lines
      if (line(1:1) == '#' .or. line(1:1) == '!' .or. &
          len(trim(line)) == 0 ) then
      else
         read(line,*) NEB
      end if
   end do
   if (NEB .lt. 0) then
      call getm_error("bdy_spec()","NEB < 0")
   end if
   if (NEB .ge. 1 ) then
      allocate(ei(NEB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (ei)'
      allocate(efj(NEB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (efj)'
      allocate(elj(NEB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (elj)'
      n=0
      do while (n .ne. NEB .and. iostat == 0)
         read(BDYINFO,'(A)',iostat=iostat,end=91,err=92) line
!        skip comments and empty lines
         if (line(1:1) == '#' .or. line(1:1) == '!' .or. &
             len(trim(line)) == 0 ) then
         else
            n=n+1
            read(line,*,END=91,ERR=92) ei(n),efj(n),elj(n), &
                                       type_2d(3,n),type_3d(3,n)
            nsbv = nsbv + (elj(n)-efj(n)+1)
         end if
      end do
      do n = 1,NEB
         if (type_2d(3,n) .eq. CLAMPED) need_2d_bdy_elev = .true.
         if (type_2d(3,n) .eq. FLATHER_ELEV) then
            need_2d_bdy_elev = .true.
            need_2d_bdy_u    = .true.
         end if
      end do
   end if

!  Southern boundary info
   do while (NSB .eq. -1 .and. iostat == 0)
      read(BDYINFO,'(A)',iostat=iostat,end=91,err=92) line
!     skip comments and empty lines
      if (line(1:1) == '#' .or. line(1:1) == '!' .or. &
          len(trim(line)) == 0 ) then
      else
         read(line,*) NSB
      end if
   end do
   if (NSB .lt. 0) then
      call getm_error("bdy_spec()","NSB < 0")
   end if
   if (NSB .ge. 1 ) then
      allocate(sj(NSB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (sj)'
      allocate(sfi(NSB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (sfi)'
      allocate(sli(NSB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (sli)'
      n=0
      do while (n .ne. NSB .and. iostat == 0)
         read(BDYINFO,'(A)',iostat=iostat,end=91,err=92) line
!        skip comments and empty lines
         if (line(1:1) == '#' .or. line(1:1) == '!' .or. &
             len(trim(line)) == 0 ) then
         else
            n=n+1
            read(line,*,END=91,ERR=92) sj(n),sfi(n),sli(n), &
                                       type_2d(4,n),type_3d(4,n)
            nsbv = nsbv + (sli(n)-sfi(n)+1)
         end if
      end do
      do n=1,NSB
         if (type_2d(4,n) .eq. CLAMPED) need_2d_bdy_elev = .true.
         if (type_2d(4,n) .eq. FLATHER_ELEV) then
            need_2d_bdy_elev = .true.
            need_2d_bdy_v    = .true.
         end if
      end do
   end if

   close(BDYINFO)

   NOB = NWB+NNB+NEB+NSB

   allocate(bdy_2d_type(NOB),stat=rc)
   if (rc /= 0) stop 'bdy_spec: Error allocating memory (bdy_2d_type)'

   allocate(bdy_3d_type(NOB),stat=rc)
   if (rc /= 0) stop 'bdy_spec: Error allocating memory (bdy_3d_type)'

   l=1
   do n = 1,NWB
      bdy_2d_type(l) = type_2d(1,n)
      bdy_3d_type(l) = type_3d(1,n)
      l=l+1
   end do
   do n = 1,NNB
      bdy_2d_type(l) = type_2d(2,n)
      bdy_3d_type(l) = type_3d(2,n)
      l=l+1
   end do
   do n = 1,NEB
      bdy_2d_type(l) = type_2d(3,n)
      bdy_3d_type(l) = type_3d(3,n)
      l=l+1
   end do
   do n = 1,NSB
      bdy_2d_type(l) = type_2d(4,n)
      bdy_3d_type(l) = type_3d(4,n)
      l=l+1
   end do

   LEVEL2 '# of open boundaries ',NOB
   LEVEL2 '# of open bdy points ',nsbv

   if(NOB .gt. 0) then
      allocate(bdy_index(NOB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (bdy_index)'

      allocate(bdy_map(nsbv,2),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (bdy_map)'

      k = 1
      l = 1
      do n=1,NWB
         bdy_index(l) = k
         l = l+1
         i = wi(n)
         do j=wfj(n),wlj(n)
            bdy_map(k,1) = i
            bdy_map(k,2) = j
            k = k+1
         end do
      end do
      do n=1,NNB
         bdy_index(l) = k
         l = l+1
         j = nj(n)
         do i=nfi(n),nli(n)
            bdy_map(k,1) = i
            bdy_map(k,2) = j
            k = k+1
         end do
      end do
      do n=1,NEB
         bdy_index(l) = k
         l = l+1
         i = ei(n)
         do j=efj(n),elj(n)
            bdy_map(k,1) = i
            bdy_map(k,2) = j
            k = k+1
         end do
      end do
      do n=1,NSB
         bdy_index(l) = k
         l = l+1
         j = sj(n)
         do i=sfi(n),sli(n)
            bdy_map(k,1) = i
            bdy_map(k,2) = j
            k = k+1
         end do
      end do
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving bdy_spec()'
   write(debug,*)
#endif

   return
90 FATAL 'can not open ',TRIM(fn)
   stop
91 STDERR 'EOF ',TRIM(fn)
92 STDERR 'Error reading ',TRIM(fn)

   return
   end subroutine bdy_spec
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
