#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: read_par_setup - test input file integrity
!
! !INTERFACE:
   subroutine read_par_setup(fn,nprocs,myid,imax,jmax,iextr,jextr,  &
                             ioff,joff,neighbours,numthreads)
!
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn ! file to read from
   integer, intent(in)                 :: nprocs,myid  ! Number of jobs actually started.
   integer, intent(in)                 :: imax,jmax,iextr,jextr
!
! !OUTPUT PARAMETERS:
   integer, intent(out)                :: ioff,joff,neighbours(8),numthreads
!
! !DESCRIPTION:
!  Test the content of a file with neighbour list information.
!  At the time of calling, the file must be opened and rewound.  
!  The routine will then rewind the file after use.
!
!  The global grid has extent 1:iextr,1:jextr
!  The local grid has extent imin:imax,jmin:jmax,
!  where it is assumed that imin=jmin=1.
!
! !REVISION HISTORY:
!  2002-02-12 Bjarne Buchmann (bjb@fomfrv.dk) Initial code
!  2009-10-22 Bjarne Buchmann (bjb@frv.dk) Add read number threads per subdomain
!
! !LOCAL VARIABLES:
   integer                   :: myid_read
   integer                   :: nprocs_read,err,ijob,ioff_read,joff_read
   integer                   :: imax_read,jmax_read,iextr_read,jextr_read
   integer                   :: iline, njob, nnjob, ineigh, nthreads_read
   integer                   :: thislineok
   integer, allocatable      :: neighbourlist(:,:), nthreads(:)
   integer                   :: neighbour_inverse(8)  = &
                                 (/5, 6, 7, 8, 1, 2, 3, 4/)
   integer, parameter        :: false_flag = -2 ! No-good PID
   integer, parameter        :: iunit=87 ! check 87 KBK
   character(len=255)        :: line
   integer                   :: iostat,iostat2
!
!EOP
!-------------------------------------------------------------------------
!BOC
! Read #jobs, test vs. actual nprocs in use

   open(unit=iunit,file=fn,iostat=iostat)
   iline   = 0 ! Index for line number in file
   thislineok = 0 ! Flag for this line read OK
   do while (thislineok .eq. 0 .and. iostat == 0)
      iline=iline+1
      read(iunit,'(A)',iostat=iostat,end=1020,err=1010) line
!        skip comments and empty lines
      if (line(1:1) == '#' .or. line(1:1) == '!' .or. &
           len(trim(line)) == 0 ) then
      else
         read(line,*,err=1010,end=1020) nprocs_read
         thislineok=1
      end if
   end do

   if (nprocs_read /= nprocs) then
      FATAL 'read_par_setup: Number of jobs do not match'
      FATAL '  Expected value ',nprocs
      FATAL '  Read value     ',nprocs_read
     stop 
   end if

   allocate(neighbourlist(0:nprocs-1,8),stat=err)
   if (err /= 0) &
      stop 'read_par_setup: Error allocating memory (neighbourlist)'
   allocate(nthreads(0:nprocs-1),stat=err)
   if (err /= 0) &
      stop 'read_par_setup: Error allocating memory (nthreads)'
!
! Flag all neighbourlists to be "unrecongnized value"
   neighbourlist(:,:) =  false_flag
   nthreads(:)        =  -1
!
   thislineok = 0
   do while (thislineok .eq. 0 .and. iostat == 0)
      iline=iline+1
      read(iunit,'(A)',iostat=iostat,end=1020,err=1010) line
      if (line(1:1) == '#' .or. line(1:1) == '!' .or. &
           len(trim(line)) == 0 ) then
      else
         read(line,*,err=1010,end=1020) &
              imax_read,jmax_read,iextr_read,jextr_read
         thislineok=1
      end if
   end do

   if (iextr_read /= iextr .OR. jextr_read /= jextr) then
      FATAL 'read_par_setup: Global grid sizes do not match'
      FATAL '  Expected ',iextr,' by ',jextr
      FATAL '  Read     ',iextr_read,' by ',jextr_read
      stop
   end if 

   if (imax_read /= imax .OR. jmax_read /= jmax) then
      FATAL 'read_par_setup: Local grid sizes do not match'
      FATAL '  Expected ',imax,' by ',jmax
      FATAL '  Read     ',imax_read,' by ',jmax_read
      stop
   end if
!
! Read following lines (one per job)
   do ijob=0,nprocs-1
      thislineok = 0
      do while (thislineok .eq. 0 .and. iostat == 0)
         iline=iline+1
         read(iunit,'(A)',iostat=iostat,end=1020,err=1010) line
         if (line(1:1) == '#' .or. line(1:1) == '!' .or. &
              len(trim(line)) == 0 ) then
         else
! Read lines of format:
!  ID IOFF JOFF 8xNEIGHs numThr
! Accept older format without "numThr" (number of threads)
!   First try new format:
            read(line,*,iostat=iostat2) &
                 myid_read,ioff_read,joff_read,neighbours,nthreads_read
            if (iostat2 .ne. 0) then
!   Try alternative old format
               read(line,*,err=1010,end=1020) &
                    myid_read,ioff_read,joff_read,neighbours
               nthreads_read=0
            end if
            thislineok=1
         end if
      end do

      if(myid_read .eq. myid) then
         ioff = ioff_read
         joff = joff_read
         numthreads = nthreads_read
      end if
!
! Perform straight-forward tests on this input:

      if ( ((ioff_read+imax) .LT. 1) .OR. ((joff_read+jmax) .LT. 1) .OR. &
          (ioff_read .GT. iextr).OR. (joff_read.GT. jextr)    ) then
         FATAL 'read_par_setup: Line ',iline
         FATAL '   Local grid fully outside global grid'
         stop
      end if

      if (myid_read < 0 .OR. minval(neighbours)<-1) then
         FATAL ' read_par_setup: Negative job ID on line ',iline
         stop 
      end if

      if (myid_read > nprocs-1 .OR. maxval(neighbours)> nprocs-1 ) then
         FATAL 'read_par_setup: Job ID appears too large on line ',iline
         stop 
      end if

      if (neighbourlist(myid_read,1) .NE.  false_flag) then
         FATAL 'read_par_setup: Line ',iline
         FATAL '  Job ID ',myid_read,' included twice '
         stop 
      end if
!
!  Store in local array
      neighbourlist(myid_read,1:8) = neighbours(1:8)
   end do
!
!  Test whether the entire neighbour array is coherent...
   do ijob = 0,nprocs-1
      if (minval(neighbourlist(ijob,:))<-1) then
         FATAL ' read_par_setup: Job ID ',ijob, &
               'does not seem to be specified'
         stop 
      end if
!
!  Test each neighbour for consistency on the basic notion that 
!  "my right-hand-side neighbour's left-hand-side neighbour 
!   should be myself". The vector neighbour_inverse contains the 
!  "right-hand-side to left-hand-side" relations.
!  The "inverse lookup list" functions so that my neighbour(i) should
!  have me as neighbour neighbour_inverse(i)
!  Neighbours are indexed clock-wise starting with WEST.
!  This makes the vector have the form 
!      (/5, 6, 7, 8, 1, 2, 3, 4/)
!
      do ineigh=1,8
         njob  = neighbourlist(ijob,ineigh)
         if (njob .gt. -1) then
            nnjob = neighbourlist(njob,neighbour_inverse(ineigh))
            if ( nnjob .ne. ijob ) then
               FATAL 'read_par_setup: Neighbourlist inconsistency'
               FATAL '  Neighbour no. ',ineigh,'of job ',ijob,' is ',njob
               FATAL '  BUT Neighbour ',neighbour_inverse(ineigh),'of job ', &
                     njob,' is NOT ',ijob
               FATAL '  (rather, it is ',nnjob,')'
               stop 
            end if
         end if
      end do
   end do

   LEVEL1 'read_par_setup: Consistency OK'
   close(iunit)

   neighbours(1:8) = neighbourlist(myid,1:8)
   return
!
! Capture errors: 

 1010 continue
   FATAL 'read_par_setup: Unexpected format at line ',iline
   stop

 1020 continue
   FATAL 'read_par_setup: Premature EOF at line',iline
   stop

end subroutine read_par_setup
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2002 - Bjarne Buchmann                                 !
!-----------------------------------------------------------------------
