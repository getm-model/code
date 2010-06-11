!$Id: time.F90,v 1.5 2007-08-24 10:43:45 frv-bjb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  time - keeps control of time
!
! !INTERFACE:
   MODULE time
!
! !DESCRIPTION:
!  This module implements a set of variables and subroutines which should
!  make it relative easy to keep track of time - both model time and real
!  time. The main principle is that times are converted to true Julian Day
!  and seconds of that day - starting at midnight. Both are stored as
!  integers. The main advantages of this approach are 1) very easy to
!  subtract and add time, 2) no problems with round off errors as if
!  real values have been used.
!  The main disadvantge is that 2 variables are used to keep track of time.
!  A number of subroutines/functions are supplied for doing specific things.
!  The module provides the 3 integer variables $julianday, secondsofday$ and
!  $yearday$ which contains the time. These variables should only be updated
!  via a call to \emph{update\_time()}. Any model component which needs
!  information about time should just have a
!  \emph{use time, only: julianday, secondsofday} to get the actual model time.
!
! !USES:
   IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
   public                              :: time_diff

! !PUBLIC DATA MEMBERS:
   integer                             :: julianday,secondsofday,yearday,month,day
   integer                             :: jul0=-1,secs0=-1
   integer                             :: juln=-1,secsn=-1
   REALTYPE                            :: fsecs,simtime
   REALTYPE                            :: timestep
   character(len=19)                   :: timestr
   character(len=19)                   :: start='2000-01-01 00:00:00',stop
   integer                             :: leapyear=0,days_in_mon(0:1,12)
   integer, parameter                  :: secsprday=86400
!
! !private DATA MEMBERS:
   logical, private          :: HasRealTime=.true.
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: time.F90,v $
!  Revision 1.5  2007-08-24 10:43:45  frv-bjb
!  Allow negative seconds in meteo nc-files input
!
!  Revision 1.4  2004-04-06 16:32:28  kbk
!  TimeDiff --> time_diff
!
!  Revision 1.3  2003/12/15 16:03:59  kbk
!  new correct algorithm in in_interval()
!
!  Revision 1.2  2003/04/23 12:02:43  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.1.1.1  2002/05/02 14:01:19  gotm
!  recovering after CVS crash
!
!  Revision 1.8  2001/10/17 08:19:42  bbh
!  jul0 and jul1 now public
!
!  Revision 1.7  2001/09/28 12:29:22  bbh
!  Added add_secs() and in_interval()
!
!  Revision 1.6  2001/07/26 13:52:41  bbh
!  Cleaning + description
!
!  Revision 1.5  2001/05/14 12:40:23  bbh
!  Removed test output in string_to_time()
!
!  Revision 1.4  2001/05/07 19:52:56  bbh
!  Included - as valid delimiter in string_to_julsecs
!
!  Revision 1.3  2001/05/07 14:37:45  bbh
!  Added string_to_julsecs() - reads time imbedded in a string
!
!  Revision 1.2  2001/05/06 18:51:55  bbh
!  Towards proper implementation of specified 2D bdy.
!
!  Revision 1.1.1.1  2001/04/17 08:43:09  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_time - initialise the time system in getm
!
! !INTERFACE:
   subroutine init_time(MinN,MaxN)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   integer, intent(out)                :: MinN,MaxN
!
! !DESCRIPTION:
!  Reads the namelist and makes calls to the init functions of the
!  various model components.
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
   integer                   :: timefmt=1
   integer                   :: jul1,secs1,jul2,secs2
   integer                   :: ndays,nsecs
   integer                   :: nfirst,nlast
   namelist /time/ timestep,timefmt,nlast,start,stop
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_time() # ',Ncall
#endif
   days_in_mon(0,:) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
   days_in_mon(1,:) = (/31,29,31,30,31,30,31,31,30,31,30,31/)
!
!  Read time specific things from the namelist.
!
   LEVEL1 'init_time'
   READ(NAMLST,NML=time)
!
!  Calculate MaxN -> MinN is 1 if not changed by HotStart
!
   MinN = 1
   MaxN = nlast
   LEVEL2 'Time step:      ',timestep,' seconds'
   LEVEL2 'Time format:    ',timefmt
   select case (timefmt)
      case (0)
!KBK
         LEVEL2 'Hopefully we will get the time from the hot start file'
      case (1)
         HasRealTime=.false.
         LEVEL2 '# of timesteps: ',MaxN
         start='2000-01-01 00:00:00'

         call String2JulSecs(start,jul1,secs1)

         nsecs = nint(MaxN*timestep) + secs1
         ndays = nsecs/86400
         jul2  = jul1 + ndays
         secs2 = mod(nsecs,86400)
         call write_time_string(jul2,secs2,stop)

         LEVEL2 'Fake start:     ',start
         LEVEL2 'Fake stop:      ',stop
      case (2)
         LEVEL2 'Start:          ',start
         LEVEL2 'Stop:           ',stop

         call String2JulSecs(start,jul1,secs1)
         call String2JulSecs(stop,jul2,secs2)

         nsecs = time_diff(jul2,secs2,jul1,secs1)
         MaxN  = nint(nsecs/timestep)

         ndays = jul2-jul1
         if (nsecs .lt. 86400 .and. jul1 .ne. jul2) ndays = ndays-1
         nsecs = nsecs - 86400*ndays
         STDERR '        ==> ',ndays,' day(s) and ',nsecs,' seconds ==> ',MaxN,' micro time steps'
      case (3)
         LEVEL2 'Start:          ',start
         LEVEL2 '# of timesteps: ',MaxN

         call String2JulSecs(start,jul1,secs1)

         nsecs = nint(MaxN*timestep) + secs1
         ndays = nsecs/86400
         jul2  = jul1 + ndays
         secs2 = mod(nsecs,86400)

         call write_time_string(jul2,secs2,stop)
         LEVEL2 'Stop:           ',stop
      case default
         STDERR 'Fatal error: A non valid input format has been chosen'
         stop 'init_time'
   end select

   jul0  = jul1
   secs0 = secs1

   juln  = jul2
   secsn = secs2

   julianday    = jul0
   secondsofday = secs0

   simtime = timestep*(MaxN-MinN+1)

   call update_time(0)

#ifdef DEBUG
   write(debug,*) 'Leaving init_time()'
   write(debug,*)
#endif
   return
   end subroutine init_time
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CalDat() - converts true Julian day to calender date
!
! !INTERFACE:
   subroutine CalDat(julian,yyyy,mm,dd)
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   integer  julian
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   integer  yyyy,mm,dd
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
   integer, parameter        :: IGREG=2299161
   integer                   :: ja,jb,jc,jd,je
   REAL                      :: x
!EOP
!-----------------------------------------------------------------------
!BOC
   if(julian .ge. IGREG ) then
      x = ((julian-1867216)-0.25)/36524.25
      ja = julian+1+int(x)-int(0.25*x)
   else
      ja = julian
   end if

   jb = ja+1524
   jc = int(6680 + ((jb-2439870)-122.1)/365.25)
   jd = int(365*jc+(0.25*jc))
   je = int((jb-jd)/30.6001)

   dd = jb-jd-int(30.6001*je)
   mm = je-1
   if (mm .gt. 12) mm = mm-12
   yyyy = jc - 4715
   if (mm .gt. 2) yyyy = yyyy-1
   if (yyyy .le. 0) yyyy = yyyy-1

   return
   end subroutine CalDat

!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  JulDay() - converts a calendar date to true Julian day
!
! !INTERFACE:
   subroutine JulDay(yyyy,mm,dd,julian)
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   integer  yyyy,mm,dd
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   integer  julian
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
   integer, parameter        :: IGREG=15+31*(10+12*1582)
   integer                   :: ja,jy,jm
!EOP
!-----------------------------------------------------------------------
!BOC
!
   jy = yyyy
   if(jy .lt. 0) jy = jy+1
   if (mm .gt. 2) then
      jm = mm+1
   else
      jy = jy-1
      jm = mm+13
   end if
   julian = int(floor(365.25*jy)+floor(30.6001*jm)+dd+1720995)
   if (dd+31*(mm+12*yyyy) .ge. IGREG) then
      ja = int(0.01*jy)
      julian = julian+2-ja+int(0.25*ja)
   end if

   return
   end subroutine JulDay

!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  update_time() - keep track of real time (Julian Days and secs)
!
! !INTERFACE:
   subroutine update_time(n)
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: n
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
   integer                   :: nsecs
   integer                   :: yy,mm,dd,jd
!EOP
!-----------------------------------------------------------------------
!BOC

   nsecs = nint(n*timestep) + secs0
   fsecs = n*timestep + secs0
   julianday    = jul0 + nsecs/86400
   secondsofday = mod(nsecs,86400)

   call CalDat(julianday,yy,mm,dd)
   month=mm
   day=dd
   call JulDay(yy,1,1,jd)
   yearday = julianday-jd+1

   leapyear=0
   if((mod(yy,4) .eq. 0 .and. mod(yy,100) .ne. 0) .or. mod(yy,400) .eq. 0) then
      leapyear=1
   end if

   return
   end subroutine update_time
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  string_to_julsecs() - converts time string to julian day and secs
!
! !INTERFACE:
   subroutine string_to_julsecs(str,jul,secs)
!
! !DESCRIPTION:
!  Searches for [:/-] in string - subtracts 4 (yyyy) and reads the time
!  yyyy/mm/dd hh:dd:ss. Fails if the first occurence of [:/] is not in the
!  date field - but before.
!
! !USES:
!
! !INPUT PARAMETERS:
   character(len=*)                    :: str
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   integer, intent(out)                :: jul,secs
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
   character                 :: c1,c2,c3,c4
   integer                   :: yy,mm,dd,hh,min,ss
   integer                   :: i
   character(len=3)          :: set='/:-'
!EOP
!-----------------------------------------------------------------------
!BOC
   i = scan(str,set)
   if(i .lt. 5) then
      FATAL "Can't read valid time in: ", trim(str)
      stop 'string_to_julsecs'
   end if
   read(str(i-4:),'(i4,a1,i2,a1,i2,1x,i2,a1,i2,a1,i2)')  &
                          yy,c1,mm,c2,dd,hh,c3,min,c4,ss
   call JulDay(yy,mm,dd,jul)
   secs = 3600*hh + 60*min + ss
#ifdef DEBUG_TIME
   STDERR 'string_to_julsecs" ',yy,mm,dd,hh,min,ss
#endif
   return
   end subroutine string_to_julsecs
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  String2JulSecs() - converts time string to julian day and secs
!
! !INTERFACE:
   subroutine String2JulSecs(timestr,jul,secs)
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   character(len=19)                   :: timestr
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   integer, intent(out)                :: jul,secs
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
   character                 :: c1,c2,c3,c4
   integer                   :: yy,mm,dd,hh,min,ss
!EOP
!-----------------------------------------------------------------------
!BOC

   READ(timestr,'(i4,a1,i2,a1,i2,1x,i2,a1,i2,a1,i2)')  &
                          yy,c1,mm,c2,dd,hh,c3,min,c4,ss
   call JulDay(yy,mm,dd,jul)
   secs = 3600*hh + 60*min + ss

   return
   end subroutine String2JulSecs

!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  write_time_string() - formats time (julian day and secs) to a time string
!
! !INTERFACE:
   subroutine write_time_string(jul,secs,str)
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   integer, intent(in), optional       :: jul,secs
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   character(len=19), optional         :: str
!
! !REVISION HISTORY:
!  See module for log.
!
! !LOCAL VARIABLES:
   integer                   :: ss,min,hh,dd,mm,yy
!EOP
!-----------------------------------------------------------------------
!BOC

   if (present(secs)) then
      hh   = secs/3600
      min  = (secs-hh*3600)/60
      ss   = secs - 3600*hh - 60*min
   else
      hh   = secondsofday/3600
      min  = (secondsofday-hh*3600)/60
      ss   = secondsofday - 3600*hh - 60*min
   end if

   if (present(jul)) then
      call CalDat(jul,yy,mm,dd)
   else
      call CalDat(julianday,yy,mm,dd)
   end if

   if(present(str)) then
      write(str,'(i4.4,a1,i2.2,a1,i2.2,1x,i2.2,a1,i2.2,a1,i2.2)')  &
                        yy,'-',mm,'-',dd,hh,':',min,':',ss
   else
      write(timestr,'(i4.4,a1,i2.2,a1,i2.2,1x,i2.2,a1,i2.2,a1,i2.2)')  &
                        yy,'-',mm,'-',dd,hh,':',min,':',ss
   end if

   return
   end subroutine write_time_string

!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  time_diff() - the time difference in seconds
!
! !INTERFACE:
   integer FUNCTION time_diff(jul1,secs1,jul2,secs2)
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: jul1,secs1,jul2,secs2
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   time_diff = 86400*(jul1-jul2) + (secs1-secs2)
   return
   end function  time_diff
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  add_secs() -
!
! !INTERFACE:
   subroutine add_secs(j1,s1,secs,j2,s2)
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: j1,s1,secs
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   integer, intent(out)                :: j2,s2
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   s2 = s1 + mod(secs,86400)
   j2 = j1 + secs/86400
   if (s2 .gt. 86400) then
      s2 = s2 - 86400
      j2 = j2 +1
   else if (s2 .lt. 0) then
      s2 = s2 + 86400
      j2 = j2 -1
   end if
   return
   end subroutine  add_secs
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  in_interval() -
!
! !INTERFACE:
   logical function in_interval(j1,s1,j,s,j2,s2)
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: j1,s1,j,s,j2,s2
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
   logical         :: before,after
!
!EOP
!-----------------------------------------------------------------------
!BOC

   before = (j .lt. j1) .or. ( j .eq. j1 .and. (s .lt. s1) )
   after  = (j .gt. j2) .or. ( j .eq. j2 .and. (s .gt. s2) )

   in_interval = ( .not. before ) .and. ( .not. after )
   return
   end function in_interval
!EOC

!-----------------------------------------------------------------------

   end module time

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
