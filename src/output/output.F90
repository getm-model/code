!$Id: output.F90,v 1.5 2003-09-16 07:45:30 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  output - output specifications
!
! !INTERFACE:
   module output
!
! !DESCRIPTION:
!
! !USES:
   use time, only: write_time_string,timestep,timestr
   use ncdf_out
   use ascii_out
#ifdef TEST_NESTING
   use nesting
#endif
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   integer                             :: out_fmt=NETCDF
   character(LEN = PATH_MAX)           :: in_dir='.'
   character(LEN = PATH_MAX)           :: out_dir='.'
   character(LEN = PATH_MAX)           :: out_f_2d
   character(LEN = PATH_MAX)           :: out_f_3d
   character(LEN = PATH_MAX)           :: hot_out

   logical                             :: save_meteo=.false.
   logical                             :: save_2d=.true.
   logical                             :: save_3d=.true.
   logical                             :: save_vel=.true.
   logical                             :: save_strho=.true.
   logical                             :: save_s=.true.
   logical                             :: save_t=.true.
   logical                             :: save_rho=.true.
   logical                             :: save_turb=.true.
   logical                             :: save_tke=.true.
   logical                             :: save_eps=.true.
   logical                             :: save_num=.true.
   logical                             :: save_nuh=.true.
   logical                             :: save_spm=.false.
   integer                             :: first_2d=1
   integer                             :: step_2d=1
   integer                             :: first_3d=1
   integer                             :: step_3d=1
   integer                             :: hotout=-1
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: output.F90,v $
!  Revision 1.5  2003-09-16 07:45:30  kbk
!  additional info written when hotstart time mismatch
!
!  Revision 1.4  2003/09/03 05:55:13  kbk
!  continuous=.false. - allows change of dt
!
!  Revision 1.3  2003/04/23 12:07:12  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 12:32:58  kbk
!  parallel support + NO_3D, NO_BAROCLINIC
!
!  Revision 1.1.1.1  2002/05/02 14:01:52  gotm
!  recovering after CVS crash
!
!  Revision 1.8  2001/10/23 07:37:17  bbh
!  Saving spm - if calc_spm and save_spm are both true
!
!  Revision 1.7  2001/10/18 07:14:57  bbh
!  Resolved conflicts
!
!  Revision 1.6  2001/09/18 20:08:11  bbh
!  Added save_nuh
!
!  Revision 1.5  2001/09/13 14:50:02  bbh
!  Cleaner and smaller NetCDF implementation + better axis support
!
!  Revision 1.4  2001/07/26 13:57:14  bbh
!  Meteo working - needs some polishing
!
!  Revision 1.3  2001/06/04 13:09:53  bbh
!  Includes - dryrun - in call to init_output()
!
!  Revision 1.2  2001/04/24 08:24:58  bbh
!  Use runtype instead of macro
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_output - initialise all external files and units
!
! !INTERFACE:
   subroutine init_output(runid,title,starttime,runtype,dryrun,myid)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: runid,title,starttime
   integer, intent(in)                 :: runtype,myid
   logical, intent(in)                 :: dryrun
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  See log for module
!
! !LOCAL VARIABLES:
   namelist /io_spec/ &
             out_fmt, &
             in_dir,out_dir, &
             save_2d,save_3d,save_vel, &
             save_strho,save_s,save_t,save_rho, &
             save_turb,save_tke,save_eps,save_num,save_nuh, &
             save_spm, &
             first_2d,step_2d,first_3d,step_3d,hotout, &
             save_meteo
!   logical :: nesting=.true.
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_output() # ',Ncall
#endif

   LEVEL1 'init_output'

   read(NAMLST, nml=io_spec)

   if (runtype .eq. 1) then
      save_3d = .false.
      save_vel = .false.
      save_strho = .false.
      save_turb = .false.
      save_spm = .false.
   end if

   call file_names(runid,myid)

   if(save_2d) then
      LEVEL2 '2D results: ',trim(out_f_2d)
      LEVEL2 'First=',first_2d,' step=',step_2d
      if(save_meteo) then
         LEVEL2 'Saving meteo forcing in ',trim(out_f_2d)
      end if
   end if
#ifndef NO_3D
   if(save_3d) then
      LEVEL2 '3D results: ',trim(out_f_3d)
      LEVEL2 'First=',first_3d,' step=',step_3d
   end if
#endif

   if( .not. dryrun) then
      select case (out_fmt)
         case (ASCII)
#if 0
            if (save_2d) call init_2d_ascii(out_f_2d,title,starttime)
            if (save_3d) call init_3d_ascii(out_f_3d,title,starttime)
#else
            STDERR 'ASCII output - not coded yet'
            stop 'init_output'
#endif
         case (NETCDF)
            if (save_2d) call init_2d_ncdf(out_f_2d,title,starttime)
#ifndef NO_3D
            if (save_3d) call init_3d_ncdf(out_f_3d,title,starttime)
#endif
         case (GRADS)
         case DEFAULT
           STDERR 'Fatal error: A non valid input format has been chosen'
           stop 'init_output'
      end select
   end if
#ifdef TEST_NESTING
!   if(nesting) then
      call init_nesting()
!   end if
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving init_output()'
   write(debug,*)
#endif
   return
   end subroutine init_output

!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_output - write model results to file(s)
!
! !INTERFACE:
   subroutine do_output(runtype,n,timestep)
!
! !DESCRIPTION:
!  Writes calculated fields to files.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: runtype,n
   REALTYPE, intent(in)                :: timestep
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
   REALTYPE                  :: secs
   logical                   :: write_2d,write_3d
   integer                   :: dummy
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_output() # ',Ncall
#endif

   write_2d = save_2d .and. n .ge. first_2d .and. mod(n,step_2d).eq.0
   write_3d = save_3d .and. n .ge. first_3d .and. mod(n,step_3d).eq.0

   if (write_2d .or. write_3d) then
      call write_time_string()
      LEVEL2 'Saving.... ',timestr
!      call divergence()
      secs = n*timestep
      select case (out_fmt)
         case (ASCII)
#if 0
            if (write_2d) call save_2d_ascii(secs,save_meteo)
            if (write_3d) call save_3d_ascii(secs)
#else
            STDERR 'ASCII output - not coded yet'
            stop 'do_output'
#endif
         case (NETCDF)
            if (write_2d) call save_2d_ncdf(secs)
#ifndef NO_3D
            if (write_3d) call save_3d_ncdf(secs)
#endif
         case DEFAULT
           STDERR 'Fatal error: A non valid input format has been chosen'
           stop 'do_output'
      end select
   end if

!  Restart file
   if (hotout .gt. 0 .and. mod(n,hotout) .eq. 0) then
      dummy = n
      call restart_file(WRITING,trim(hot_out),dummy,runtype)
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving do_output()'
   write(debug,*)
#endif
   return
   end subroutine do_output

!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: restart_file - read/write the restart file
!
! !INTERFACE:
   subroutine restart_file(mode,fname,loop,runtype)
!
! !DESCRIPTION:
!  This routine write the variables necesaary for a 'hot' model start.
!
! !USES:
   use time, only: timestep,julianday,secondsofday
   use m2d,  only: z,zo,U,fU,zu,zub,SlUx,Slru,V,fV,zv,zvb,SlVx,Slrv,Uint,Vint
#ifndef NO_3D
   use m3d,  only: uu,vv,tke,eps,num,nuh,ssen,ssun,ssvn
#ifndef NO_BAROCLINIC
   use m3d,  only: T,S
#endif
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: mode
   character(len=*), intent(in)        :: fname
   integer, intent(in)                 :: runtype
!
! !INPUT/OUTPUT PARAMETERS:
   integer, intent(inout)              :: loop
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES
   integer, save             :: n=0
   logical, save             :: continuous=.false.
   integer                   :: n,jd,secs 
   REALTYPE                  :: dt
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'restart_file() # ',Ncall
#endif

   if (mode .eq. WRITING) then
      n = n + 1
      LEVEL2 'Saving hotstart file # ',n,' as ',fname
      open(RESTART,file=fname,status='unknown',form='unformatted')
      LEVEL3 'saving loop, julianday, secondsofday and timestep'
      write(RESTART) loop,julianday,secondsofday,timestep
      LEVEL3 'saving basic variables'
      write(RESTART) z,zo,U,zu,zub,SlUx,Slru,V,zv,zvb,SlVx,Slrv
#ifndef NO_3D
      if (runtype .ge. 2)  then
         LEVEL3 'saving 3D barotropic variables'
!kbk         write(RESTART) Uint,Vint
         write(RESTART) uu,vv
         write(RESTART) tke,eps
         write(RESTART) num,nuh
#ifndef NO_BAROCLINIC
         if(runtype .ge. 3) then
            LEVEL3 'saving 3D baroclinic variables'
            write(RESTART) T,S
         end if
#endif
      end if
#endif
      close(RESTART)
   end if

   if (mode .eq. READING) then
      LEVEL2 'Reading hotstart file: '
      LEVEL3 trim(fname)
      open(RESTART,file=fname,status='unknown',form='unformatted')
      LEVEL3 'reading loop, julianday, secondsofday and timestep'
!KBK      read(RESTART) loop,julianday,secondsofday,timestep
      read(RESTART) n,jd,secs,dt
      LEVEL3 'reading basic variables'
      read(RESTART) z,zo,U,zu,zub,SlUx,Slru,V,zv,zvb,SlVx,Slrv
#ifndef NO_3D
      if (runtype .ge. 2)  then
!KBK This needs to be changed !!!! KBK
!Only works because E2DFIELD = I2DFIELD
      ssen=z
      ssun=zu
      ssvn=zv
         LEVEL3 'reading 3D barotropic variables'
!kbk         read(RESTART) Uint,Vint
         read(RESTART) uu,vv
         read(RESTART) tke,eps
         read(RESTART) num,nuh
#ifndef NO_BAROCLINIC
         if(runtype .ge. 3) then
            LEVEL3 'reading 3D baroclinic variables'
            read(RESTART) T,S
         end if
#endif
      end if
#endif
      close(RESTART)
!     make some sanity checks
      if (continuous) then
         loop=n; 
         julianday=jd; secondsofday=secs; timestep=dt;
      else
         if (jd .ne. julianday .or. secs .ne. secondsofday) then
            FATAL 'start time given in getm.inp does not match time'
            FATAL 'read from hot-start file'
            FATAL 'from getm.inp: ',julianday,secondsofday
            FATAL 'from hotstart: ',jd,secs
            stop 'restart_file()'
         end if
         if (dt .ne. timestep) then
            LEVEL3 ''
            LEVEL3 'INFO:'
            LEVEL3 'time step changed between hotstart file and value '
            LEVEL3 'given in getm.inp (this is OK - but beware when post-'
            LEVEL3 'processing)'
            LEVEL3 ''
         end if
         loop = 0
         julianday=jd; secondsofday=secs
      end if
   end if
#ifdef DEBUG
   write(debug,*) 'Leaving restart_file()'
   write(debug,*)
#endif
   return
   end subroutine restart_file
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clean_output - cleans up after run
!
! !INTERFACE:
   subroutine clean_output()
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Writes calculated fields to files.
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
   REALTYPE :: dummy
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'clean_output() # ',Ncall
#endif

   select case (out_fmt)
      case (NETCDF)
         dummy = -_ONE_*1.
         call save_2d_ncdf(dummy)
         call ncdf_close()
      case DEFAULT
         STDERR 'Fatal error: A non valid input format has been chosen'
         stop 'clean_output'
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving clean_output()'
   write(debug,*)
#endif
   return
   end subroutine clean_output

!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: file_names - setup output file names
!
! !INTERFACE:
   subroutine file_names(runid,myid)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        ::  runid
   integer, intent(in)                 ::  myid
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
   character(len=3)                    :: buf
   character(len=16)                   :: pid,ext
   character(len=PATH_MAX)             :: fname
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (myid .ge. 0) then
      write(buf,'(I3.3)') myid
      pid = '.' // trim(buf)
   else
      pid = ''
   end if

   hot_out = trim(out_dir) //'/'// 'restart' // trim(pid) // '.out'
   if (out_fmt .eq. NETCDF) then
      ext = 'nc'
      out_f_2d =  &
         trim(out_dir) //'/'// trim(runid) // '.2d' // trim(pid) // '.' // ext
      out_f_3d =  &
         trim(out_dir) //'/'// trim(runid) // '.3d' // trim(pid) // '.' // ext
   end if

   return
   end subroutine file_names
!EOC

!-----------------------------------------------------------------------

   end module output

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
