!$Id: output.F90,v 1.1 2002-05-02 14:01:52 gotm Exp $
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
   use commhalo, only: myid
   use time, only: write_time_string,timestep,timestr
   use ncdf_out
   use ascii_out
#ifdef TEST_NESTING
   use nesting
#endif
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   integer			:: out_fmt=NETCDF
   character(LEN = PATH_MAX)	:: in_dir='.'
   character(LEN = PATH_MAX)	:: out_dir='.'
   character(LEN = PATH_MAX)	:: out_f_2d
   character(LEN = PATH_MAX)	:: out_f_3d

   logical	:: save_meteo=.false.
   logical	:: save_2d=.true.
   logical	:: save_3d=.true.
   logical	:: save_vel=.true.
   logical	:: save_strho=.true.
   logical	:: save_s=.true.
   logical	:: save_t=.true.
   logical	:: save_rho=.true.
   logical	:: save_turb=.true.
   logical	:: save_tke=.true.
   logical	:: save_eps=.true.
   logical	:: save_num=.true.
   logical	:: save_nuh=.true.
   logical	:: save_spm=.false.
   integer	:: first_2d=1
   integer	:: step_2d=1
   integer	:: first_3d=1
   integer	:: step_3d=1
   integer	:: hotout=-1
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: output.F90,v $
!  Revision 1.1  2002-05-02 14:01:52  gotm
!  Initial revision
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
   subroutine init_output(runid,title,starttime,runtype,dryrun)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)		:: runid,title,starttime
   integer, intent(in)			:: runtype
   logical, intent(in)			:: dryrun
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
   namelist /io_spec/ out_fmt,					&
                     in_dir,out_dir,				&
                     save_2d,save_3d,save_vel,			&
                     save_strho,save_s,save_t,save_rho,		&
                     save_turb,save_tke,save_eps,save_num,save_nuh,	&
                     save_spm,					&
                     first_2d,step_2d,first_3d,step_3d,hotout,	&
		     save_meteo
!   logical	:: nesting=.true.
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

   call file_names(runid)

   if(save_2d) then
      LEVEL2 '2D results: ',TRIM(out_f_2d)
      LEVEL2 'First=',first_2d,' step=',step_2d
      if(save_meteo) then
         LEVEL2 'Saving meteo forcing in ',trim(out_f_2d)
      end if
   end if
   if(save_3d) then
      LEVEL2 '3D results: ',TRIM(out_f_3d)
      LEVEL2 'First=',first_3d,' step=',step_3d
   end if

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
            if (save_3d) call init_3d_ncdf(out_f_3d,title,starttime)
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
   integer, intent(in)	:: runtype,n
   REALTYPE, intent(in)	:: timestep
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
   REALTYPE	:: secs
   logical	:: write_2d,write_3d
   integer	:: dummy
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
      call divergence()
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
            if (write_3d) call save_3d_ncdf(secs)
         case DEFAULT
           STDERR 'Fatal error: A non valid input format has been chosen'
           stop 'do_output'
      end select
   end if

!  Restart file
   if (hotout .gt. 0 .and. mod(n,hotout) .eq. 0) then
      dummy = n
      call restart_file(WRITING,trim(out_dir) // '/' // 'restart.out',dummy,runtype)
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
   use m3d,  only: uu,vv,tke,eps,num,nuh,T,S,ssen,ssun,ssvn
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)		:: mode
   character(len=*), intent(in)	:: fname
   integer, intent(in)		:: runtype
!
! !INPUT/OUTPUT PARAMETERS:
   integer, intent(inout)	:: loop
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES
   integer, save	:: n=0
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
      if (runtype .gt. 1)  then
         LEVEL3 'saving 3D barotropic variables'
!kbk         write(RESTART) Uint,Vint
         write(RESTART) uu,vv
         write(RESTART) tke,eps
         write(RESTART) num,nuh
         if(runtype .ge. 3) then
            LEVEL3 'saving 3D baroclinic variables'
            write(RESTART) T,S
         end if
      end if
      close(RESTART)
   end if

   if (mode .eq. READING) then
      LEVEL2 'Reading hotstart file: '
      LEVEL3 trim(fname)
      open(RESTART,file=fname,status='unknown',form='unformatted')
      LEVEL3 'reading loop, julianday, secondsofday and timestep'
      read(RESTART) loop,julianday,secondsofday,timestep
      LEVEL3 'reading basic variables'
      read(RESTART) z,zo,U,zu,zub,SlUx,Slru,V,zv,zvb,SlVx,Slrv
      if (runtype .gt. 1)  then
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
         if(runtype .ge. 3) then
            LEVEL3 'reading 3D baroclinic variables'
            read(RESTART) T,S
         end if
      end if
      close(RESTART)
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
         call save_2d_ncdf( -_ONE_ )
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
   subroutine file_names(runid)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)		::  runid
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
   logical 				:: TO_FILE=.false.
   character(len=3)			:: buf
   character(len=16)			:: pid,ext

   character(len=PATH_MAX)		:: fname
!
!EOP
!-----------------------------------------------------------------------
!BOC

   if (myid.ge.0) then
      write(buf,'(I3.3)') myid
      pid = '.' // TRIM(buf)
   else
      pid = ''
   end if

!  stdin, stdout, stderr and debug files
   if (TO_FILE) then
      ext   = 'stderr'
      fname = TRIM(runid) // TRIM(pid) // '.' // ext
      open(stderr,file=Fname)

      ext   = 'stdout'
      fname = TRIM(runid) // TRIM(pid) // '.' // ext
      open(stdout,file=Fname)
   end if

   if (out_fmt .eq. NETCDF) then
      ext = 'nc'
      out_f_2d =  &
         TRIM(out_dir) //'/'// TRIM(runid) // '.2d' // TRIM(pid) // '.' // ext
      out_f_3d =  &
         TRIM(out_dir) //'/'// TRIM(runid) // '.3d' // TRIM(pid) // '.' // ext
   end if

   return
   end subroutine file_names
!EOC

!-----------------------------------------------------------------------

   end module output

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
