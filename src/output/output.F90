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
   use ascii_out
   use m2d, only: sealevel_check
#ifndef NO_3D
   use m3d, only: calc_salt,calc_temp
   use variables_3d, only: do_numerical_analyses
#endif
#ifdef TEST_NESTING
   use nesting
#endif
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   integer                             :: out_fmt=NETCDF
   integer                             :: hotin_fmt=NETCDF
   integer                             :: hotout_fmt=NETCDF
   character(LEN = PATH_MAX)           :: in_dir='.'
   character(LEN = PATH_MAX)           :: out_dir='.'
   character(LEN = PATH_MAX)           :: out_f_2d
   character(LEN = PATH_MAX)           :: out_f_3d
   character(LEN = PATH_MAX)           :: out_f_mean
   character(LEN = PATH_MAX)           :: hot_out

   logical                             :: save_metrics=.false.
   logical                             :: save_masks=.false.
   logical                             :: save_2d=.true.
   logical                             :: save_meteo=.false.
   logical                             :: save_3d=.true.
   logical                             :: save_mean=.false.
   logical                             :: save_vel=.true.
   logical                             :: destag=.false.
   logical                             :: save_strho=.false.
   logical                             :: save_s=.false.
   logical                             :: save_t=.false.
   logical                             :: save_rho=.false.
   logical                             :: save_rad=.false.
   logical                             :: save_turb=.true.
   logical                             :: save_tke=.true.
   logical                             :: save_eps=.true.
   logical                             :: save_num=.true.
   logical                             :: save_nuh=.true.
   logical                             :: save_ss_nn=.false.
   logical                             :: save_taub=.false.
   integer                             :: first_2d=1
   integer                             :: step_2d=1
   integer                             :: sync_2d=1
   integer                             :: first_3d=1
   integer                             :: step_3d=1
   integer                             :: sync_3d=1
   integer                             :: hotout(3)=-1
   integer                             :: mean0=0
   integer                             :: meanout=-1
   logical                             :: save_numerical_analyses=.false.
   logical,private                     :: save_restart
   integer,private                     :: firstN=-1
   integer,private                     :: lastN=-1
   logical,private                     :: save_init=.false.

!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_output - initialise all external files and units
!
! !DESCRIPTION:
!
! !INTERFACE:
   subroutine init_output(runid,title,starttime,runtype,dryrun,myid,MinN,MaxN,save_initial)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: runid,title,starttime
   integer, intent(in)                 :: runtype,myid,MinN,MaxN
   logical, intent(in)                 :: dryrun,save_initial
!
! !REVISION HISTORY:
!
!  See log for module
!
! !LOCAL VARIABLES:
   namelist /io_spec/ &
             out_fmt,hotin_fmt,hotout_fmt, &
             in_dir,out_dir, save_metrics, save_masks, &
             save_2d,save_3d,save_vel,destag, &
             save_strho,save_s,save_t,save_rho,save_rad, &
             save_turb,save_tke,save_eps,save_num,save_nuh, &
             save_ss_nn,save_taub, &
             first_2d,step_2d,sync_2d,first_3d,step_3d,sync_3d,hotout, &
             meanout, save_meteo, save_numerical_analyses
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

   firstN = MinN
   lastN = MaxN
   save_init = save_initial

   read(NAMLST, nml=io_spec)

   if (hotin_fmt .ne. NETCDF) then
     LEVEL2 'WARNING: Support of non-netcdf restart files will be stopped.'
     LEVEL2 '         Do a zero-length simulation to convert your restart files to netcdf!'
   end if

   LEVEL2 'save_nuh',save_nuh
   LEVEL2 'save_num',save_num
   LEVEL2 'save_tke',save_tke

   if (runtype .eq. 1) then
      save_3d = .false.
      save_vel = .false.
      save_strho = .false.
      save_turb = .false.
   end if

   if (runtype .eq. 2) then
      save_strho = .false.
   end if

   if (.not. save_strho) then
      save_s = .false.
      save_t = .false.
      save_rho = .false.
   end if

#ifndef NO_3D
   if (.not. calc_salt) then
      save_s = .false.
   end if
   if (.not. calc_temp) then
      save_t = .false.
   end if
#endif

   if(destag) then
      LEVEL2 'de-stag velocities to T-points'
   else
      LEVEL2 'keeping velocities on calculation grid'
   end if

   call file_names(runid,myid)

   if(save_2d) then
      LEVEL2 '2D results: ',trim(out_f_2d)
      if (sync_2d .lt. 0) sync_2d=0
      LEVEL2 'First=',first_2d,' step=',step_2d,' sync_2d= ',sync_2d
      if(save_meteo) then
         LEVEL2 'Saving meteo forcing in ',trim(out_f_2d)
      end if
   end if
#ifndef NO_3D
   if(save_3d) then
      LEVEL2 '3D results: ',trim(out_f_3d)
      if (sync_3d .lt. 0) sync_3d=0
      LEVEL2 'First=',first_3d,' step=',step_3d,' sync_3d= ',sync_3d
   end if
#endif

   save_restart = (hotout(1).ge.0)
   if (save_restart) then
      if (hotout_fmt .ne. NETCDF) then
        STDERR 'Writing of non-netcdf restart files not supported anymore!'
        stop
      end if
   end if
   if ( hotout(1) .gt. 0 .and. hotout(2) .lt. hotout(1) ) then
      if ( hotout(2) .eq. -1 ) then
         hotout(2) = 2147483647
         hotout(3) = hotout(1)
      else
         hotout(2) = hotout(1)
         hotout(3) = 1
      end if
   end if
   if ( hotout(1) .eq. hotout(2)) then
      hotout(3) = 1
   end if

#ifndef NO_3D
   save_mean=(meanout .ge. 0 .and. runtype .gt. 1)
   if ( save_mean ) then
      LEVEL2 'Mean fields in: ',trim(out_f_mean)
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
!KB            if (save_mean)  &
!KB                     call init_mean_ncdf(out_f_mean,title,starttime)
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

#ifndef NO_3D
   if (save_numerical_analyses) then
      LEVEL2 "calculate and save mixing analysis"
      do_numerical_analyses=.true.
   end if
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
   use getm_timers, only: tic, toc, TIM_OUTPUT
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: runtype,n
   REALTYPE, intent(in)                :: timestep
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
   REALTYPE                  :: secs
   logical                   :: write_2d,write_3d,write_mean=.false.,write_restart=.false.
   integer                   :: dummy
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_output() # ',Ncall
#endif
   call tic(TIM_OUTPUT)

   write_2d = save_2d .and. n .ge. first_2d .and. mod(n,step_2d).eq.0
   write_3d = save_3d .and. n .ge. first_3d .and. mod(n,step_3d).eq.0

!  TODO: Presently save_init can only be used to switch off initial output.
!        Maybe we want to extend this so that save_init=T in any case
!        (independent of the checks above) causes initial output?

   if (.not.save_init .and. n.eq.firstN-1) then
      write_2d = .false.
      write_3d = .false.
   end if

#ifndef NO_3D
   if (save_mean .and. n.gt.mean0) then
      if (meanout .eq. 0) then
         write_mean = (n.eq.lastN)
      else
         write_mean = (mod(n,meanout).eq.0)
         !write_mean = (mod(n-mean0,meanout).eq.0)
      end if
!KB      call calc_mean_fields(n,write_mean)
   end if
#endif

   if (write_2d .or. write_3d .or. write_mean) then
      call write_time_string()
      if (write_2d)   LEVEL3 timestr, ': saving 2D .... '
      if (write_3d)   LEVEL3 timestr, ': saving 3D .... '
      if (write_mean) LEVEL3 timestr, ': saving mean fields .... '

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
!KB            if (write_mean) call save_mean_ncdf(secs)
#endif
         case DEFAULT
           STDERR 'Fatal error: A non valid input format has been chosen'
           stop 'do_output'
      end select
   end if

!  Restart file
   if (save_restart) then
      if (hotout(1) .eq. 0) then
!        Save last restart file
!        also works for zero-length simulations (called from init_model)
         write_restart = (n.eq.lastN)
      else if (firstN .le. n) then ! avoid recreating just read restart file
         write_restart = hotout(1).le.n .and. n.le.hotout(2) .and. mod(n,hotout(3)).eq.0
      end if
      if (write_restart) then
         if ( sealevel_check .ne. 0 ) then
            LEVEL2 'Checking for NANs before saving hotstart file...'
            call sealevel_nan_check()
         end if
         dummy = n
         call restart_file(WRITING,trim(hot_out),dummy,runtype)
      end if
   end if

   call toc(TIM_OUTPUT)
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
   subroutine restart_file(mode,fname,loop,runtype,use_epoch)
!
! !DESCRIPTION:
!  This routine write the variables necesaary for a 'hot' model start.
!
! !USES:
   use time, only: timestep,julianday,secondsofday
   use variables_2d, only: z,zo
   use domain, only: au,av
   use domain, only: imin,imax,jmin,jmax
#ifdef ZUB_ZVB
   use variables_2d, only: U,fU,res_du,SlUx,Slru
   use variables_2d, only: V,fV,res_dv,SlVx,Slrv
#else
   use variables_2d, only: U,fU,SlUx,Slru
   use variables_2d, only: V,fV,SlVx,Slrv
#endif
   use variables_2d, only: Uinto,Vinto
#ifndef NO_3D
   use variables_3d, only: ssen,ssun,ssvn
   use variables_3d, only: sseo,ssuo,ssvo
   use variables_3d, only: uu,vv,ww
   use variables_3d, only: uuEx,vvEx
   use variables_3d, only: tke,eps,num,nuh
   use variables_3d, only: hn
#ifndef NO_BAROCLINIC
   use variables_3d, only: T,S
#endif
#endif
#ifdef SPM
  use suspended_matter, only: spm_calc,spm_hotstart
  use variables_3d, only: spm,spm_pool
#endif
#ifdef GETM_BIO
  use bio, only: bio_calc
  use bio_var, only: numc
  use getm_bio, only: bio_init_method
  use variables_3d,  only: cc3d
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: mode
   character(len=*), intent(in)        :: fname
   integer, intent(in)                 :: runtype
   logical, optional, intent(in)       :: use_epoch
!
! !INPUT/OUTPUT PARAMETERS:
   integer, intent(inout)              :: loop
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES
   integer, save             :: n=0
   integer                   :: i,j
   integer                   :: jd,secs
   character(len=19)         :: timestr_out
   REALTYPE                  :: dt
#ifdef _HOT_ZU_ZV_
   REALTYPE,dimension(E2DFIELD) :: wrk2d
#endif
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
      select case (hotout_fmt)
         case(NETCDF)
            if (n .eq. 1) then
               call create_restart_ncdf(fname,loop,runtype)
            end if
            call write_restart_ncdf(runtype,_ZERO_,loop,julianday,secondsofday)
         case(BINARY)
            open(RESTART,file=fname,status='unknown',form='unformatted')
            LEVEL3 'saving loop, julianday, secondsofday and timestep'
            write(RESTART) loop,julianday,secondsofday,timestep
            LEVEL3 'saving basic variables'
            write(RESTART) z,zo,U
#ifdef _HOT_ZU_ZV_
            LEVEL3 'obsolete saving of placeholder for zu'
            write(RESTART) wrk2d
#endif
#ifdef ZUB_ZVB
            LEVEL3 'obsolete saving of res_du'
            write(RESTART) res_du
#endif
            write(RESTART) SlUx,Slru,V
#ifdef _HOT_ZU_ZV_
            LEVEL3 'obsolete saving of placeholder for zv'
            write(RESTART) wrk2d
#endif
#ifdef ZUB_ZVB
            LEVEL3 'obsolete saving of res_dv'
            write(RESTART) res_dv
#endif
            write(RESTART) SlVx,Slrv
#ifndef NO_3D
            if (runtype .ge. 2)  then
               LEVEL3 'saving 3D barotropic variables'
               write(RESTART) ssen,ssun,ssvn
               write(RESTART) sseo,ssuo,ssvo
               write(RESTART) Uinto,Vinto
               write(RESTART) uu,vv,ww
               write(RESTART) uuEx,vvEx
               write(RESTART) tke,eps
               write(RESTART) num,nuh
               write(RESTART) hn
#ifndef NO_BAROCLINIC
               if (runtype .ge. 3) then
                  LEVEL3 'saving 3D baroclinic variables'
                  write(RESTART) T,S
               end if
#endif
#ifdef SPM
               if (spm_calc) then
                  LEVEL3 'saving spm'
                  write(RESTART) spm
                  write(RESTART) spm_pool
               end if
#endif
#ifdef GETM_BIO
               if (bio_calc) then
                  LEVEL3 'saving bio variables'
                  write(RESTART) cc3d
               end if
#endif
            end if
#endif
            close(RESTART)
         case DEFAULT
            STDERR 'Fatal error: A non valid hot format has been chosen'
            stop 'clean_output'
      end select
   end if ! WRITING

   if (mode .eq. READING) then
      LEVEL2 'Reading hotstart file: '
      LEVEL3 trim(fname)
      select case (hotin_fmt)
         case(NETCDF)
            call open_restart_ncdf(fname,runtype)
            call read_restart_ncdf(runtype,j,jd,secs,dt)
         case(BINARY)
            open(RESTART,file=fname,status='unknown',form='unformatted')
            LEVEL3 'reading loop, julianday, secondsofday and timestep'
            read(RESTART) j,jd,secs,dt
            LEVEL3 'reading basic variables'
            read(RESTART) z,zo,U
#ifdef _HOT_ZU_ZV_
            LEVEL3 'obsolete reading of placeholder for zu'
            read(RESTART) wrk2d
#endif
#ifdef ZUB_ZVB
            LEVEL3 'obsolete reading of placeholder for res_du'
            read(RESTART) res_du
            res_du=_ZERO_
#endif
            read(RESTART) SlUx,Slru,V
#ifdef _HOT_ZU_ZV_
            LEVEL3 'obsolete reading of placeholder for zv'
            read(RESTART) wrk2d
#endif
#ifdef ZUB_ZVB
            LEVEL3 'obsolete reading of placeholder for res_dv'
            read(RESTART) res_dv
            res_dv=_ZERO_
#endif
            read(RESTART) SlVx,Slrv
            where(au .eq. 0) U=_ZERO_
            where(av .eq. 0) V=_ZERO_
#ifndef NO_3D
            if (runtype .ge. 2)  then
               LEVEL3 'reading 3D barotropic variables'
               read(RESTART) ssen,ssun,ssvn
               read(RESTART) sseo,ssuo,ssvo
               read(RESTART) Uinto,Vinto
               read(RESTART) uu,vv,ww
               read(RESTART) uuEx,vvEx
               read(RESTART) tke,eps
               read(RESTART) num,nuh
               read(RESTART) hn
               forall(i=imin-HALO:imax+HALO,j=jmin-HALO:jmax+HALO, au(i,j).EQ.0) &
                    uu(i,j,:)=_ZERO_
               forall(i=imin-HALO:imax+HALO,j=jmin-HALO:jmax+HALO, av(i,j).EQ.0) &
                    vv(i,j,:)=_ZERO_
#ifndef NO_BAROCLINIC
               if(runtype .ge. 3) then
                  LEVEL3 'reading 3D baroclinic variables'
                  read(RESTART) T,S
               end if
#endif
#ifdef SPM
               if(spm_calc) then
                  if (spm_hotstart) then
                     LEVEL3 'reading spm variables'
                     read(RESTART) spm
                     read(RESTART) spm_pool
                  else
                     LEVEL3 'spm variables not read from hotstart file'
                     LEVEL3 'set spm_init_method=0 to read them from hotstart file'
                  end if
               end if
#endif

#ifdef GETM_BIO
               if(bio_calc .and. bio_init_method .eq. 0) then
                  LEVEL3 'reading bio variables'
                  read(RESTART) cc3d
               end if
#endif
            end if
#endif
            close(RESTART)
         case DEFAULT
            STDERR 'Fatal error: A non valid hot format has been chosen'
            stop 'clean_output'
      end select

!     make some sanity checks
      if (use_epoch) then
         if (jd .eq. julianday .and. secs .eq. secondsofday) then
            loop=0
         else
            loop=j;
            julianday=jd; secondsofday=secs; timestep=dt;
         end if
      else
         if (jd .ne. julianday .or. secs .ne. secondsofday) then
            FATAL 'start time given in getm.inp does not match time'
            FATAL 'read from hot-start file'
            FATAL 'from getm.inp: ',julianday,secondsofday
            call write_time_string(julianday,secondsofday,timestr_out)
            LEVEL3 timestr_out
            FATAL 'from hotstart: ',jd,secs
            call write_time_string(jd,secs,timestr_out)
            LEVEL3 timestr_out
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
      firstN = loop+1
      mean0 = loop
   end if ! READING

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
   subroutine clean_output(runtype,loop)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: runtype
!
! !INPUT/OUTPUT PARAMETERS:
   integer, intent(inout)              :: loop
!
! !DESCRIPTION:
!  Writes calculated fields to files.
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
   integer :: zero=0
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
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
   character(len=4)                    :: buf
   character(len=16)                   :: pid,ext
   character(len=PATH_MAX)             :: fname
!EOP
!-----------------------------------------------------------------------
!BOC
   if (myid .ge. 0) then
      write(buf,'(I4.4)') myid
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
      out_f_mean =  &
         trim(out_dir) //'/'// trim(runid) // '.mean' // trim(pid) // '.' // ext
   end if

   return
   end subroutine file_names
!EOC


!-----------------------------------------------------------------------

   end module output

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!----------------------------------------------------------------------
