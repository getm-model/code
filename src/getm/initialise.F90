#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  initialise - setup the entire model
!
! !INTERFACE:
   module initialise
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   public                              :: init_model
   integer                             :: runtype=1
   logical                             :: dryrun=.false.
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
! !IROUTINE: init_model - initialise getm
!
! !INTERFACE:
   subroutine init_model(dstr,tstr)
!
! !USES:
   use kurt_parallel, only: init_parallel,myid
#ifdef GETM_PARALLEL
   use halo_mpi, only: init_mpi,print_MPI_info
#endif
   use output, only: init_output,do_output,restart_file,out_dir
   use input,  only: init_input
   use domain, only: init_domain
   use domain, only: iextr,jextr,imin,imax,jmin,jmax,kmax
   use domain, only: vert_cord,maxdepth
   use time, only: init_time,update_time,write_time_string
   use time, only: start,timestr,timestep
   use m2d, only: init_2d,postinit_2d, z
   use les, only: init_les
   use getm_timers, only: init_getm_timers, tic, toc, TIM_INITIALIZE
#ifndef NO_3D
   use m3d, only: init_3d,postinit_3d, ssen,ssun,ssvn
#ifndef NO_BAROCLINIC
   use m3d, only: T
#endif
   use turbulence, only: init_turbulence
   use mtridiagonal, only: init_tridiagonal
   use rivers, only: init_rivers
   use variables_3d, only: avmback,avhback
#ifdef SPM
   use suspended_matter, only: init_spm
#endif
#ifdef _FABM_
   use getm_fabm, only: fabm_calc
   use getm_fabm, only: init_getm_fabm, postinit_getm_fabm
   use rivers, only: init_rivers_fabm
#endif
#ifdef GETM_BIO
   use bio, only: bio_calc
   use getm_bio, only: init_getm_bio
   use rivers, only: init_rivers_bio
#endif
#endif
   use meteo, only: init_meteo,do_meteo
   use integration,  only: MinN,MaxN
   use exceptions
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*)                    :: dstr,tstr
!
! !DESCRIPTION:
!  Reads the namelist and makes calls to the init functions of the
!  various model components.
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
   integer:: i,j
   character(len=8)          :: buf
   character(len=64)         :: runid
   character(len=80)         :: title
   logical                   :: parallel=.false.
   logical                   :: hotstart=.false.
   logical                   :: use_epoch=.false.
   logical                   :: save_initial=.false.
#if (defined GETM_PARALLEL && defined INPUT_DIR)
   character(len=PATH_MAX)   :: input_dir=INPUT_DIR
#else
   character(len=PATH_MAX)   :: input_dir='./'
#endif
   character(len=PATH_MAX)   :: hot_in=''

   namelist /param/ &
             dryrun,runid,title,parallel,runtype,  &
             hotstart,use_epoch,save_initial
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_model() # ',Ncall
#endif
#ifndef NO_TIMERS
   call init_getm_timers()
#endif
   ! Immediately start to time (rest of) init:
   call tic(TIM_INITIALIZE)

   ! We need to pass info about the input directory
#if 0
   call getarg(1,base_dir)
   if(len_trim(base_dir) .eq. 0) then
      call getenv("base_dir",base_dir)
   end if
   if(len_trim(base_dir) .gt. 0) then
      base_dir = trim(base_dir) // '/'
   end if
#endif

!
! In parallel mode it is imperative to let the instances
! "say hello" right away. For MPI this changes the working directory,
! so that input files can be read.
!
#ifdef GETM_PARALLEL
   call init_mpi()
#endif

#if (defined GETM_PARALLEL && defined INPUT_DIR)
   STDERR 'input_dir:'
   STDERR trim(input_dir)
#endif
!
! Open the namelist file to get basic run parameters.
!
   title='A descriptive title can be specified in the param namelist'
   open(NAMLST,status='unknown',file=trim(input_dir) // "/getm.inp")
   read(NAMLST,NML=param)

#ifdef NO_BAROCLINIC
   if(runtype .ge. 3) then
      FATAL 'getm not compiled for baroclinic runs'
      stop 'init_model()'
   end if
#endif

#ifdef NO_3D
   if(runtype .ge. 2) then
      FATAL 'getm not compiled for 3D runs'
      stop 'init_model()'
   end if
#endif

! call all modules init_ ... routines

   if (parallel) then
#ifdef GETM_PARALLEL
      call init_parallel(runid,input_dir)
#else
      STDERR 'You must define GETM_PARALLEL and recompile'
      STDERR 'in order to run in parallel'
      stop 'init_model()'
#endif
   end if

#if (defined GETM_PARALLEL && defined SLICE_MODEL)
    call getm_error('init_model()', &
         'SLICE_MODEL does not work with GETM_PARALLEL - for now')
#endif

   STDERR LINE
   STDERR 'getm ver. ',RELEASE,': Started on  ',dstr,' ',tstr
   STDERR LINE
   STDERR 'Initialising....'
   STDERR LINE
   LEVEL1 'the run id is: ',trim(runid)
   LEVEL1 'the title is:  ',trim(title)

   select case (runtype)
      case (1)
         LEVEL1 '2D run (hotstart=',hotstart,')'
      case (2)
         LEVEL1 '3D run - no density (hotstart=',hotstart,')'
      case (3)
         LEVEL1 '3D run - frozen density (hotstart=',hotstart,')'
      case (4)
         LEVEL1 '3D run - full (hotstart=',hotstart,')'
      case default
         FATAL 'A non valid runtype has been specified.'
         stop 'initialise()'
   end select

   call init_time(MinN,MaxN)
   if(use_epoch) then
      LEVEL2 'using "',start,'" as time reference'
   end if

   call init_domain(input_dir)

   call init_meteo(hotstart)

#ifndef NO_3D
   call init_rivers(hotstart)
#endif

   call init_2d(runtype,timestep,hotstart)

#ifndef NO_3D
   if (runtype .gt. 1) then
      call init_3d(runtype,timestep,hotstart)
#ifndef CONSTANT_VISCOSITY
      call init_turbulence(60,trim(input_dir) // 'gotmturb.nml',kmax)
#else
      LEVEL3 'turbulent viscosity and diffusivity set to constant (-DCONSTANT_VISCOSITY)'
#endif
      LEVEL2 'background turbulent viscosity set to',avmback
      LEVEL2 'background turbulent diffusivity set to',avhback
      call init_tridiagonal(kmax)

#ifdef SPM
      call init_spm(trim(input_dir) // 'spm.inp',runtype)
#endif
#ifdef _FABM_
      call init_getm_fabm(trim(input_dir) // 'getm_fabm.inp')
      call init_rivers_fabm
#endif
#ifdef GETM_BIO
      call init_getm_bio(trim(input_dir) // 'getm_bio.inp')
      call init_rivers_bio
#endif
   end if
#endif

   call init_les(runtype)

   call init_output(runid,title,start,runtype,dryrun,myid,MinN,MaxN,save_initial)

   close(NAMLST)

#if 0
   call init_waves(hotstart)
   call init_biology(hotstart)
#endif

   if (hotstart) then
      LEVEL1 'hotstart'
      if (myid .ge. 0) then
         write(buf,'(I3.3)') myid
         buf = '.' // trim(buf) // '.in'
      else
         buf = '.in'
      end if
      hot_in = trim(out_dir) //'/'// 'restart' // trim(buf)
      call restart_file(READING,trim(hot_in),MinN,runtype,use_epoch)
      LEVEL3 'MinN adjusted to ',MinN
      call update_time(MinN)
      call write_time_string()
      LEVEL3 timestr
      MinN = MinN+1
   end if

   call postinit_2d(runtype,timestep,hotstart)
#ifndef NO_3D
   if (runtype .gt. 1) then
      call postinit_3d(runtype,MinN-1,hotstart)
#ifdef _FABM_
      if (fabm_calc) call postinit_getm_fabm()
#endif
   end if
#endif

   call init_input(input_dir,MinN)

   call toc(TIM_INITIALIZE)
   ! The rest is timed with meteo and output.

   if(runtype .le. 2) then
      call do_meteo(MinN)
#ifndef NO_3D
#ifndef NO_BAROCLINIC
   else
      call do_meteo(MinN,T(:,:,kmax))
#endif
#endif
   end if

   if (.not. dryrun) then
      call do_output(runtype,MinN-1,timestep)
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving init_model()'
   write(debug,*)
#endif
   return
   end subroutine init_model
!EOC

!-----------------------------------------------------------------------

   end module initialise

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
