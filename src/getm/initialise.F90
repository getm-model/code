!$Id: initialise.F90,v 1.3 2003-04-23 12:03:46 kbk Exp $
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
!  $Log: initialise.F90,v $
!  Revision 1.3  2003-04-23 12:03:46  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 16:39:16  kbk
!  parallel support, NO_3D
!
!  Revision 1.1.1.1  2002/05/02 14:01:25  gotm
!  recovering after CVS crash
!
!  Revision 1.10  2001/10/22 12:24:07  bbh
!  Typo
!
!  Revision 1.9  2001/10/22 11:45:27  bbh
!  Only save initial fields if not dryrun
!
!  Revision 1.8  2001/09/18 17:52:03  bbh
!  save_initial now in namelist
!
!  Revision 1.7  2001/09/18 17:48:32  bbh
!  Added algoritm for rivers - getting river data still missing
!
!  Revision 1.6  2001/09/13 15:00:28  bbh
!  Use the new ncdf-scheme
!
!  Revision 1.5  2001/07/26 13:57:14  bbh
!  Meteo working - needs some polishing
!
!  Revision 1.4  2001/06/04 13:09:53  bbh
!  Includes - dryrun - in call to init_output()
!
!  Revision 1.3  2001/05/03 20:11:11  bbh
!  Use runtype in init_3d
!
!  Revision 1.2  2001/04/24 08:24:58  bbh
!  Use runtype instead of macro
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
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
! !IROUTINE: init_model - initialise getm
!
! !INTERFACE:
   subroutine init_model(dstr,tstr)
!
! !USES:
   use kurt_parallel, only: init_parallel,myid
   use output, only: init_output,do_output,restart_file
   use input,  only: init_input
   use domain, only: init_domain
   use domain, only: iextr,jextr,imin,imax,jmin,jmax
   use domain, only: iimin,iimax,jjmin,jjmax,kmax
   use domain, only: vert_cord
   use time, only: init_time,update_time,write_time_string
   use time, only: start,timestr,timestep
   use m2d, only: init_2d,z,zu,zv
#ifndef NO_3D
   use m3d, only: cord_relax,init_3d,ssen,ssun,ssvn
#ifndef NO_BAROCLINIC
   use m3d, only: T
#endif
   use turbulence, only: init_turbulence
   use mtridiagonal, only: init_tridiagonal
   use rivers, only: init_rivers
#endif
   use meteo, only: init_meteo,do_meteo
   use integration,  only: MinN,MaxN
#ifndef NO_BAROCLINIC
   use eqstate, only: do_eqstate
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*)                    :: dstr,tstr
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
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
   logical                   :: save_initial=.false.
#if (defined PARALLEL && defined INPUT_DIR)
   character(len=PATH_MAX)   :: input_dir=INPUT_DIR
#else
   character(len=PATH_MAX)   :: input_dir='./'
#endif
   character(len=PATH_MAX)   :: hot_in=''

   namelist /param/ &
             dryrun,runid,title,parallel,runtype,  &
             hotstart,save_initial
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_model() # ',Ncall
#endif

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

#if (defined PARALLEL && defined INPUT_DIR)
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
#ifdef PARALLEL
      call init_parallel(runid,input_dir)
#else
      STDERR 'You must define GETM_PARALLEL and recompile'
      STDERR 'in order to run in parallel'
      stop 'init_model()'
#endif
   end if
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

   call init_domain(input_dir)

!KBK-2003-02-10   call init_output(runid,title,start,runtype,dryrun,myid)

   call init_meteo()

#ifndef NO_3D
   call init_rivers()
#endif

   call init_output(runid,title,start,runtype,dryrun,myid)

   call init_2d(runtype,timestep,hotstart)

#ifndef NO_3D
   if (runtype .gt. 1) then
      call init_3d(runtype,timestep,hotstart)
      call init_turbulence(60,trim(input_dir) // 'gotmturb.inp',kmax)
      call init_tridiagonal(kmax)
   end if
#endif

#if 0
   call init_waves(hotstart)
   call init_biology(hotstart)
#endif

   if (hotstart) then
      LEVEL1 'hotstart'
      if (myid .ge. 0) then
         write(buf,'(I3.3)') myid
         buf = '.' // trim(buf) // '.in'
!STDERR buf
!stop 'kbk: initialise'
      else
         buf = '.in'
      end if
      hot_in = trim(input_dir) //'/'// 'restart' // trim(buf)
      call restart_file(READING,trim(hot_in),MinN,runtype)
#ifndef NO_3D
      if (runtype .gt. 1) then
         call start_macro()
         call coordinates(vert_cord,cord_relax)
      end if
#endif
      call depth_update
#ifndef NO_BAROCLINIC
      if (runtype .ge. 3) call do_eqstate()
#endif
      call update_time(MinN)
      call write_time_string()
      LEVEL3 timestr
      MinN = MinN+1
   end if

#ifndef NO_3D
   if (runtype .ge. 2) then
      do j=jjmin-1,jjmax
         do i=iimin-1,iimax
            ssen(i,j)=z(i,j)
            ssun(i,j)=zu(i,j)
            ssvn(i,j)=zv(i,j)
         end do
      end do
   end if
#endif

   call init_input(input_dir,MinN)

   if(runtype .le. 2) then
      call do_meteo(MinN)
#ifndef NO_BAROCLINIC
   else
      call do_meteo(MinN,T(:,:,kmax))
#endif
   end if

   if (save_initial .and. .not. dryrun) call do_output(runtype,0,_ZERO_)

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
