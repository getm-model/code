!$Id: initialise.F90,v 1.1.1.1 2002-05-02 14:01:25 gotm Exp $
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
   public	:: init_model
   integer	:: runtype=1
   logical	:: dryrun=.false.
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: initialise.F90,v $
!  Revision 1.1.1.1  2002-05-02 14:01:25  gotm
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
   subroutine init_model()
!
! !USES:
   use commhalo,     only: init_mpi,print_mpi_info,myid,nprocs
   use output,       only: init_output,do_output,restart_file
   use input,        only: init_input
   use domain,       only: init_domain,iimin,iimax,jjmin,jjmax,kmax
   use domain,       only: vert_cord
   use time,         only: init_time,update_time,write_time_string
   use time,         only: start,timestr,timestep
   use m2d,          only: init_2d,z,zu,zv
   use m3d,          only: cord_relax,init_3d,ssen,ssun,ssvn,T
   use turbulence,   only: init_turbulence
   use mtridiagonal, only: init_tridiagonal
   use meteo,        only: init_meteo,do_meteo
   use rivers,       only: init_rivers
   use integration,  only: MinN,MaxN
   use eqstate,      only: do_eqstate
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
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
   integer		:: i,j
   character(len=64)	:: runid
   character(len=80)	:: title
   logical		:: parallel=.false.
   logical		:: hotstart=.false.
   logical		:: save_initial=.false.
   namelist /param/ dryrun,runid,title,parallel,runtype,	&
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

   STDERR LINE
   STDERR 'Initialising....'
   STDERR LINE
!
! Open the namelist file to get basic run parameters.
!
   title='A descriptive title can be specified in the param namelist'
   open(NAMLST,status='unknown',file='getm.inp')
   read(NAMLST,NML=param)

   LEVEL1 'The run id is: ',trim(runid)
   LEVEL1 'The title is:  ',trim(title)
!
! call all modules Init... routines
!
   if (parallel) then
     call init_mpi()
     call print_mpi_info()
   else
     LEVEL1 'OK - we are making a sequential run'
     myid = -1; nprocs = 1
   end if

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

! This is not at all ready yet - requires change in sequence of namelists!!!!
! Maybe it should not be done - main reason is handling of - hotstart.

   call init_time(MinN,MaxN)

   call init_domain()

   call init_output(runid,title,start,runtype,dryrun)

   call init_meteo()

   call init_rivers()

   call init_2d(runtype,timestep,hotstart)

   if (runtype .gt. 1) then
      call init_3d(runtype,timestep,hotstart)
      call init_turbulence(60,'gotmturb.inp',kmax)
      call init_tridiagonal(kmax)
   end if

#if 0
   call init_waves(hotstart)
   call init_biology(hotstart)
#endif

   if (hotstart) then
      LEVEL1 'hotstart'
      call restart_file(READING,'restart.in',MinN,runtype)
      if (runtype .gt. 1) then
         call start_macro()
         call coordinates(vert_cord,cord_relax)
      end if
      call depth_update
      if (runtype .gt. 1) call do_eqstate()
      call update_time(MinN)
      call write_time_string()
      LEVEL3 timestr
      MinN = MinN+1
   end if

   if (runtype .gt. 1) then
      do j=jjmin-1,jjmax
         do i=iimin-1,iimax
            ssen(i,j)=z(i,j)
            ssun(i,j)=zu(i,j)
            ssvn(i,j)=zv(i,j)
         end do
      end do
   end if

   call init_input(MinN)

   call do_meteo(MinN,T(:,:,kmax))

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
