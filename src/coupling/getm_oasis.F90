#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  getm_oasis - routines to integrate GETM into OASIS
!
! !INTERFACE:
   module getm_oasis
!
! !DESCRIPTION:
!
! !USES:
   use mod_oasis

   IMPLICIT NONE
   private
!
! !PUBLIC DATA MEMBERS:
   public do_getm_oasis
!
! !PRIVATE DATA MEMBERS:
   integer :: compid,il_part_id
!
! !REVISION HISTORY:
!  Original author(s): Xaver Lange & Knut Klingbeil
!
!EOP
!-----------------------------------------------------------------------

    contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: do_getm_oasis
!
! !INTERFACE:
   subroutine do_getm_oasis()
!
! !DESCRIPTION:
!
! !USES:
   use initialise , only: init_model,dryrun,runtype
   use integration, only: time_step,MinN,MaxN
   use getm_timers, only: write_getm_timers,tic,toc,TIM_OASIS
#ifdef GETM_PARALLEL
   use mpi
   use halo_mpi   , only: comm_getm
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES
   character(len=8)                    :: datestr
   character(len=10)                   :: timestr
#ifdef GETM_PARALLEL
   character(len=MPI_MAX_ERROR_STRING) :: mpierrmsg
#endif
   integer                             :: local_comm
   integer                             :: n,length
   integer                             :: ierror,rc
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_getm_oasis() # ',Ncall
#endif

   call tic(TIM_OASIS)

   call oasis_init_comp(compid,"getm",ierror)
   if (ierror /= 0) call oasis_abort(compid,"InitializeP0()","cannot init component")

   call oasis_get_localcomm(local_comm,ierror)
   if (ierror /= 0) call oasis_abort(compid,"InitializeP1()","cannot get localcomm")

#ifdef GETM_PARALLEL
   call MPI_COMM_DUP(local_comm,comm_getm,ierror)
   if (ierror .ne. MPI_SUCCESS) then
!     need depends on specified mpi error handler (i.e. not MPI_ERRORS_ARE_FATAL)
      call MPI_ERROR_STRING(ierror,mpierrmsg,length,rc)
      call oasis_abort(compid,"InitializeP1()","cannot duplicate communicator")
   end if
#endif

   call toc(TIM_OASIS)
   call date_and_time(datestr,timestr)
   call init_model(datestr,timestr)
   call tic(TIM_OASIS)

   call define_grid()
   !call declare_fields()

   call oasis_enddef(ierror)
   if (ierror /= 0) call oasis_abort(compid,"getmComp_init_grid()","cannot end definition phase")

   call toc(TIM_OASIS)

   if (.not.dryrun) then

      STDERR LINE
      LEVEL1 'integrating....'
      STDERR LINE

      do n=MinN,MaxN
         call tic(TIM_OASIS)
!        exchange data...
         call toc(TIM_OASIS)
         call time_step(runtype,n)
      end do

   end if

   call clean_up(dryrun,runtype,MaxN)
#ifndef NO_TIMERS
   STDERR LINE
   call write_getm_timers
#endif

   call tic(TIM_OASIS)

   call oasis_terminate(ierror)

   call toc(TIM_OASIS)

#ifdef DEBUG
   write(debug,*) 'Leaving do_getm_oasis()'
   write(debug,*)
#endif
   return

   end subroutine do_getm_oasis
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: define_grid -
!
! !INTERFACE:
   subroutine define_grid()
!
! !DESCRIPTION:
!  Defines partition and prepares grid data files.
!
! !USES:
   use domain    , only: ioff,imin,imax,iextr,joff,jmin,jmax,jextr
   use domain    , only: lonc,latc,lonx,latx,areac,az,grid_type
   use halo_zones, only: nprocs
#ifdef GETM_PARALLEL
   use mpi
   use halo_mpi  , only: comm_getm
#endif
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES
   REALTYPE,dimension(imin:imax,jmin:jmax,4) :: clon,clat
   integer,dimension(imin:imax,jmin:jmax)    :: mask
   integer                                   :: myedges(4),alledges(4*nprocs)
   integer                                   :: i0,j0,nx_global,ny_global
   integer                                   :: i,j,idx(1),length
   integer                                   :: ierror,rc
#ifdef GETM_PARALLEL
   character(len=MPI_MAX_ERROR_STRING)       :: mpierrmsg
#endif
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'define_grid() # ',Ncall
#endif

   myedges = (/ ioff , joff , imax , jmax /)

#ifdef GETM_PARALLEL
   call MPI_ALLGATHER(myedges,4,MPI_INTEGER,alledges,4,MPI_INTEGER,comm_getm,ierror)
   if (ierror .ne. MPI_SUCCESS) then
!     need depends on specified mpi error handler (i.e. not MPI_ERRORS_ARE_FATAL)
      call MPI_ERROR_STRING(ierror,mpierrmsg,length,rc)
      call oasis_abort(compid,"getmComp_init_grid()","cannot mpi_allgather")
   end if
#else
   alledges = myedges
#endif

   i0 = minval(alledges(1::4))
   j0 = minval(alledges(2::4))

   idx = maxloc(alledges(1::4))
   nx_global = alledges(1+4*(idx(1)-1)) + alledges(3+4*(idx(1)-1)) - i0

   idx = maxloc(alledges(2::4))
   ny_global = alledges(2+4*(idx(1)-1)) + alledges(4+4*(idx(1)-1)) - j0

!  Box partition
   call oasis_def_partition(il_part_id,                                              &
                            (/2,(joff-j0)*nx_global+(ioff-i0),imax,jmax,nx_global/), &
                            ierror,iextr*jextr,"getm")
   if (ierror /= 0) call oasis_abort(compid,"getmComp_init_grid()","cannot define partition")

   if (grid_type.ne.2 .and. grid_type.ne.4) then
!      call oasis_abort(compid,"getmComp_init_grid()","grid_type not supported")
   end if

   do i=imax,imax
      do j=jmin,jmax
         clon(i,j,1) = lonx(i  ,j  )
         clon(i,j,2) = lonx(i-1,j  )
         clon(i,j,3) = lonx(i-1,j  )
         clon(i,j,4) = lonx(i  ,j  )
         clat(i,j,1) = latx(i  ,j  )
         clat(i,j,2) = latx(i  ,j  )
         clat(i,j,3) = latx(i  ,j-1)
         clat(i,j,4) = latx(i  ,j-1)
      end do
   end do
   mask = ( az(imin:imax,jmin:jmax) .eq. 0 )

   call oasis_start_grids_writing(ierror)
   call oasis_write_grid  ("getm",nx_global,ny_global,lonc(imin:imax,jmin:jmax),latc(imin:imax,jmin:jmax),il_part_id)
   call oasis_write_corner("getm",nx_global,ny_global,4,clon,clat,il_part_id)
   call oasis_write_mask  ("getm",nx_global,ny_global,mask,il_part_id)
   call oasis_write_area  ("getm",nx_global,ny_global,areac(imin:imax,jmin:jmax),il_part_id)
   call oasis_terminate_grids_writing()

#ifdef DEBUG
   write(debug,*) 'Leaving define_grid()'
   write(debug,*)
#endif
   return

   end subroutine define_grid
!EOC
!-----------------------------------------------------------------------
#if 0
!BOP
!
! !ROUTINE: InitializeP1 -
!
! !INTERFACE:
   subroutine InitializeP1()
!
! !DESCRIPTION:
!
! !USES:
   use time       , only: init_time,start,timestep
   use initialise , only: init_model,init_initialise,do_initialise,dryrun
   use integration, only: MinN,MaxN
   use getm_timers, only: 
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES
   integer                    :: localrc,comm,length,getmRunTimeStepCount
   logical                    :: abort,clockIsPresent
   character(len=ESMF_MAXSTR) :: name="getm"
   character(len=19)          :: TimeStrISOFrac,start_external,stop_external
!
!EOP
!-----------------------------------------------------------------------
!BOC
   return

   end subroutine InitializeP1
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: InitializeP2 -
!
! !INTERFACE:
   subroutine InitializeP2(getmComp,importState,exportState,clock,rc)
!
! !DESCRIPTION:
!  Here neccessary import and export fields are allocated.
!
!  Note: states and clock are uninitialized if the toplevel component
!        did not provide corresponding arguments to
!        ESMF_GridCompInitialize(getmComp).
!        Status can be checked with [State|Clock]IsCreated().
!  The toplevel component can inquire rc via optional keyword argument
!  userRc to ESMF_GridCompInitialize().
!
! !USES:
   use initialise, only: dryrun
   use getm_timers, only: tic,toc,TIM_OASIS
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp) :: getmComp
   type(ESMF_State)    :: importState,exportState ! may be uninitialized
   type(ESMF_Clock)    :: clock                   ! may be uninitialized
!
! !OUTPUT PARAMETERS:
   integer,intent(out) :: rc
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES
   logical                    :: abort
   character(len=ESMF_MAXSTR) :: name="getm"
!
!EOP
!-----------------------------------------------------------------------
!BOC
   return

   end subroutine InitializeP2
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: RunP1 -
!
! !INTERFACE:
   subroutine RunP1(getmComp,importState,exportState,clock,rc)
!
! !DESCRIPTION:
!  Note: states and clock are uninitialized if the toplevel component
!        did not provide corresponding arguments to
!        ESMF_GridCompRun(getmComp).
!        Status can be checked with [State|Clock]IsCreated().
!  The toplevel component can inquire rc via optional keyword argument
!  userRc to ESMF_GridCompRun().
!
! !USES:


   use getm_timers, only: tic,toc,TIM_OASIS
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp) :: getmComp
   type(ESMF_State)    :: importState,exportState ! may be uninitialized
   type(ESMF_Clock)    :: clock                   ! may be uninitialized
!
! !OUTPUT PARAMETERS:
   integer,intent(out) :: rc
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES
   type(ESMF_Clock)           :: getmClock
   type(ESMF_Time)            :: getmTime,nextTime
   type(ESMF_TimeInterval)    :: getmTimeStep
   integer(ESMF_KIND_I8)      :: advanceCount
   integer(kind=kind(MinN))   :: n
   logical                    :: abort
   character(len=ESMF_MAXSTR) :: name="getm"
!
!EOP
!-----------------------------------------------------------------------
!BOC
   return

   end subroutine RunP1
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: FinalizeP1 -
!
! !INTERFACE:
   subroutine FinalizeP1(getmComp,importState,exportState,clock,rc)
!
! !DESCRIPTION:
!  Note: states and clock are uninitialized if the toplevel component
!        did not provide corresponding arguments to
!        ESMF_GridCompFinalize(getmComp).
!        Status can be checked with [State|Clock]IsCreated().
!  The toplevel component can inquire rc via optional keyword argument
!  userRc to ESMF_GridCompFinalize().
!
! !USES:
   use initialise , only: runtype,dryrun
   use integration, only: MaxN

   use getm_timers, only: tic,toc,TIM_OASIS
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp) :: getmComp
   type(ESMF_State)    :: importState,exportState ! may be uninitialized
   type(ESMF_Clock)    :: clock                   ! may be uninitialized
!
! !OUTPUT PARAMETERS:
   integer,intent(out) :: rc
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES
   type(ESMF_Clock)           :: getmClock
   type(ESMF_Grid)            :: getmGrid
   type(ESMF_DistGrid)        :: getmDistGrid3D
   type(type_getmInternalState)               :: gis
   type(type_getmInternalStateStruct),pointer :: isd
   logical                    :: abort,clockIsPresent,GridIsPresent
   character(len=ESMF_MAXSTR) :: name="getm"
!
!EOP
!-----------------------------------------------------------------------
!BOC
   return

   end subroutine FinalizeP1
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: getmComp_init_grid - Creates Grid
!
! !INTERFACE:
   subroutine getmComp_init_grid(getmComp,getmGrid2D)
!
! !DESCRIPTION:
!
! !USES:
   use initialise, only: runtype
   use domain    , only: ioff,joff,imin,jmin,imax,jmax,kmax
   use domain    , only: grid_type,have_lonlat
   use domain    , only: az,ax,arcd1
   use domain    , only: xcord,ycord,xc,yc,lonc,latc
   use domain    , only: xxcord,yxcord,xx,yx,lonx,latx
#ifndef NO_3D
   use variables_3d, only: zwn,zcn
#endif
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp) :: getmComp
!
! !OUTPUT PARAMETERS:
   type(ESMF_Grid),intent(out) :: getmGrid2D
!
! !REVISION HISTORY:
!  Original Author(s): Knut Klingbeil
!
! !LOCAL VARIABLES
   real(ESMF_KIND_R8),pointer :: xc1D  (:)    =>NULL()
   real(ESMF_KIND_R8),pointer :: yc1D  (:)    =>NULL()
   real(ESMF_KIND_R8),pointer :: xx1D  (:)    =>NULL()
   real(ESMF_KIND_R8),pointer :: yx1D  (:)    =>NULL()
   real(ESMF_KIND_R8),pointer :: xc2D  (:,:)  =>NULL()
   real(ESMF_KIND_R8),pointer :: yc2D  (:,:)  =>NULL()
   real(ESMF_KIND_R8),pointer :: xx2D  (:,:)  =>NULL()
   real(ESMF_KIND_R8),pointer :: yx2D  (:,:)  =>NULL()
   real(ESMF_KIND_R8),pointer :: lonc1D(:)    =>NULL()
   real(ESMF_KIND_R8),pointer :: latc1D(:)    =>NULL()
   real(ESMF_KIND_R8),pointer :: lonx1D(:)    =>NULL()
   real(ESMF_KIND_R8),pointer :: latx1D(:)    =>NULL()
   real(ESMF_KIND_R8),pointer :: lonc2D(:,:)  =>NULL()
   real(ESMF_KIND_R8),pointer :: latc2D(:,:)  =>NULL()
   real(ESMF_KIND_R8),pointer :: lonx2D(:,:)  =>NULL()
   real(ESMF_KIND_R8),pointer :: latx2D(:,:)  =>NULL()
   real(ESMF_KIND_R8),pointer :: zw    (:,:,:)=>NULL()
   real(ESMF_KIND_R8),pointer :: zc    (:,:,:)=>NULL()
   real(ESMF_KIND_R8),pointer :: zx    (:,:,:)=>NULL()
   real(ESMF_KIND_R8),pointer :: areaC (:,:)
   real(ESMF_KIND_R8),pointer :: p2dr  (:,:)
   real(ESMF_KIND_R8),dimension(:,:,:),allocatable,target :: areaW3D
   integer(ESMF_KIND_I4),pointer :: maskC(:,:)
   integer(ESMF_KIND_I4),pointer :: maskX(:,:)
   integer(ESMF_KIND_I4),dimension(:,:,:),allocatable,target :: maskC3D
   integer(ESMF_KIND_I4),dimension(:,:,:),allocatable,target :: maskX3D
   integer(ESMF_KIND_I4),pointer                         :: p2di(:,:)
   REALTYPE                 :: getmreal
   integer                  :: rc,pet,i0,j0,ilen,jlen,klen,k
   integer                  :: petCount
   integer,dimension(3)     :: coordDimCount
   integer,dimension(3,3)   :: coordDimMap
   integer,dimension(:,:,:),allocatable                  :: deBlockList
   logical                    :: abort,noKindMatch
!
!EOP
!-----------------------------------------------------------------------
   return

   end subroutine getmComp_init_grid
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: getmComp_update_grid -
!
! !INTERFACE:
   subroutine getmComp_update_grid(getmComp)
!
! !DESCRIPTION:
!
! !USES:
   use initialise  , only: runtype
   use domain      , only: imin,imax,jmin,jmax,kmax,az,au,av,H
   use variables_2d, only: z
#ifndef NO_3D
   use variables_3d, only: zwn,zcn
#endif
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp) :: getmComp
!
! !REVISION HISTORY:
!  Original Author(s): Knut Klingbeil
!
! !LOCAL VARIABLES
   type(ESMF_Grid)                             :: getmGrid3D
   real(ESMF_KIND_R8),dimension(:,:,:),pointer :: zw,zc,zx
   REALTYPE,dimension(E2DFIELD)                :: zwu
   REALTYPE                                    :: getmreal
   integer                                     :: rc,i,j,k,klen
   logical                                     :: abort,noKindMatch
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'getmComp_update_grid() # ',Ncall
#endif

   call ESMF_GridCompGet(getmComp,grid=getmGrid3D,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridGetCoord(getmGrid3D,3,farrayPtr=zw,                   &
                          staggerloc=ESMF_StaggerLoc_CENTER_VFACE,     &
                          rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridGetCoord(getmGrid3D,3,farrayPtr=zc,                   &
                          staggerloc=ESMF_StaggerLoc_CENTER_VCENTER,   &
                          rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridGetCoord(getmGrid3D,3,farrayPtr=zx,                   &
                          staggerloc=ESMF_StaggerLoc_CORNER,           &
                          rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   if (runtype .eq. 1) then
      klen = 1
   else
      klen = kmax
   end if

   noKindMatch = ( kind(getmreal) .ne. ESMF_KIND_R8 )

   if (noKindMatch .or. runtype.eq.1) then
      if (runtype .eq. 1) then
         zw(:,:,0)    = -H
         zw(:,:,klen) = z
         zc(:,:,klen) = 0.5d0 * ( zw(:,:,klen-1) + zw(:,:,klen) )
      else
#ifndef NO_3D
         zw = zwn
         zc = zcn
#endif
      end if
   end if

   do k=0,klen
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO-1
            if (au(i,j) .gt. 0) then
               zwu(i,j) = 0.5d0 * ( zw(i,j,k) + zw(i+1,j,k) )
            else
               if (az(i,j) .gt. 0) then
                  zwu(i,j) = zw(i,j,k)
               else if (az(i+1,j) .gt. 0) then
                  zwu(i,j) = zw(i+1,j,k)
               end if
            end if
         end do
      end do
      do j=jmin-HALO,jmax+HALO-1
         do i=imin-HALO,imax+HALO-1
            if (av(i,j).gt.0 .or. av(i+1,j).gt.0) then
               zx(i,j,k) = 0.5d0 * ( zwu(i,j) + zwu(i,j+1) )
            else
               if (az(i,j).gt.0 .or. az(i+1,j).gt.0) then
                  zx(i,j,k) = zwu(i,j)
               else if (az(i,j+1).gt.0 .or. az(i+1,j+1).gt.0) then
                  zx(i,j,k) = zwu(i,j+1)
               end if
            end if
         end do
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving getmComp_update_grid()'
   write(debug,*)
#endif
   return

   end subroutine getmComp_update_grid
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_importStateP1 -
!
! !INTERFACE:
   subroutine init_importStateP1(getmComp,getmGrid2D,importState)
!
! !DESCRIPTION:
!
! !USES:
   use meteo,only: met_method,calc_met,METEO_FROMEXT
   use waves,only: waveforcing_method,WAVES_FROMEXT
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type(ESMF_Grid),intent(in)   :: getmGrid2D
!
! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp) :: getmComp
   type(ESMF_State)    :: importState
!
! !REVISION HISTORY:
!  Original Author(s): Knut Klingbeil
!
! !LOCAL VARIABLES
   type(ESMF_Grid) :: getmGrid3D
   integer         :: rc
   logical         :: abort
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_importStateP1() # ',Ncall
#endif

   call ESMF_GridCompGet(getmComp,grid=getmGrid3D,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call StateAddGridSetField(importState,"gridSetField2D",getmGrid2D)

   if (met_method .eq. METEO_FROMEXT) then
   !call StateAddGridSetField(importState,trim(name_swr    ),getmGrid2D,units="W m-2")
   if (calc_met) then
   call StateAddGridSetField(importState,trim(name_windU  ),getmGrid2D,units="m s-1")
   call StateAddGridSetField(importState,trim(name_windV  ),getmGrid2D,units="m s-1")
   end if ! calc_met
   end if ! meteo

   if (waveforcing_method .eq. WAVES_FROMEXT) then
   call StateAddGridSetField(importState,trim(name_waveDir),getmGrid2D,units="rad")
   call StateAddGridSetField(importState,trim(name_waveH  ),getmGrid2D,units="m")
   call StateAddGridSetField(importState,trim(name_waveK  ),getmGrid2D,units="m-1")
   call StateAddGridSetField(importState,trim(name_waveT  ),getmGrid2D,units="s")
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving init_importStateP1()'
   write(debug,*)
#endif
   return

   end subroutine init_importStateP1
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: read_importState -
!
! !INTERFACE:
   subroutine read_importState(getmComp,importState)
!
! !DESCRIPTION:
!
! !USES:
   use domain         ,only: grid_type,convc,cosconv,sinconv
   use meteo          ,only: met_method,calc_met,METEO_FROMEXT,new_meteo
   use meteo          ,only: u10,v10
   use waves          ,only: waveforcing_method,WAVES_FROMEXT,new_waves
   use variables_waves,only: coswavedir,sinwavedir,waveH,waveK,waveT
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp) :: getmComp
   type(ESMF_State)    :: importState
!
! !REVISION HISTORY:
!  Original Author(s): Knut Klingbeil
!
! !LOCAL VARIABLES
   real(ESMF_KIND_R8),pointer :: p2dr(:,:),windU(:,:),windV(:,:)
   logical                    :: frc
   REALTYPE, parameter :: pi=3.1415926535897932384626433832795029d0
   REALTYPE, parameter :: deg2rad=pi/180
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'read_importState() # ',Ncall
#endif

   if (met_method .eq. METEO_FROMEXT) then

      new_meteo = .true. ! KK-TODO: should be set by coupler
      !call StateCompleteConnectedField(importState,trim(name_swr  ),farray2d=swr)

      if (calc_met) then
!        force reading if grid rotation needs to be removed
         frc = (grid_type .ne. 2)
         call StateReadCompleteField(importState,trim(name_windU  ),   &
                                     farray2d=u10,p2dr=windU,frc=frc)
         call StateReadCompleteField(importState,trim(name_windV  ),   &
                                     farray2d=v10,p2dr=windV,frc=frc)
         if (frc) then
            u10 = cosconv*windU - sinconv*windV
            v10 = sinconv*windU + cosconv*windV
         end if
      end if ! calc_met

   end if ! meteo


   if (waveforcing_method .eq. WAVES_FROMEXT) then

         new_waves = .true. ! KK-TODO: should be set by coupler

         call StateReadCompleteField(importState,trim(name_waveDir),   &
                                     p2dr=p2dr,frc=.true.)
         coswavedir = cos( p2dr + convc*deg2rad )
         sinwavedir = sin( p2dr + convc*deg2rad )

         call StateReadCompleteField(importState,trim(name_waveH  ),   &
                                     farray2d=waveH)
         call StateReadCompleteField(importState,trim(name_waveK  ),   &
                                     farray2d=waveK)
         call StateReadCompleteField(importState,trim(name_waveT  ),   &
                                     farray2d=waveT)

   end if

#ifdef DEBUG
   write(debug,*) 'Leaving read_importState()'
   write(debug,*)
#endif
   return

   end subroutine read_importState
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_exportStateP1 -
!
! !INTERFACE:
   subroutine init_exportStateP1(getmComp,getmGrid2D,exportState)
!
! !DESCRIPTION:
!
! !USES:
   use initialise     ,only: runtype
#ifndef NO_3D
   use m3d            ,only: calc_salt,calc_temp
#endif
   use meteo          ,only: met_method,METEO_CONST,METEO_FROMFILE
   use waves          ,only: waveforcing_method
   use waves          ,only: WAVES_FROMWIND,WAVES_FROMFILE
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type(ESMF_Grid),intent(in)   :: getmGrid2D
!
! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp) :: getmComp
   type(ESMF_State)    :: exportState
!
! !REVISION HISTORY:
!  Original Author(s): Knut Klingbeil
!
! !LOCAL VARIABLES
   type(ESMF_Grid) :: getmGrid3D
   integer         :: rc
   logical         :: abort
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_exportStateP1() # ',Ncall
#endif

   call ESMF_GridCompGet(getmComp,grid=getmGrid3D,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call StateAddGridSetField(exportState,"gridSetField2D",getmGrid2D)

#ifdef ALLEXPORT
   call StateAddGridSetField(exportState,trim(name_depth  ),getmGrid2D,units="m")
   call StateAddGridSetField(exportState,trim(name_h3D    ),getmGrid3D,units="m")
   call StateAddGridSetField(exportState,trim(name_hbot   ),getmGrid2D,units="m")
   call StateAddGridSetField(exportState,trim(name_U2D    ),getmGrid2D,units="m s-1")
   call StateAddGridSetField(exportState,trim(name_U3D    ),getmGrid3D,units="m s-1")
   call StateAddGridSetField(exportState,trim(name_Ubot   ),getmGrid2D,units="m s-1")
   call StateAddGridSetField(exportState,trim(name_V2D    ),getmGrid2D,units="m s-1")
   call StateAddGridSetField(exportState,trim(name_V3D    ),getmGrid3D,units="m s-1")
   call StateAddGridSetField(exportState,trim(name_Vbot   ),getmGrid2D,units="m s-1")
   if (runtype .ge. 2) then
#ifndef NO_3D
   call StateAddGridSetField(exportState,trim(name_eps3D  ),getmGrid3D,units="m2 s-3",staggerloc=ESMF_STAGGERLOC_CENTER_VFACE)
   call StateAddGridSetField(exportState,trim(name_epsbot ),getmGrid2D,units="m2 s-3")
   call StateAddGridSetField(exportState,trim(name_num3D  ),getmGrid3D,units="m2 s-1",staggerloc=ESMF_STAGGERLOC_CENTER_VFACE)
   call StateAddGridSetField(exportState,trim(name_numbot ),getmGrid2D,units="m2 s-1")
   call StateAddGridSetField(exportState,trim(name_SS3D   ),getmGrid3D,units="s-2",staggerloc=ESMF_STAGGERLOC_CENTER_VFACE)
   call StateAddGridSetField(exportState,trim(name_taubmax),getmGrid2D,units="Pa")
   call StateAddGridSetField(exportState,trim(name_tke3D  ),getmGrid3D,units="m2 s-2",staggerloc=ESMF_STAGGERLOC_CENTER_VFACE)
   call StateAddGridSetField(exportState,trim(name_tkebot ),getmGrid2D,units="m2 s-2")
   if (runtype .ge. 3) then
#ifndef NO_BAROCLINIC
   call StateAddGridSetField(exportState,trim(name_NN3D   ),getmGrid3D,units="s-2",staggerloc=ESMF_STAGGERLOC_CENTER_VFACE)
   call StateAddGridSetField(exportState,trim(name_nuh3D  ),getmGrid3D,units="m",staggerloc=ESMF_STAGGERLOC_CENTER_VFACE)
   if (calc_salt) then
   call StateAddGridSetField(exportState,trim(name_S3D    ),getmGrid3D,units="")
   end if
   if (calc_temp) then
   call StateAddGridSetField(exportState,trim(name_T3D    ),getmGrid3D,units="degC")
   call StateAddGridSetField(exportState,trim(name_Tbot   ),getmGrid2D,units="degC")
   end if
#endif
   end if
#endif
   end if
   if (met_method.eq.METEO_CONST .or. met_method.eq.METEO_FROMFILE) then
   call StateAddGridSetField(exportState,trim(name_swr    ),getmGrid2D,units="W m-2")
   call StateAddGridSetField(exportState,trim(name_windU  ),getmGrid2D,units="m s-1")
   call StateAddGridSetField(exportState,trim(name_windV  ),getmGrid2D,units="m s-1")
   end if
   if (     waveforcing_method.eq.WAVES_FROMWIND                       &
       .or. waveforcing_method.eq.WAVES_FROMFILE) then
   call StateAddGridSetField(exportState,trim(name_waveDir),getmGrid2D,units="rad")
   call StateAddGridSetField(exportState,trim(name_waveH  ),getmGrid2D,units="m")
   call StateAddGridSetField(exportState,trim(name_waveK  ),getmGrid2D,units="m-1")
   call StateAddGridSetField(exportState,trim(name_waveT  ),getmGrid2D,units="s")
   end if
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving init_exportStateP1()'
   write(debug,*)
#endif
   return

   end subroutine init_exportStateP1
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: update_exportState -
!
! !INTERFACE:
   subroutine update_exportState(getmComp,exportState)
!
! !DESCRIPTION:
!
! !USES:
   use initialise     ,only: runtype
   use domain         ,only: kmax
   use domain         ,only: grid_type,convc,cosconv,sinconv
   use parameters     ,only: rho_0
   use variables_2d   ,only: D,Dvel
   use variables_2d   ,only: velx,vely
#ifndef NO_3D
   use m3d            ,only: calc_salt,calc_temp
   use variables_3d   ,only: Dn,Dveln,hn,hvel
   use variables_3d   ,only: velx3d,vely3d,velx2dadv,vely2dadv
   use variables_3d   ,only: taubmax_3d
   use variables_3d   ,only: eps,nuh,num,SS,tke
#ifndef NO_BAROCLINIC
   use variables_3d   ,only: NN,S,T
#endif
#endif
   use meteo          ,only: met_method,METEO_CONST,METEO_FROMFILE
   use meteo          ,only: swr,u10,v10
   use waves          ,only: waveforcing_method
   use waves          ,only: NO_WAVES,WAVES_FROMWIND,WAVES_FROMFILE
   use waves          ,only: waves_method,WAVES_NOSTOKES
   use variables_waves,only: coswavedir,sinwavedir,waveH,waveK,waveT
   use variables_waves,only: UStokesC,VStokesC
   use variables_waves,only: UStokesCadv,VStokesCadv,uuStokesC,vvStokesC
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp) :: getmComp
   type(ESMF_State)    :: exportState
!
! !REVISION HISTORY:
!  Original Author(s): Knut Klingbeil
!
! !LOCAL VARIABLES
   real(ESMF_KIND_R8),pointer :: p2dr(:,:),p3dr(:,:,:)
   integer                    :: k
   logical                    :: frc
   REALTYPE, parameter :: pi=3.1415926535897932384626433832795029d0
   REALTYPE, parameter :: deg2rad=pi/180
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'update_exportState() # ',Ncall
#endif

#ifdef ALLEXPORT
!  force allocation of new memory if Stokes drift needs to be removed
   frc = (waveforcing_method.ne.NO_WAVES .and. waves_method.ne.WAVES_NOSTOKES)

   if (runtype .eq. 1) then

!     2D-only and general 3D variables

         p2dr => NULL()
         call StateSetCompleteField(exportState,trim(name_depth  ),    &
                                    farray2d=D,p2dr=p2dr,frc=.true.)
         p3dr => NULL()
         call StateReadCompleteField(exportState,trim(name_h3D    ),   &
                                     p3dr=p3dr,frc=.true.,ign=.true.)
      if (associated(p3dr)) then
         p3dr(:,:,1) = D
!        if hbot present, then already pointing to h3D(:,:,1)
      else if (associated(p2dr)) then
!        if hbot present, then already pointing to depth
      else
         call StateSetCompleteField(exportState,trim(name_depth  ),    &
                                    farray2d=D)
      end if

         p2dr => NULL()
         call StateSetCompleteField(exportState,trim(name_U2D    ),    &
                                    farray2d=velx,p2dr=p2dr,frc=.true.)
      if (associated(p2dr)) then
         if (frc) p2dr = velx - (  cosconv*UStokesC + sinconv*VStokesC )/Dvel
      end if

         p3dr => NULL()
         call StateReadCompleteField(exportState,trim(name_U3D    ),   &
                                     p3dr=p3dr,frc=.true.,ign=.true.)
      if (associated(p3dr)) then
         if (frc) then
            p3dr(:,:,1) = velx - (  cosconv*UStokesC + sinconv*VStokesC )/Dvel
         else
            p3dr(:,:,1) = velx
         end if
!        if Ubot present, then already pointing to U3D(:,:,1)
      else if (associated(p2dr)) then
!        if Ubot present, then already pointing to U2D
      else
         p2dr => NULL()
         call StateSetCompleteField(exportState,trim(name_Ubot   ),    &
                                    farray2d=velx,p2dr=p2dr,frc=.true.)
      if (associated(p2dr)) then
         if (frc) p2dr = velx - (  cosconv*UStokesC + sinconv*VStokesC )/Dvel
      end if
      end if

         p2dr => NULL()
         call StateSetCompleteField(exportState,trim(name_V2D    ),    &
                                    farray2d=vely,p2dr=p2dr,frc=.true.)
      if (associated(p2dr)) then
         if (frc) p2dr = vely - ( -sinconv*UStokesC + cosconv*VStokesC )/Dvel
      end if

         p3dr => NULL()
         call StateReadCompleteField(exportState,trim(name_V3D    ),   &
                                     p3dr=p3dr,frc=.true.,ign=.true.)
      if (associated(p3dr)) then
         if (frc) then
            p3dr(:,:,1) = vely - ( -sinconv*UStokesC + cosconv*VStokesC )/Dvel
         else
            p3dr(:,:,1) = vely
         end if
!        if Vbot present, then already pointing to V3D(:,:,1)
      else if (associated(p2dr)) then
!        if Vbot present, then already pointing to V2D
      else
         p2dr => NULL()
         call StateSetCompleteField(exportState,trim(name_Vbot   ),    &
                                    farray2d=vely,p2dr=p2dr,frc=.true.)
      if (associated(p2dr)) then
         if (frc) p2dr = vely - ( -sinconv*UStokesC + cosconv*VStokesC )/Dvel
      end if
      end if

   else ! runtype.gt.1

#ifndef NO_3D

         call StateSetCompleteField(exportState,trim(name_depth  ),    &
                                    farray2d=Dn)
         p3dr => NULL()
         call StateSetCompleteField(exportState,trim(name_h3D    ),    &
                                    farray3D=hn,p3dr=p3dr,frc=.true.)
      if (associated(p3dr)) then
!        if hbot present, then already pointing to h3D(:,:,1)
      else
         call StateSetCompleteField(exportState,trim(name_hbot   ),    &
                                    farray2D=hn(:,:,1))
      end if

         p2dr => NULL()
         call StateSetCompleteField(exportState,trim(name_U2D    ),    &
                                    farray2D=velx2dadv,p2dr=p2dr,frc=frc)
      if (associated(p2dr)) then
         if (frc) p2dr = velx2dadv - (  cosconv*UStokesCadv + sinconv*VStokesCadv )/Dveln
      end if

         p3dr => NULL()
         call StateSetCompleteField(exportState,trim(name_U3D    ),    &
                                    farray3D=velx3d,p3dr=p3dr,frc=frc)
      if (associated(p3dr)) then
         if (frc) then
            do k=1,kmax
            p3dr(:,:,k) = velx3d(:,:,k) - (  cosconv*uuStokesC(:,:,k) + sinconv*vvStokesC(:,:,k) )/hvel(:,:,k)
            end do
          end if
!        if Ubot present, then already pointing to U3D(:,:,1)
      else
         p2dr => NULL()
         call StateSetCompleteField(exportState,trim(name_Ubot   ),    &
                                    farray2D=velx3d(:,:,1),p2dr=p2dr,frc=frc)
         if (associated(p2dr)) then
         if (frc) p2dr = velx3d(:,:,1) - (  cosconv*uuStokesC(:,:,1) + sinconv*vvStokesC(:,:,1) )/hvel(:,:,1)
         end if
      end if

         p2dr => NULL()
         call StateSetCompleteField(exportState,trim(name_V2D    ),    &
                                    farray2D=vely2dadv,p2dr=p2dr,frc=frc)
         if (associated(p2dr)) then
         if (frc) p2dr = vely2dadv - ( -sinconv*UStokesCadv + cosconv*VStokesCadv )/Dveln
         end if

         p3dr => NULL()
         call StateSetCompleteField(exportState,trim(name_V3D    ),    &
                                    farray3D=vely3d,p3dr=p3dr,frc=frc)
      if (associated(p3dr)) then
         if (frc) then
            do k=1,kmax
            p3dr(:,:,k) = vely3d(:,:,k) - ( -sinconv*uuStokesC(:,:,k) + cosconv*vvStokesC(:,:,k) )/hvel(:,:,k)
            end do
         end if
!        if Vbot present, then already pointing to V3D(:,:,1)
      else
         p2dr => NULL()
         call StateSetCompleteField(exportState,trim(name_Vbot   ),    &
                                    farray2D=vely3d(:,:,1),p2dr=p2dr,frc=frc)
         if (associated(p2dr)) then
         if (frc) p2dr = vely3d(:,:,1) - ( -sinconv*uuStokesC(:,:,1) + cosconv*vvStokesC(:,:,1) )/hvel(:,:,1)
         end if
      end if

!     additional 3D-only variables

         p3dr => NULL()
         call StateReadCompleteField(exportState,trim(name_eps3D  ),   &
                                     p3dr=p3dr,frc=.true.,ign=.true.)
      if (associated(p3dr)) then
         p3dr = eps
!        if epsbot present, then already pointing to eps3D(:,:,1)
      else
         p2dr => NULL()
         call StateReadCompleteField(exportState,trim(name_epsbot ),   &
                                     p2dr=p2dr,frc=.true.,ign=.true.)
         if (associated(p2dr)) p2dr = eps(:,:,1)
      end if

         p3dr => NULL()
         call StateReadCompleteField(exportState,trim(name_num3D  ),   &
                                     p3dr=p3dr,frc=.true.,ign=.true.)
      if (associated(p3dr)) then
         p3dr = num
!        if numbot present, then already pointing to num3D(:,:,1)
      else
         p2dr => NULL()
         call StateReadCompleteField(exportState,trim(name_numbot ),   &
                                     p2dr=p2dr,frc=.true.,ign=.true.)
         if (associated(p2dr)) p2dr = num(:,:,1)
      end if

         p3dr => NULL()
         call StateReadCompleteField(exportState,trim(name_SS3D   ),   &
                                     p3dr=p3dr,frc=.true.,ign=.true.)
         if (associated(p3dr)) p3dr = SS

         p2dr => NULL()
         call StateReadCompleteField(exportState,trim(name_taubmax),   &
                                     p2dr=p2dr,frc=.true.,ign=.true.)
         if (associated(p2dr)) p2dr = rho_0 * taubmax_3d

         p3dr => NULL()
         call StateReadCompleteField(exportState,trim(name_tke3D  ),   &
                                     p3dr=p3dr,frc=.true.,ign=.true.)
      if (associated(p3dr)) then
         p3dr = tke
!        if tkebot present, then already pointing to tke3D(:,:,1)
      else
         p2dr => NULL()
         call StateReadCompleteField(exportState,trim(name_tkebot ),   &
                                     p2dr=p2dr,frc=.true.,ign=.true.)
         if (associated(p2dr)) p2dr = tke(:,:,1)
      end if

      if (runtype .ge. 3) then

#ifndef NO_BAROCLINIC

         p3dr => NULL()
         call StateReadCompleteField(exportState,trim(name_NN3D   ),   &
                                     p3dr=p3dr,frc=.true.,ign=.true.)
         if (associated(p3dr)) p3dr = NN

         call StateSetCompleteField(exportState,trim(name_nuh3D  ),    &
                                    farray3D=nuh)

      if (calc_salt) then
         call StateSetCompleteField(exportState,trim(name_S3D    ),    &
                                    farray3D=S)
      end if ! calc_salt

      if (calc_temp) then
         p3dr => NULL()
         call StateSetCompleteField(exportState,trim(name_T3D    ),    &
                                    farray3D=T,p3dr=p3dr,frc=.true.)
      if (associated(p3dr)) then
!        if Tbot present, then already pointing to T3D(:,:,1)
      else
         call StateSetCompleteField(exportState,trim(name_Tbot   ),    &
                                    farray2d=T(:,:,1))
      end if
      end if ! calc_temp

      end if ! runtype.ge.3

#endif

#endif

   end if ! runtype


   if (met_method.eq.METEO_CONST .or. met_method.eq.METEO_FROMFILE) then
!     force update of pointer because of pointer swap within GETM
         call StateSetCompleteField(exportState,trim(name_swr    ),    &
                                    farray2D=swr,frc=.true.)
!     check whether grid rotation needs to be removed
         frc = (grid_type .ne. 2)
         p2dr => NULL()
         call StateSetCompleteField(exportState,trim(name_windU  ),    &
                                    farray2D=u10,p2dr=p2dr,frc=frc)
      if (associated(p2dr)) then
         if (frc) p2dr =  cosconv*u10 + sinconv*v10
      end if

         p2dr => NULL()
         call StateSetCompleteField(exportState,trim(name_windV  ),    &
                                    farray2D=v10,p2dr=p2dr,frc=frc)
      if (associated(p2dr)) then
         if (frc) p2dr = -sinconv*u10 + cosconv*v10
      end if

   end if ! meteo

   if (     waveforcing_method.eq.WAVES_FROMWIND                       &
       .or. waveforcing_method.eq.WAVES_FROMFILE) then
         p2dr => NULL()
         call StateReadCompleteField(exportState,trim(name_waveDir),   &
                                     p2dr=p2dr,frc=.true.,ign=.true.)
         if (associated(p2dr)) p2dr = atan2(sinwavedir,coswavedir) - convc*deg2rad

!      force update of pointer because of pointer swap within GETM
         call StateSetCompleteField(exportState,trim(name_waveH  ),    &
                                    farray2D=waveH,frc=.true.)
         call StateSetCompleteField(exportState,trim(name_waveK  ),    &
                                    farray2D=waveK)
         call StateSetCompleteField(exportState,trim(name_waveT  ),    &
                                    farray2D=waveT)
   end if ! waves
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving update_exportState()'
   write(debug,*)
#endif
   return

   end subroutine update_exportState
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: TimeStringISOFrac2ESMFtime - converts timestring to ESMF_Time
!
! !INTERFACE:
   subroutine TimeStringISOFrac2ESMFtime(TimeStrISOFrac,ESMFtime)
!
! !DESCRIPTION:
!  So far missing extension to ESMF_TimeSet().
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)  :: TimeStrISOFrac
!
! !OUTPUT PARAMETERS:
   type(ESMF_Time) , intent(out) :: ESMFtime
!
! !REVISION HISTORY:
!  Original Author(s): Carsten Lemmen and Richard Hofmeister
!
! !LOCAL VARIABLES
   integer :: yy,mm,dd,h,m,s,rc
   logical :: abort
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'TimeStringISOFrac2ESMFtime() # ',Ncall
#endif

   read(TimeStrISOFrac( 1: 4),'(i4)') yy
   read(TimeStrISOFrac( 6: 7),'(i2)') mm
   read(TimeStrISOFrac( 9:10),'(i2)') dd
   read(TimeStrISOFrac(12:13),'(i2)') h
   read(TimeStrISOFrac(15:16),'(i2)') m
   read(TimeStrISOFrac(18:19),'(i2)') s

   call ESMF_TimeSet(ESMFtime,yy=yy,mm=mm,dd=dd,h=h,m=m,s=s,           &
                     calkindflag=ESMF_CALKIND_GREGORIAN,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

#ifdef DEBUG
   write(debug,*) 'Leaving TimeStringISOFrac2ESMFtime()'
   write(debug,*)
#endif
   return

   end subroutine TimeStringISOFrac2ESMFtime
!EOC
#endif
!-----------------------------------------------------------------------

   end module getm_oasis

!-----------------------------------------------------------------------
! Copyright (C) 2017 - Xaver Lange and Knut Klingbeil                  !
!-----------------------------------------------------------------------
