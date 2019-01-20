#include "cppdefs.h"
#define _GETM_NUOPC_
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  getm_esmf - routines to integrate GETM into ESMF
!
! !INTERFACE:
   module getm_esmf
!
! !DESCRIPTION:
!  Example for GriddedComponent can be found in
!  ESMFDIR/src/Superstructure/Component/examples/ESMF_GCompEx.F90.
!
! !USES:
   use esmf
   use NUOPC

   IMPLICIT NONE
   private
!
! !PUBLIC DATA MEMBERS:
   public SetServices ! (NUOPC requires this name), must be public
   !public SetVM       ! optional
   public do_getm_esmf
!
! !PRIVATE DATA MEMBERS:
   type type_getmInternalStateStruct
      sequence
   end type

   type type_getmInternalState
      sequence
      type(type_getmInternalStateStruct),pointer :: wrap
   end type

!  Note (KK): __FILE__ includes full path, thus too long for log
   character(len=*),parameter :: FILENAME="getm_esmf.F90"

   character(len=*),parameter :: name_depth  = &
      "water_depth_at_soil_surface"
   character(len=*),parameter :: name_eps3D  = &
      "dissipation_of_tke_in_water"
   character(len=*),parameter :: name_epsbot = &
      "dissipation_of_tke_at_soil_surface"
   character(len=*),parameter :: name_h3D    = &
      "layer_height_in_water"
   character(len=*),parameter :: name_hbot   = &
      "layer_height_at_soil_surface"
   character(len=*),parameter :: name_NN3D   = &
      "buoyancy_frequency_squared_in_water"
   character(len=*),parameter :: name_nuh3D  = &
      "turbulent_diffusivity_of_heat_in_water"
   character(len=*),parameter :: name_num3D  = &
      "turbulent_diffusivity_of_momentum_in_water"
   character(len=*),parameter :: name_numbot = &
      "turbulent_diffusivity_of_momentum_at_soil_surface"
   character(len=*),parameter :: name_S3D    = &
      "salinity_in_water"
   character(len=*),parameter :: name_slp    = &
      "air_pressure_at_sea_level"
   character(len=*),parameter :: name_SS3D   = &
      "shear_frequency_squared_in_water"
   character(len=*),parameter :: name_swr    = &
      "surface_downwelling_photosynthetic_radiative_flux"
   character(len=*),parameter :: name_T3D    = &
      "temperature_in_water"
   character(len=*),parameter :: name_taubmax= &
      "maximum_bottom_stress"
   character(len=*),parameter :: name_Tbot   = &
      "temperature_at_soil_surface"
   character(len=*),parameter :: name_tke3D  = &
      "turbulent_kinetic_energy_in_water"
   character(len=*),parameter :: name_tkebot = &
      "turbulent_kinetic_energy_at_soil_surface"
   character(len=*),parameter :: name_U2D    = &
      "depth_averaged_x_velocity_in_water"
   character(len=*),parameter :: name_U3D    = &
      "x_velocity_in_water"
   character(len=*),parameter :: name_Ubot   = &
      "x_velocity_at_soil_surface"
   character(len=*),parameter :: name_V2D    = &
      "depth_averaged_y_velocity_in_water"
   character(len=*),parameter :: name_V3D    = &
      "y_velocity_in_water"
   character(len=*),parameter :: name_Vbot   = &
      "y_velocity_at_soil_surface"
   character(len=*),parameter :: name_waveDir= &
      "wave_direction"
   character(len=*),parameter :: name_waveH  = &
      "wave_height"
   character(len=*),parameter :: name_waveK  = &
      "wave_number"
   character(len=*),parameter :: name_waveT  = &
      "wave_period"
   character(len=*),parameter :: name_windU  = &
      "wind_x_velocity_at_10m"
   character(len=*),parameter :: name_windV  = &
      "wind_y_velocity_at_10m"
   character(len=*),parameter :: name_airT2  = &
      "air_temperature_at_2m"
   character(len=*),parameter :: name_humr   = &
      "relative_humidity"
   character(len=*),parameter :: name_hums   = &
      "specific_humidity"
   character(len=*),parameter :: name_dev2   = &
      "dew_point_temperature"
   character(len=*),parameter :: name_tcc    = &
      "total_cloud_cover"
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
!EOP
!-----------------------------------------------------------------------

    contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: do_getm_esmf
!
! !INTERFACE:
   subroutine do_getm_esmf()
!
! !DESCRIPTION:
!
! !USES:
   use getm_timers, only: tic,toc,TIM_ESMF
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES
   type(ESMF_Clock)        :: mainClock,getmClock
   type(ESMF_GridComp)     :: getmComp
   type(ESMF_State)        :: importState,exportState
   type(ESMF_TimeInterval) :: runDuration
   integer                 :: phase,phase0,phaseCount
   logical                 :: phaseZeroFlag
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_getm_esmf() # ',Ncall
#endif

   call tic(TIM_ESMF)

   call ESMF_Initialize(defaultLogFileName='getmESMF_LogFile',         &
                        defaultCalKind=ESMF_CALKIND_GREGORIAN)
   call ESMF_LogSet(flush=.true.,                                      &
                    logmsgList=(/ESMF_LOGMSG_WARNING,ESMF_LOGMSG_ERROR/))

   !mainClock = ESMF_ClockCreate(...)

!  Note (KK): contextflag=ESMF_CONTEXT_OWN_VM (default) is required for
!             an individual petList!
!             By default ESMF_GridComp[SetServices|Initialize|Run|Finalize]()
!             are skipped for PETs not in the petList provided
!             to ESMF_GridCompCreate().
   getmComp = ESMF_GridCompCreate(contextflag=ESMF_CONTEXT_PARENT_VM,name="getm")

   importState = ESMF_StateCreate(stateintent=ESMF_STATEINTENT_IMPORT, &
                                  name="getmImportState")
   exportState = ESMF_StateCreate(stateintent=ESMF_STATEINTENT_EXPORT, &
                                  name="getmExportState")

   call ESMF_GridCompSetServices(getmComp,SetServices)

!  Note (KK): check whether a NUOPC-Zero phase is available, inquire
!             InitializePhaseMap (see e.g. NUOPC_Comp.F90::NUOPC_SearchPhaseMap())
   call ESMF_GridCompGetEPPhaseCount(getmComp,ESMF_METHOD_INITIALIZE,  &
                                     phaseCount,phaseZeroFlag)
   phase0=1
   if (phaseZeroFlag) phase0=0

   call toc(TIM_ESMF)

   do phase=phase0,phaseCount
!     Note (KK): The parent component (i.e. this routine) can provide a
!                clock with [start|stop]Time to the getmComp, either by
!                creating the getmComp with this clock or by providing
!                a clock to the InitService.
      call ESMF_GridCompInitialize(getmComp,importState=importState,   &
                                   exportState=exportState,phase=phase)
   end do

   call tic(TIM_ESMF)

!  Note (KK): The InitService of the GETM component created a clock.
!             For PETs not part of GETM's petList, this check is mandatory
!             and will return False!
   !call ESMF_GridCompGet(getmComp,clockIsPresent=clockIsPresent)
   !if (clockIsPresent) then
!     Note (KK): Usually the already created mainClock has only to be
!                adapted (don't forget to update currTime!). However,
!                here the mainClock has not even been created yet.
      call ESMF_GridCompGet(getmComp,clock=getmClock)
      !call ESMF_ClockGet(getmClock,startTime=startTime, &
      !                   stopTime=stopTime,runDuration=runDuration)
      !call ESMF_ClockSet(mainClock,startTime=startTime, &
      !                   stopTime=stopTime,timeStep=runDuration, &
      !                   currTime=startTime)
      mainClock = ESMF_ClockCreate(getmClock)
      call ESMF_ClockGet(mainClock,runDuration=runDuration)
      call ESMF_ClockSet(mainClock,name="mainClock",timeStep=runDuration)
   !end if ! if(clockIsPresent)

   call toc(TIM_ESMF)

   call ESMF_GridCompRun(getmComp,importState=importState,             &
                         exportState=exportState,clock=mainClock)
   call ESMF_GridCompFinalize(getmComp,importState=importState,        &
                              exportState=exportState,clock=mainClock)

   call tic(TIM_ESMF)

   call ESMF_GridCompDestroy(getmComp)

!  Note (KK): Usually mainClock is created by all PETs of the actual VM.
!             However, here it was only created by the PETs in GETM's
!             petList. In case that not all PETs are in GETM's petList,
!             a check of ESMF_GridCompIsPetLocal(getmComp) would be
!             mandatory!
   call ESMF_ClockDestroy(mainClock)

   call ESMF_Finalize()

   call toc(TIM_ESMF)

#ifdef DEBUG
   write(debug,*) 'Leaving do_getm_esmf()'
   write(debug,*)
#endif
   return

   end subroutine do_getm_esmf
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: SetServices - register GriddedComponent GETM
!
! !INTERFACE:
   subroutine SetServices(getmComp,rc)
!
! !DESCRIPTION:
!  Register user-code subroutines.
!  Set the entry points for standard ESMF Component methods
!  ESMF_GridComp[Initialize|Run|Finalize](), called by toplevel component.
!  The toplevel component requires this sub for its mandatory call to
!  ESMF_GridCompSetServices(getmComp,userRoutine=SetServices).
!  For interface see ESMFDIR/src/Superstructure/Component/src/ESMF_GridComp.F90.
!  For the NUOPC Layer this routine must be named SetServices().
!  The toplevel component can inquire rc via optional keyword argument
!  userRc to ESMF_GridCompSetServices().
!
! !USES:
   use getm_timers, only: tic,toc,TIM_ESMF
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp) :: getmComp
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
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'SetServices() # ',Ncall
#endif

   call tic(TIM_ESMF)

   call ESMF_GridCompGet(getmComp,name=name,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)

   call ESMF_LogWrite(trim(name)//"::SetServices...",ESMF_LOGMSG_TRACE)

!  Optional keyword argument "phase" (default: 1) promotes multi-phase user-code.

   call ESMF_GridCompSetEntryPoint(getmComp,ESMF_METHOD_INITIALIZE,    &
                                   userRoutine=InitializeP0,phase=0,   &
                                   rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridCompSetEntryPoint(getmComp,ESMF_METHOD_INITIALIZE,    &
                                   userRoutine=InitializeP1,phase=1,   &
                                   rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridCompSetEntryPoint(getmComp,ESMF_METHOD_INITIALIZE,    &
                                   userRoutine=InitializeP2,phase=2,   &
                                   rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridCompSetEntryPoint(getmComp,ESMF_METHOD_RUN,           &
                                   userRoutine=RunP1,                  &
                                   rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridCompSetEntryPoint(getmComp,ESMF_METHOD_FINALIZE,      &
                                   userRoutine=FinalizeP1,             &
                                   rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!  Optional registration of additional routines for checkpoint and
!  restart functions (ESMF_METHOD_[WRITE|READ]RESTART).
!  ...

   call toc(TIM_ESMF)

#ifdef DEBUG
   write(debug,*) 'Leaving SetServices()'
   write(debug,*)
#endif
   return

   end subroutine SetServices
!EOC
!-----------------------------------------------------------------------
#if 0
!BOP
!
! !ROUTINE: SetVM -
!
! !INTERFACE:
   subroutine SetVM(getmComp,rc)
!
! !DESCRIPTION:
!  The toplevel component requires this sub for the optional call to
!  ESMF_GridCompSetVM(getmComp,userRoutine=SetVM), which has to
!  be done *before* the call to ESMF_GridCompSetServices()!
!  For interface see ESMFDIR/src/Superstructure/Component/src/ESMF_GridComp.F90.
!  The toplevel component can inquire rc via optional keyword argument
!  userRc to ESMF_GridCompSetVM().
!
! !USES:
   use getm_timers, only: tic,toc,TIM_ESMF
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp) :: getmComp
!
! !OUTPUT PARAMETERS:
   integer,intent(out) :: rc
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES
   type(ESMF_VM) :: gVM
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'SetVM() # ',Ncall
#endif

   call ESMF_VMGetGlobal(gVM,rc=rc)
!  Calls to ESMF_VMGet, ESMF_GridCompSetVMMaxPEs,
!  ESMF_GridCompSetVM[Min|Max]Threads to modify VM of component

#ifdef DEBUG
   write(debug,*) 'Leaving SetVM()'
   write(debug,*)
#endif
   return

   end subroutine SetVM
!EOC
#endif
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: InitializeP0 -
!
! !INTERFACE:
   subroutine InitializeP0(getmComp,importState,exportState,clock,rc)
!
! !DESCRIPTION:
!  NUOPC initialises all components in various phases. By default NUOPC
!  assumes IPDv00 and thus requires userRoutines for init phases 1 and 2
!  (other phases are added/executed by NUOPC). If the userCode provides
!  different phases, this needs to be communicated to NUOPC via
!  InitializePhaseMap in an init phase 0.
!
!                                               | IPDv00 | IPDv01 | IPDv02
!  =============================================|========|========|========
!   add InitializePhaseMap attribute            | 0 (*?) | 0 (*?) | 0 (*?)
!  =============================================|========|========|========
!   advertise import and export fields          | p1 (!) | p1     | p1
!  =============================================|========|========|========
!   realize import and export fields (allocate) | p2 (!) | p3     | p3
!  =============================================|========|========|========
!   check field' Connected status               | p3     | p4     | p4
!   set internal clock (*)                      |        |        |
!  =============================================|========|========|========
!   initialize fields (*)                       | p4     | p5     | p5
!  =============================================|========|========|========
!
!  (!) has to be done by the user, (*) optional by the user
!
!  Note: states and clock are uninitialized if the toplevel component
!        did not provide corresponding arguments to
!        ESMF_GridCompInitialize(getmComp).
!  The toplevel component can inquire rc via optional keyword argument
!  userRc to ESMF_GridCompInitialize().
!
! !USES:
   use getm_timers, only: tic,toc,TIM_ESMF
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
   logical                                   :: abort
   character(len=ESMF_MAXSTR)                :: name="getm"
   character(len=NUOPC_PhaseMapStringLength) :: InitializePhaseMap(2)
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'InitializeP0() # ',Ncall
#endif

   call tic(TIM_ESMF)

   call ESMF_GridCompGet(getmComp,name=name,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)

   call ESMF_LogWrite(trim(name)//"::InitializeP0...",ESMF_LOGMSG_TRACE)

   InitializePhaseMap(1) = "IPDv00p1=1"
   InitializePhaseMap(2) = "IPDv00p2=2"

#ifdef _GETM_NUOPC_
   call NUOPC_CompAttributeAdd(getmComp,(/"InitializePhaseMap"/),rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call NUOPC_CompAttributeSet(getmComp,"InitializePhaseMap",          &
                               InitializePhaseMap,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
#else
!  Note (KK): NUOPC attributes are purpose="Instance"
   call ESMF_AttributeAdd(getmComp,convention="NUOPC",purpose="General", &
                          attrList=(/"InitializePhaseMap"/),rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(getmComp,name="InitializePhaseMap",          &
                          valueList=InitializePhaseMap,                &
                          convention="NUOPC",purpose="General",rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
#endif

   call toc(TIM_ESMF)

#ifdef DEBUG
   write(debug,*) 'Leaving InitializeP0()'
   write(debug,*)
#endif
   return

   end subroutine InitializeP0
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: InitializeP1 -
!
! !INTERFACE:
   subroutine InitializeP1(getmComp,importState,exportState,clock,rc)
!
! !DESCRIPTION:
!  Note: states and clock are uninitialized if the toplevel component
!        did not provide corresponding arguments to
!        ESMF_GridCompInitialize(getmComp).
!        Status can be checked with [State|Clock]IsCreated().
!  The toplevel component can inquire rc via optional keyword argument
!  userRc to ESMF_GridCompInitialize().
!
! !USES:
#ifdef GETM_PARALLEL
   use mpi
   use halo_mpi   , only: comm_getm
#endif
   use time       , only: init_time,start,timestep
   use initialise , only: init_initialise,do_initialise,dryrun
   use integration, only: MinN,MaxN
   use getm_timers, only: tic,toc,TIM_ESMF
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
   type(ESMF_Time)            :: getmRefTime,getmStartTime,getmStopTime
   type(ESMF_TimeInterval)    :: getmTimeStep
   type(ESMF_VM)              :: vm
   type(type_getmInternalState)               :: gis
   type(type_getmInternalStateStruct),pointer :: isd=>NULL()
   integer                    :: localrc,comm,length,getmRunTimeStepCount
   logical                    :: abort,clockIsPresent
   character(len=ESMF_MAXSTR) :: name="getm"
   character(len=8)           :: datestr
   character(len=10)          :: timestr
   character(len=19)          :: TimeStrISOFrac,start_external,stop_external
#ifdef GETM_PARALLEL
    character(len=MPI_MAX_ERROR_STRING) :: mpierrmsg
#endif
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'InitializeP1() # ',Ncall
#endif

   call tic(TIM_ESMF)

   call ESMF_GridCompGet(getmComp,name=name,                           &
                         clockIsPresent=clockIsPresent,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_LogWrite(trim(name)//"::InitializeP1...",ESMF_LOGMSG_TRACE)

!  Set Internal State
   allocate(isd)
   gis%wrap => isd
   call ESMF_GridCompSetInternalState(getmComp,gis,rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!  Optional Creation and Initialization of child components
!  (Create, SetServices, Initialize)
!  ...

!  Check whether toplevel component called ESMF_GridCompCreate() with clock.
   if (clockIsPresent) then
      call ESMF_GridCompGet(getmComp,clock=getmClock,rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

      clockIsPresent = ESMF_ClockIsCreated(getmClock,rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
   end if
   if (.not. clockIsPresent) then
      clockIsPresent = ESMF_ClockIsCreated(clock,rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

      if (clockIsPresent) then
         getmClock = ESMF_ClockCreate(clock,rc=rc)
         abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
         if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

         call ESMF_ClockSet(getmClock,name=trim(name)//"Clock",rc=rc)
         abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
         if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

         call ESMF_GridCompSet(getmComp,clock=getmClock,rc=rc)
         abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
         if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
   end if


!  This is where the model specific setup code goes
!  (allocation, open files, initial conditions).

#ifdef GETM_PARALLEL
   call ESMF_GridCompGet(getmComp,vm=vm,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_VMGet(vm,mpiCommunicator=comm,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call MPI_COMM_DUP(comm,comm_getm,rc)
   if (rc .ne. MPI_SUCCESS) then
!     need depends on specified mpi error handler (i.e. not MPI_ERRORS_ARE_FATAL)
      call MPI_ERROR_STRING(rc,mpierrmsg,length,localrc)
      call ESMF_LogWrite(trim(mpierrmsg),                              &
                         ESMF_LOGMSG_ERROR,line=__LINE__,file=FILENAME)
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
   end if
#endif

   call date_and_time(datestr,timestr)


   if (clockIsPresent) then

!     Use startTime and stopTime from already initialised getmClock.
      call ESMF_ClockGet(getmClock,startTime=getmStartTime,            &
                         stopTime=getmStopTime,rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

      call ESMF_TimeGet(getmStartTime,timeStringISOFrac=start_external,rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

      call ESMF_TimeGet(getmStopTime,timeStringISOFrac=stop_external,rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

      call toc(TIM_ESMF)

      call init_initialise(datestr,timestr)
      call init_time(MinN,MaxN,start_external=start_external, &
                     stop_external=stop_external)
      call do_initialise()

      call tic(TIM_ESMF)

!     use internal GETM time step
      call ESMF_TimeIntervalSet(getmTimeStep,                          &
                                s_r8=real(timestep,kind=ESMF_KIND_R8), &
                                rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

      call ESMF_ClockSet(getmClock,timeStep=getmTimeStep,rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   else

      call toc(TIM_ESMF)

      ! set up clock based on internal GETM specifications
      call init_initialise(datestr,timestr)
      call init_time(MinN,MaxN)
      call do_initialise()

      call tic(TIM_ESMF)

!     reference time
      TimeStrISOFrac=start(1:10)//"T"//start(12:19)
      !call ESMF_TimeSet(refTime,timeStringISOFrac=TimeStrISOFrac)
      call TimeStringISOFrac2ESMFtime(TimeStrISOFrac,getmRefTime)

!     time step
      call ESMF_TimeIntervalSet(getmTimeStep,                          &
                                s_r8=real(timestep,kind=ESMF_KIND_R8), &
                                rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

      getmStartTime = getmRefTime + (MinN-1)*getmTimeStep
      getmRunTimeStepCount = MaxN - MinN + 1

!     create clock based on internal GETM specifications
      getmClock = ESMF_ClockCreate(getmTimeStep,getmStartTime,            &
                                   runTimeStepCount=getmRunTimeStepCount, &
                                   refTime=getmRefTime,                   &
                                   name='getmClock',rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

      call ESMF_GridCompSet(getmComp,clock=getmClock,rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   end if

   call getmComp_init_grid(getmComp)
   call NUOPC_FieldDictionarySetAutoAdd(.true.)
   call init_importStateP1(getmComp,importState)
   call init_exportStateP1(getmComp,exportState)

   call toc(TIM_ESMF)

#ifdef DEBUG
   write(debug,*) 'Leaving InitializeP1()'
   write(debug,*)
#endif
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
   use initialise, only: finalise_initialise,dryrun
   use getm_timers, only: tic,toc,TIM_ESMF
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
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'InitializeP2() # ',Ncall
#endif

   call tic(TIM_ESMF)

   call ESMF_GridCompGet(getmComp,name=name,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)

   call ESMF_LogWrite(trim(name)//"::InitializeP2...",ESMF_LOGMSG_TRACE)

   call init_exportStateP2(getmComp,exportState)

!  If the initial Export state needs to be filled, do it here.
   call getmComp_update_grid(getmComp)
   call update_exportState(getmComp,exportState)

   call toc(TIM_ESMF)

   call finalise_initialise()

   if (.not.dryrun) then
      STDERR LINE
      LEVEL1 'integrating....'
      STDERR LINE
   end if

   rc = ESMF_SUCCESS

#ifdef DEBUG
   write(debug,*) 'Leaving InitializeP2()'
   write(debug,*)
#endif
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
   use initialise , only: runtype,dryrun
   use integration, only: time_step,MinN,MaxN
   use getm_timers, only: tic,toc,TIM_ESMF
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
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'RunP1() # ',Ncall
#endif

   call tic(TIM_ESMF)

   call ESMF_GridCompGet(getmComp,name=name,clock=getmClock,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_LogWrite(trim(name)//"::RunP1...",ESMF_LOGMSG_TRACE)

!  call ESMF_StateGet(), etc to get fields, bundles, arrays from import state
   call read_importState(getmComp,importState)

   call ESMF_ClockGet(getmClock,timeStep=getmTimeStep,currtime=getmTime, &
                      advanceCount=advanceCount,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!  use clock to do determine time of calling routine
   call ESMF_ClockGetNextTime(clock,nextTime,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   do while (getmTime + 0.5d0*getmTimeStep < nextTime )

      n = MinN + int(advanceCount,kind=kind(MinN))
      if (n .gt. MaxN ) then
         call ESMF_LogWrite('already reached MaxN',                    &
                            ESMF_LOGMSG_WARNING,line=__LINE__,file=FILENAME)
         !call ESMF_Finalize(endflag=ESMF_END_ABORT)
         rc = ESMF_RC_NOT_SET
         return
      end if

      call toc(TIM_ESMF)

!     optional Run of child components
!     ...

!     This is where the model specific computation goes.
      if (.not.dryrun) call time_step(runtype,n)

      call tic(TIM_ESMF)

      call ESMF_ClockAdvance(getmClock,rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

      call ESMF_ClockGet(getmClock,currtime=getmTime,                  &
                         advanceCount=advanceCount,rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   end do

!  Fill export state here using ESMF_StateAdd(), etc
   call getmComp_update_grid(getmComp)
   call update_exportState(getmComp,exportState)

   call toc(TIM_ESMF)

#ifdef DEBUG
   write(debug,*) 'Leaving RunP1()'
   write(debug,*)
#endif
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
   use getm_timers, only: write_getm_timers
   use getm_timers, only: tic,toc,TIM_ESMF
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
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'FinalizeP1() # ',Ncall
#endif

   call tic(TIM_ESMF)

   call ESMF_GridCompGet(getmComp,name=name,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)

   call ESMF_LogWrite(trim(name)//"::FinalizeP1...",ESMF_LOGMSG_TRACE)

!  optional Finalize and Destruction of child components
!  ...

!  Add whatever code here needed (deallocation,close files,flush results)
   call toc(TIM_ESMF)
   call clean_up(dryrun,runtype,MaxN)
#ifndef NO_TIMERS
   STDERR LINE
   call write_getm_timers
#endif
   call tic(TIM_ESMF)

!  Get Internal State
   call ESMF_GridCompGetInternalState(getmComp,gis,rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (.not. abort) then
!     Deallocate private data block
      isd => gis%wrap
      !deallocate(isd%waveH)
      deallocate(isd)
   end if

   call ESMF_GridCompGet(getmComp,clockIsPresent=clockIsPresent,       &
                                  gridIsPresent=GridIsPresent,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (.not. abort) then
      if (GridIsPresent) then
         call ESMF_GridCompGet(getmComp,grid=getmGrid,rc=rc)
         abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
         if (.not. abort) then
            call ESMF_GridGet(getmGrid,distgrid=getmDistGrid3D,rc=rc)
            abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
            if (.not. abort) then
               !call ESMF_DistGridDestroy(getmDistGrid2D,rc=rc)
               !abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
               call ESMF_DistGridDestroy(getmDistGrid3D,rc=rc)
               abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
            end if
            !call ESMF_GridGetCoord(getmGrid,coordDim=...,staggerloc=...,array=array)
            !call ESMF_ArrayDestroy(array)
            call ESMF_GridDestroy(getmGrid,rc=rc)
            abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
         end if
      end if
      if (clockIsPresent) then
         call ESMF_GridCompGet(getmComp,clock=getmClock,rc=rc)
         abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
         if (.not. abort) then
            call ESMF_ClockDestroy(getmClock,rc=rc)
            abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
         end if
      end if
   end if

   call toc(TIM_ESMF)

   rc = ESMF_SUCCESS

#ifdef DEBUG
   write(debug,*) 'Leaving FinalizeP1()'
   write(debug,*)
#endif
   return

   end subroutine FinalizeP1
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: getmComp_init_grid - Creates Grid
!
! !INTERFACE:
   subroutine getmComp_init_grid(getmComp)
!
! !DESCRIPTION:
!
! !USES:
   use initialise, only: runtype
   use domain    , only: ioff,joff,imin,jmin,imax,jmax,kmax
   use domain    , only: grid_type,have_lonlat
   use domain    , only: az,ax,areac_getm=>areac
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
! !REVISION HISTORY:
!  Original Author(s): Knut Klingbeil
!
! !LOCAL VARIABLES
   type(ESMF_Array),target  :: xcArray2D,ycArray2D,xxArray2D,yxArray2D
   type(ESMF_Array),target  :: xcArray3D,ycArray3D,xxArray3D,yxArray3D
   type(ESMF_Array),target  :: loncArray2D,latcArray2D,lonxArray2D,latxArray2D
   type(ESMF_Array),target  :: loncArray3D,latcArray3D,lonxArray3D,latxArray3D
   type(ESMF_Array),pointer :: coordArrayXC2D=>NULL()
   type(ESMF_Array),pointer :: coordArrayXC3D=>NULL()
   type(ESMF_Array),pointer :: coordArrayYC2D=>NULL()
   type(ESMF_Array),pointer :: coordArrayYC3D=>NULL()
   type(ESMF_Array),pointer :: coordArrayXX2D=>NULL()
   type(ESMF_Array),pointer :: coordArrayXX3D=>NULL()
   type(ESMF_Array),pointer :: coordArrayYX2D=>NULL()
   type(ESMF_Array),pointer :: coordArrayYX3D=>NULL()
   type(ESMF_Array)         :: array
   type(ESMF_CoordSys_Flag) :: coordSys
   type(ESMF_DistGrid)      :: getmDistGrid2D,getmDistGrid3D,distgrid
   type(ESMF_Grid)          :: getmGrid3D,getmGrid2D
   type(ESMF_StaggerLoc)    :: staggerloc
   type(ESMF_VM)            :: vm
!  Note (KK): ESMF_ARRAY's are deep classes, that persist after return.
!             (even without save attribute), same for allocated pointers.
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
   integer(ESMF_KIND_I4),dimension(:),allocatable,target :: alledges
   integer(ESMF_KIND_I4),dimension(4),target             :: myedges
   integer(ESMF_KIND_I4),pointer                         :: p2di(:,:)
   REALTYPE                 :: getmreal
   integer                  :: rc,pet,i0,j0,ilen,jlen,klen,k
   integer                  :: petCount
   integer,dimension(3)     :: coordDimCount
   integer,dimension(3,3)   :: coordDimMap
   integer,dimension(:,:,:),allocatable                  :: deBlockList
   logical                    :: abort,noKindMatch
   character(len=ESMF_MAXSTR) :: name="getm"
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'getmComp_init_grid() # ',Ncall
#endif

   call ESMF_GridCompGet(getmComp,name=name,vm=vm,petCount=petCount,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   myedges = (/ ioff , joff , imax , jmax /)
   allocate(alledges(4*petCount))

!  syncflag=ESMF_SYNC_BLOCKING (default)
   call ESMF_VMAllGather(vm,myedges,alledges,4,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   if (runtype .eq. 1) then
      klen = 1
   else
      klen = kmax
   end if

   allocate(deBlockList(3,2,petCount))
   do pet = 0,petCount-1
      i0   = 1 + alledges(1+4*pet)
      j0   = 1 + alledges(2+4*pet)
      ilen =     alledges(3+4*pet)
      jlen =     alledges(4+4*pet)
      deBlockList(:,1,1+pet) = (/ i0        , j0        , 1        /)
      deBlockList(:,2,1+pet) = (/ i0+ilen-1 , j0+jlen-1 , 1+klen-1 /)
   end do

!  indexflag=ESMF_INDEX_DELOCAL (default) starting at 1
!  (for ESMF_INDEX_USER [grid|stagger]MemLBound can be set)
#if 1
!  Single-tile DistGrid (1 subdomain = 1 DE)
!  internal call to ESMF_DistGridCreateDB()
   getmDistGrid2D = ESMF_DistGridCreate(minval(deBlockList(1:2,1,:),2), &
                                        maxval(deBlockList(1:2,2,:),2), &
                                        int(deBlockList(1:2,:,:)),rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
   getmDistGrid3D = ESMF_DistGridCreate(minval(deBlockList(:,1,:),2), &
                                        maxval(deBlockList(:,2,:),2), &
                                        deBlockList,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
#else
!  Multi-tile DistGrid (1 subdomain = 1 tile = 1 DE) by specification of
!  [min|max]IndexPTile.
!  Note (KK): int() intrinsic routines are needed, because ESMF does not
!             accept subarrays as arguments
!  internal call to ESMF_DistGridCreateRDT()
   getmDistGrid2D = ESMF_DistGridCreate(int(deBlockList(1:2,1,:)), &
                                        int(deBlockList(1:2,2,:)),rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
   getmDistGrid3D = ESMF_DistGridCreate(int(deBlockList(:,1,:)), &
                                        int(deBlockList(:,2,:)),rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
#endif

   noKindMatch = ( kind(getmreal) .ne. ESMF_KIND_R8 )

!  KK-TODO: get 2d farrayPtr from created persistent deep 3D Array
   allocate(maskC3D(E2DFIELD,0:klen))
   allocate(maskX3D(E2DFIELD,0:klen))
   allocate(areaW3D(E2DFIELD,0:klen))

   p2di => maskC3D(:,:,klen) ; maskC(imin-HALO:,jmin-HALO:) => p2di
   p2di => maskX3D(:,:,klen) ; maskX(imin-HALO:,jmin-HALO:) => p2di
   p2dr => areaW3D(:,:,klen) ; areaC(imin-HALO:,jmin-HALO:) => p2dr

   maskC = az
   maskX = ax

#ifdef SLICE_MODEL
   maskC(:,jmax/2+1) = 0
   maskX(:,jmax/2+1) = 0
#endif

   areaC = areac_getm

   do k=0,klen
      maskC3D(:,:,k) = maskC
      maskX3D(:,:,k) = maskX
      areaW3D(:,:,k) = areaC
   end do
   maskC3D(:,:,0) = 0

   allocate(zx(E2DXFIELD,0:klen))
   if (noKindMatch .or. runtype.eq.1) then
      allocate(zw(E2DFIELD,0:klen))
      allocate(zc(E2DFIELD,0:klen))
   else
#ifndef NO_3D
      zw => zwn
      zc => zcn
#endif
   end if

   select case (grid_type)

!        coordDimMap: for each grid dimension i map array dimension j
!                     to grid dimension coordDimMap(i,j)
!                     i=1,dimCount ; j=1,coordDimCount(i)

      case(1)

         coordSys = ESMF_COORDSYS_CART
         coordDimCount = (/ 1 , 1 , 3 /)     ! rectilinear horizontal coordinates
         coordDimMap = reshape( (/1,2,1,0,0,2,0,0,3/) , (/3,3/) )

      case(2)

         coordSys = ESMF_COORDSYS_SPH_DEG    ! (default)
         coordDimCount = (/ 1 , 1 , 3 /)     ! rectilinear horizontal coordinates
         coordDimMap = reshape( (/1,2,1,0,0,2,0,0,3/) , (/3,3/) )

      case(3)

         coordSys = ESMF_COORDSYS_CART
         coordDimCount = (/ 2 , 2 , 3 /)
         coordDimMap = reshape( (/1,1,1,2,2,2,0,0,3/) , (/3,3/) ) ! (default)

      case(4)

         coordSys = ESMF_COORDSYS_SPH_DEG                         ! (default)
         coordDimCount = (/ 2 , 2 , 3 /)
         coordDimMap = reshape( (/1,1,1,2,2,2,0,0,3/) , (/3,3/) ) ! (default)

   end select

!  Note (KK): gridAlign specifies which corner point in a grid cell
!             shares the center indices [ default=(/-1,...,-1/) ].
!             gridEdge[L|U]Width only affect DE's at the edge of tiles
!             (thus it matters whether a single- or multi-tile DistGrid
!              was created). If gridEdgeWidth's are not set, they are set
!             automatically based on gridAlign.
!             For single staggered grid items, default specification for
!             gridAlign and gridEdgesWidth's can be overwritten by
!             staggerAlign and staggerEdgeWidth's.
!  internal call to ESMF_GridCreateFrmDistGrid()
   getmGrid2D = ESMF_GridCreate(getmDistGrid2D,name=trim(name)//"Grid2D", &
#if 0
! bug in ESMF
                                gridAlign=(/1,1/),                        &
#else
                                gridEdgeLWidth=(/1,1/),                   &
#endif
                                coordSys=coordSys,                        &
                                coordDimCount=int(coordDimCount(1:2)),    &
                                coordDimMap=int(coordDimMap(1:2,1:2)),    &
                                rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   getmGrid3D = ESMF_GridCreate(getmDistGrid3D,name=trim(name)//"Grid3D", &
#if 0
! bug in ESMF
                                gridAlign=(/1,1,1/),                      &
#else
                                gridEdgeLWidth=(/1,1,1/),                 &
#endif

                                coordSys=coordSys,                        &
                                coordDimCount=coordDimCount,              &
                                coordDimMap=coordDimMap,                  &
                                rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   select case (grid_type)

      case(1)

         if (noKindMatch) then
            allocate(xc1D(   _IRANGE_HALO_)) ; xc1D = xcord
            allocate(yc1D(   _JRANGE_HALO_)) ; yc1D = ycord
            allocate(xx1D(-1+_IRANGE_HALO_)) ; xx1D = xxcord
            allocate(yx1D(-1+_JRANGE_HALO_)) ; yx1D = yxcord
         else
            xc1D => xcord
            yc1D => ycord
            xx1D => xxcord
            yx1D => yxcord
         end if

!        put pointers into persistent deep Arrays
         call createCoordArrays1D(xc1d,yc1d,xx1d,yx1d,                     &
                                  getmGrid2D,getmGrid3D,                   &
                                  xcArray2D,ycArray2D,xxArray2D,yxArray2D, &
                                  xcArray3D,ycArray3D,xxArray3D,yxArray3D)
         call setCoordUnits("m","m",                                   &
                            xcArray2D,ycArray2D,xxArray2D,yxArray2D,   &
                            xcArray3D,ycArray3D,xxArray3D,yxArray3D)

         coordArrayXC2D => xcArray2D ; coordArrayYC2D => ycArray2D
         coordArrayXX2D => xxArray2D ; coordArrayYX2D => yxArray2D
         coordArrayXC3D => xcArray3D ; coordArrayYC3D => ycArray3D
         coordArrayXX3D => xxArray3D ; coordArrayYX3D => yxArray3D

      case(2)

         if (noKindMatch) then
            allocate(lonc1D(   _IRANGE_HALO_)) ; lonc1D = xcord
            allocate(latc1D(   _JRANGE_HALO_)) ; latc1D = ycord
            allocate(lonx1D(-1+_IRANGE_HALO_)) ; lonx1D = xxcord
            allocate(latx1D(-1+_JRANGE_HALO_)) ; latx1D = yxcord
         else
            lonc1D => xcord
            latc1D => ycord
            lonx1D => xxcord
            latx1D => yxcord
         end if

!        put pointers into persistent deep Arrays
         call createCoordArrays1D(lonc1d,latc1d,lonx1d,latx1d,                     &
                                  getmGrid2D,getmGrid3D,                           &
                                  loncArray2D,latcArray2D,lonxArray2D,latxArray2D, &
                                  loncArray3D,latcArray3D,lonxArray3D,latxArray3D)
         call setCoordUnits("degrees_east","degrees_north",                  &
                            loncArray2D,latcArray2D,lonxArray2D,latxArray2D, &
                            loncArray3D,latcArray3D,lonxArray3D,latxArray3D)

         coordArrayXC2D => loncArray2D ; coordArrayYC2D => latcArray2D
         coordArrayXX2D => lonxArray2D ; coordArrayYX2D => latxArray2D
         coordArrayXC3D => loncArray3D ; coordArrayYC3D => latcArray3D
         coordArrayXX3D => lonxArray3D ; coordArrayYX3D => latxArray3D

      case(3)

         if (noKindMatch) then
            allocate(xx2D(E2DXFIELD)) ; xx2D = xx
            allocate(yx2D(E2DXFIELD)) ; yx2D = yx
            allocate(xc2D(E2DFIELD )) ; xc2D = xc
            allocate(yc2D(E2DFIELD )) ; yc2D = yc
            if (have_lonlat) then
               allocate(lonx2D(E2DXFIELD)) ; lonx2D = lonx
               allocate(latx2D(E2DXFIELD)) ; latx2D = latx
               allocate(lonc2D(E2DFIELD )) ; lonc2D = lonc
               allocate(latc2D(E2DFIELD )) ; latc2D = latc
            end if
         else
            xx2D => xx
            yx2D => yx
            xc2D => xc
            yc2D => yc
            if (have_lonlat) then
               lonx2D => lonx
               latx2D => latx
               lonc2D => lonc
               latc2D => latc
            end if
         end if

!        put pointers into persistent deep Arrays
         call createCoordArrays2D(xc2d,yc2d,xx2d,yx2d,                     &
                                  getmGrid2D,getmGrid3D,                   &
                                  xcArray2D,ycArray2D,xxArray2D,yxArray2D, &
                                  xcArray3D,ycArray3D,xxArray3D,yxArray3D)
         call setCoordUnits("m","m",                                   &
                            xcArray2D,ycArray2D,xxArray2D,yxArray2D,   &
                            xcArray3D,ycArray3D,xxArray3D,yxArray3D)

         if (have_lonlat) then
            call createCoordArrays2D(lonc2d,latc2d,lonx2d,latx2d,                     &
                                     getmGrid2D,getmGrid3D,                           &
                                     loncArray2D,latcArray2D,lonxArray2D,latxArray2D, &
                                     loncArray3D,latcArray3D,lonxArray3D,latxArray3D)
            call setCoordUnits("degrees_east","degrees_north",                  &
                               loncArray2D,latcArray2D,lonxArray2D,latxArray2D, &
                               loncArray3D,latcArray3D,lonxArray3D,latxArray3D)
         end if

         coordArrayXC2D => xcArray2D ; coordArrayYC2D => ycArray2D
         coordArrayXX2D => xxArray2D ; coordArrayYX2D => yxArray2D
         coordArrayXC3D => xcArray3D ; coordArrayYC3D => ycArray3D
         coordArrayXX3D => xxArray3D ; coordArrayYX3D => yxArray3D

         if (have_lonlat) then
            coordArrayXC2D => loncArray2D ; coordArrayYC2D => latcArray2D
            coordArrayXX2D => lonxArray2D ; coordArrayYX2D => latxArray2D
            coordArrayXC3D => loncArray3D ; coordArrayYC3D => latcArray3D
            coordArrayXX3D => lonxArray3D ; coordArrayYX3D => latxArray3D
         end if

      case(4)

         if (noKindMatch) then
            allocate(lonx2D(E2DXFIELD)) ; lonx2D = lonx
            allocate(latx2D(E2DXFIELD)) ; latx2D = latx
            allocate(lonc2D(E2DFIELD )) ; lonc2D = lonc
            allocate(latc2D(E2DFIELD )) ; latc2D = latc
         else
            lonx2D => lonx
            latx2D => latx
            lonc2D => lonc
            latc2D => latc
         end if

!        put pointers into persistent deep Arrays
         call createCoordArrays2D(lonc2d,latc2d,lonx2d,latx2d,                     &
                                  getmGrid2D,getmGrid3D,                           &
                                  loncArray2D,latcArray2D,lonxArray2D,latxArray2D, &
                                  loncArray3D,latcArray3D,lonxArray3D,latxArray3D)
         call setCoordUnits("degrees_east","degrees_north",                  &
                            loncArray2D,latcArray2D,lonxArray2D,latxArray2D, &
                            loncArray3D,latcArray3D,lonxArray3D,latxArray3D)

         coordArrayXC2D => loncArray2D ; coordArrayYC2D => latcArray2D
         coordArrayXX2D => lonxArray2D ; coordArrayYX2D => latxArray2D
         coordArrayXC3D => loncArray3D ; coordArrayYC3D => latcArray3D
         coordArrayXX3D => lonxArray3D ; coordArrayYX3D => latxArray3D

   end select

!  2D grid

   staggerloc = ESMF_STAGGERLOC_CENTER ! (default)
   call ESMF_GridGet(getmGrid2D, staggerloc, distgrid=distgrid, rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridSetCoord(getmGrid2D,1,array=coordArrayXC2D,           &
                          staggerloc=staggerloc,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridSetCoord(getmGrid2D,2,array=coordArrayYC2D,           &
                          staggerloc=staggerloc,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!  2D center mask
   array = ESMF_ArrayCreate(distgrid,maskC,                            &
                            indexflag=ESMF_INDEX_DELOCAL,              &
                            datacopyflag=ESMF_DATACOPY_VALUE,          &
                            rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!  KK-TODO: add attribute with un-/mask value
   call ESMF_GridSetItem(getmGrid2D,ESMF_GRIDITEM_MASK,array=array,    &
                         staggerloc=staggerloc,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!  2D center area
   array = ESMF_ArrayCreate(distgrid,areaC,                            &
                            indexflag=ESMF_INDEX_DELOCAL,              &
                            datacopyflag = ESMF_DATACOPY_VALUE,        &
                            rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(array,'units','m**2',rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridSetItem(getmGrid2D,ESMF_GRIDITEM_AREA,array=array,    &
                         staggerloc=StaggerLoc, rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   staggerloc = ESMF_STAGGERLOC_CORNER
   call ESMF_GridGet(getmGrid2D, staggerloc, distgrid=distgrid, rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridSetCoord(getmGrid2D,1,array=coordArrayXX2D,           &
                          staggerloc=staggerloc,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridSetCoord(getmGrid2D,2,array=coordArrayYX2D,           &
                          staggerloc=staggerloc,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!  2D corner mask
   array = ESMF_ArrayCreate(distgrid,maskX,                            &
                            indexflag=ESMF_INDEX_DELOCAL,              &
                            datacopyflag=ESMF_DATACOPY_VALUE,          &
                            rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!  KK-TODO: add attribute with un-/mask value
   call ESMF_GridSetItem(getmGrid2D,ESMF_GRIDITEM_MASK,array=array,    &
                         staggerloc=staggerloc,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

#if 0
call ESMF_OutputScripGridFile('GETM_Grid_Update.nc', getmGrid2D, rc)
abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
#endif

!  3D grid

   staggerloc = ESMF_STAGGERLOC_CENTER_VCENTER ! (=CENTER, default for 3D)
   call ESMF_GridGet(getmGrid3D, staggerloc, distgrid=distgrid, rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridSetCoord(getmGrid3D,1,array=coordArrayXC3D,           &
                          staggerloc=staggerloc,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridSetCoord(getmGrid3D,2,array=coordArrayYC3D,           &
                          staggerloc=staggerloc,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!  3D center position
   array = ESMF_ArrayCreate(distgrid,zc,                               &
                            indexflag=ESMF_INDEX_DELOCAL,              &
                            totalLWidth=(/HALO,HALO,1/),               &
                            totalUWidth=(/HALO,HALO,0/),               &
                            rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(array,'units','m',rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridSetCoord(getmGrid3D,3,array=array,                    &
                          staggerloc=staggerloc,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!  3D center mask
!  KK-TODO: get rid of k=0 layer (implemented by MOSSCO commit 5e04bd4)
   array = ESMF_ArrayCreate(distgrid,maskC3D,                          &
                            indexflag=ESMF_INDEX_DELOCAL,              &
                            datacopyflag=ESMF_DATACOPY_VALUE,          &
                            totalLWidth=(/HALO,HALO,1/),               &
                            totalUWidth=(/HALO,HALO,0/),               &
                            rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!  KK-TODO: add attribute with un-/mask value
   call ESMF_GridSetItem(getmGrid3D,ESMF_GRIDITEM_MASK,array=array,    &
                         staggerloc=staggerloc,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!  3D center area
   !array = ESMF_ArrayCreate(getmDistGrid3D,areaC,indexflag=ESMF_INDEX_DELOCAL,rc=rc)
!  KK-TODO: k=0 layer needed because of areaW3D (needed by MOSSCO)
   array = ESMF_ArrayCreate(distgrid,areaW3D,                          &
                            indexflag=ESMF_INDEX_DELOCAL,              &
                            datacopyflag = ESMF_DATACOPY_VALUE,        &
                            totalLWidth=(/HALO,HALO,1/),               &
                            totalUWidth=(/HALO,HALO,0/),               &
                            rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(array,'units','m**2',rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridSetItem(getmGrid3D,ESMF_GRIDITEM_AREA,array=array,    &
                         staggerloc=staggerloc, rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   staggerloc = ESMF_STAGGERLOC_CORNER_VFACE
   call ESMF_GridGet(getmGrid3D, staggerloc, distgrid=distgrid, rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridSetCoord(getmGrid3D,1,array=coordArrayXX3D,           &
                          staggerloc=staggerloc,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridSetCoord(getmGrid3D,2,array=coordArrayYX3D,           &
                          staggerloc=staggerloc,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!  3D corner position
!  KK-TODO: this call defines zx[1:kmax+1],
!           [0:kmax] not possible because exclusive region is spanned
!           from 1:kmax+1
   array = ESMF_ArrayCreate(distgrid,zx,                               &
                            indexflag=ESMF_INDEX_DELOCAL,              &
                            datacopyflag = ESMF_DATACOPY_VALUE,        &
                            rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(array,'units','m',rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridSetCoord(getmGrid3D,3,array=array,                    &
                          staggerloc=staggerloc,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!  3D corner mask
   !array = ESMF_ArrayCreate(getmDistGrid3D,maskX,indexflag=ESMF_INDEX_DELOCAL,rc=rc)
!  Note (KK): no HALO+1 in totalLWidth because based on E2DFIELD
!  KK-TODO: this call defines maskX3D[1:kmax+1],
!           [0:kmax] not possible because exclusive region is spanned
!           from 1:kmax+1
   array = ESMF_ArrayCreate(distgrid,maskX3D,                          &
                            indexflag=ESMF_INDEX_DELOCAL,              &
                            datacopyflag = ESMF_DATACOPY_VALUE,        &
                            rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!  KK-TODO: add attribute with un-/mask value
   call ESMF_GridSetItem(getmGrid3D,ESMF_GRIDITEM_MASK,array=array,    &
                         staggerloc=staggerloc,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   staggerloc = ESMF_STAGGERLOC_CENTER_VFACE
   call ESMF_GridGet(getmGrid3D, staggerloc, distgrid=distgrid, rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!  3D interface position
!  KK-TODO: this call defines zw[1:kmax+1],
!           [0:kmax] not possible because exclusive region is spanned
!           from 1:kmax+1
   array = ESMF_ArrayCreate(distgrid,zw,                               &
                            indexflag=ESMF_INDEX_DELOCAL,              &
                            rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(array,'units','m',rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridSetCoord(getmGrid3D,3,array=array,                    &
                          staggerloc=staggerloc,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!  3D interface area
!  KK-TODO: this call defines areaw3D[1:kmax+1],
!           [0:kmax] not possible because exclusive region is spanned
!           from 1:kmax+1
   !array = ESMF_ArrayCreate(getmDistGrid3D,areaC,indexflag=ESMF_INDEX_DELOCAL,rc=rc)
   array = ESMF_ArrayCreate(distgrid,areaW3D,                          &
                            indexflag=ESMF_INDEX_DELOCAL,              &
                            datacopyflag = ESMF_DATACOPY_VALUE,        &
                            rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeSet(array,'units','m**2',rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridSetItem(getmGrid3D,ESMF_GRIDITEM_AREA,array=array,    &
                         staggerloc=staggerloc, rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)


   call ESMF_GridCompSet(getmComp,gridList=(/getmGrid3D,getmGrid2D/),rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call getmComp_update_grid(getmComp)

   deallocate(alledges)
   deallocate(deBlockList)

#ifdef DEBUG
   write(debug,*) 'Leaving getmComp_init_grid()'
   write(debug,*)
#endif
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
   integer                                     :: rc,i,j,k,klen,k0
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
                          staggerloc=ESMF_STAGGERLOC_CENTER_VFACE,     &
                          rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridGetCoord(getmGrid3D,3,farrayPtr=zc,                   &
                          staggerloc=ESMF_STAGGERLOC_CENTER_VCENTER,   &
                          rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridGetCoord(getmGrid3D,3,farrayPtr=zx,                   &
                          staggerloc=ESMF_STAGGERLOC_CORNER_VFACE,     &
                          rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   if (runtype .eq. 1) then
      klen = 1
   else
      klen = kmax
   end if

   k0 = lbound(zx,3)

   noKindMatch = ( kind(getmreal) .ne. ESMF_KIND_R8 )

   if (noKindMatch .or. runtype.eq.1) then
      if (runtype .eq. 1) then
         zw(:,:,k0+0)    = -H
         zw(:,:,k0+klen) = z
         zc(:,:,klen) = 0.5d0 * ( zw(:,:,k0+klen-1) + zw(:,:,k0+klen) )
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
               zwu(i,j) = 0.5d0 * ( zw(i,j,k0+k) + zw(i+1,j,k0+k) )
            else
               if (az(i,j) .gt. 0) then
                  zwu(i,j) = zw(i,j,k0+k)
               else if (az(i+1,j) .gt. 0) then
                  zwu(i,j) = zw(i+1,j,k0+k)
               end if
            end if
         end do
      end do
      do j=jmin-HALO,jmax+HALO-1
         do i=imin-HALO,imax+HALO-1
            if (av(i,j).gt.0 .or. av(i+1,j).gt.0) then
               zx(i,j,k0+k) = 0.5d0 * ( zwu(i,j) + zwu(i,j+1) )
            else
               if (az(i,j).gt.0 .or. az(i+1,j).gt.0) then
                  zx(i,j,k0+k) = zwu(i,j)
               else if (az(i,j+1).gt.0 .or. az(i+1,j+1).gt.0) then
                  zx(i,j,k0+k) = zwu(i,j+1)
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
   subroutine init_importStateP1(getmComp,importState)
!
! !DESCRIPTION:
!
! !USES:
   use initialise     ,only: runtype
   use domain         ,only: grid_type
   use meteo          ,only: met_method,calc_met,METEO_FROMEXT
   use meteo          ,only: airp,u10,v10,t2,hum,tcc
   use waves          ,only: waveforcing_method,WAVES_FROMEXT
   use waves          ,only: waves_ramp
   use variables_waves,only: waveH,waveK,waveT
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
   type(ESMF_Grid) :: getmGrid3D,getmGrid2D
   type(ESMF_Grid), allocatable :: gridList(:)
   integer         :: rc
   logical         :: abort,frc
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_importStateP1() # ',Ncall
#endif

   call ESMF_GridCompGet(getmComp, gridList=gridList, rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   getmGrid3D = gridList(1)
   getmGrid2D = gridList(2)

   if (met_method .eq. METEO_FROMEXT) then
      !call StateAddField(importState,trim(name_swr    ),getmGrid2D,units="W m-2")
      call StateAddField(importState,trim(name_slp    ),getmGrid2D,    &
                         farray2D=airp,units="Pa")
      if (calc_met) then
!        force allocation of new memory if grid rotation needs to be removed
         frc = (grid_type .ne. 2)
         call StateAddField(importState,trim(name_windU  ),getmGrid2D, &
                            farray2D=u10,units="m s-1",frc=frc)
         call StateAddField(importState,trim(name_windV  ),getmGrid2D, &
                            farray2D=v10,units="m s-1",frc=frc)
         if (runtype .gt. 2) then
         call StateAddField(importState,trim(name_airT2  ),getmGrid2D, &
                            farray2D=t2,units="K")
         call StateAddField(importState,trim(name_hums   ),getmGrid2D, &
                            farray2D=hum,units="kg/kg")
         call StateAddField(importState,trim(name_humr   ),getmGrid2D, &
                            farray2D=hum,units="%")
         call StateAddField(importState,trim(name_dev2   ),getmGrid2D, &
                            farray2D=hum,units="K")
         call StateAddField(importState,trim(name_tcc    ),getmGrid2D, &
                            farray2D=hum,units="")
         end if
      end if ! calc_met
   end if ! meteo

   if (waveforcing_method .eq. WAVES_FROMEXT) then
      call StateAddField(importState,trim(name_waveDir),getmGrid2D,    &
                         units="rad",frc=.true.)
!     force allocation of new memory if waveH is ramped
      frc = (waves_ramp .gt. 1)
      call StateAddField(importState,trim(name_waveH  ),getmGrid2D,    &
                         farray2D=waveH,units="m",frc=frc)
      call StateAddField(importState,trim(name_waveK  ),getmGrid2D,    &
                         farray2D=waveK,units="m-1")
      call StateAddField(importState,trim(name_waveT  ),getmGrid2D,    &
                         farray2D=waveT,units="s")
   end if

   !call ESMF_StatePrint(importState)

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
   use initialise     ,only: runtype
   use domain         ,only: grid_type,convc,cosconv,sinconv
   use meteo          ,only: met_method,calc_met,METEO_FROMEXT,new_meteo
   use meteo          ,only: hum_method,RELATIVE_HUM,WET_BULB,DEW_POINT,SPECIFIC_HUM
   use meteo          ,only: airp,u10,v10,t2,hum,tcc
   use waves          ,only: waveforcing_method,WAVES_FROMEXT,new_waves
   use waves          ,only: waves_ramp
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

         call StateReadCompleteField(importState,trim(name_slp    ),   &
                                     farray2d=airp)

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
         if (runtype .gt. 2) then
            call StateReadCompleteField(importState,trim(name_airT2  ),&
                                        farray2d=t2)
            p2dr => NULL()
            call StateReadCompleteField(importState,trim(name_hums   ),&
                                        p2dr=p2dr,frc=.true.,ign=.true.)
            if (associated(p2dr)) then
            hum_method = SPECIFIC_HUM
            else
            call StateReadCompleteField(importState,trim(name_humr   ),&
                                        p2dr=p2dr,frc=.true.,ign=.true.)
            if (associated(p2dr)) then
            hum_method = RELATIVE_HUM
            else
            call StateReadCompleteField(importState,trim(name_dev2   ),&
                                        p2dr=p2dr,frc=.true.,ign=.true.)
            if (associated(p2dr)) then
            hum_method = DEW_POINT
            else
            hum_method = -1
            call ESMF_LogWrite('hum_method=-1',                        &
                               ESMF_LOGMSG_WARNING,line=__LINE__,file=FILENAME)

            end if
            end if
            end if
            call StateReadCompleteField(importState,trim(name_tcc    ),&
                                        farray2d=tcc)
         end if
      end if ! calc_met

   end if ! meteo


   if (waveforcing_method .eq. WAVES_FROMEXT) then

         new_waves = .true. ! KK-TODO: should be set by coupler

         call StateReadCompleteField(importState,trim(name_waveDir),   &
                                     p2dr=p2dr,frc=.true.)
         coswavedir = cos( p2dr + convc*deg2rad )
         sinwavedir = sin( p2dr + convc*deg2rad )

         frc = (waves_ramp .gt. 1)
         call StateReadCompleteField(importState,trim(name_waveH  ),   &
                                     farray2d=waveH,frc=frc)
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
   subroutine init_exportStateP1(getmComp,exportState)
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
! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_GridComp) :: getmComp
   type(ESMF_State)    :: exportState
!
! !REVISION HISTORY:
!  Original Author(s): Knut Klingbeil
!
! !LOCAL VARIABLES
   type(ESMF_Grid) :: getmGrid3D,getmGrid2D
   type(ESMF_Grid), allocatable :: gridList(:)
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

   call ESMF_GridCompGet(getmComp, gridList=gridList, rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   getmGrid3D = gridList(1)
   getmGrid2D = gridList(2)

   !call StateAddField(exportState,trim(name_depth  ),getmGrid2D,units="m")
   !call StateAddField(exportState,trim(name_h3D    ),getmGrid3D,units="m")
   !call StateAddField(exportState,trim(name_hbot   ),getmGrid2D,units="m")
   !call StateAddField(exportState,trim(name_U2D    ),getmGrid2D,units="m s-1")
   !call StateAddField(exportState,trim(name_U3D    ),getmGrid3D,units="m s-1")
   !call StateAddField(exportState,trim(name_Ubot   ),getmGrid2D,units="m s-1")
   !call StateAddField(exportState,trim(name_V2D    ),getmGrid2D,units="m s-1")
   !call StateAddField(exportState,trim(name_V3D    ),getmGrid3D,units="m s-1")
   !call StateAddField(exportState,trim(name_Vbot   ),getmGrid2D,units="m s-1")
   if (runtype .ge. 2) then
#ifndef NO_3D
   !call StateAddField(exportState,trim(name_eps3D  ),getmGrid3D,units="m2 s-3",staggerloc=ESMF_STAGGERLOC_CENTER_VFACE)
   !call StateAddField(exportState,trim(name_epsbot ),getmGrid2D,units="m2 s-3")
   !call StateAddField(exportState,trim(name_num3D  ),getmGrid3D,units="m2 s-1",staggerloc=ESMF_STAGGERLOC_CENTER_VFACE)
   !call StateAddField(exportState,trim(name_numbot ),getmGrid2D,units="m2 s-1")
   !call StateAddField(exportState,trim(name_SS3D   ),getmGrid3D,units="s-2",staggerloc=ESMF_STAGGERLOC_CENTER_VFACE)
   !call StateAddField(exportState,trim(name_taubmax),getmGrid2D,units="Pa")
   !call StateAddField(exportState,trim(name_tke3D  ),getmGrid3D,units="m2 s-2",staggerloc=ESMF_STAGGERLOC_CENTER_VFACE)
   !call StateAddField(exportState,trim(name_tkebot ),getmGrid2D,units="m2 s-2")
   if (runtype .ge. 3) then
#ifndef NO_BAROCLINIC
   !call StateAddField(exportState,trim(name_NN3D   ),getmGrid3D,units="s-2",staggerloc=ESMF_STAGGERLOC_CENTER_VFACE)
   !call StateAddField(exportState,trim(name_nuh3D  ),getmGrid3D,units="m",staggerloc=ESMF_STAGGERLOC_CENTER_VFACE)
   if (calc_salt) then
   !call StateAddField(exportState,trim(name_S3D    ),getmGrid3D,units="")
   end if
   if (calc_temp) then
   !call StateAddField(exportState,trim(name_T3D    ),getmGrid3D,units="degC")
   !call StateAddField(exportState,trim(name_Tbot   ),getmGrid2D,units="degC")
   end if
#endif
   end if
#endif
   end if
   if (met_method.eq.METEO_CONST .or. met_method.eq.METEO_FROMFILE) then
   !call StateAddField(exportState,trim(name_swr    ),getmGrid2D,units="W m-2")
   !call StateAddField(exportState,trim(name_windU  ),getmGrid2D,units="m s-1")
   !call StateAddField(exportState,trim(name_windV  ),getmGrid2D,units="m s-1")
   end if
   if (     waveforcing_method.eq.WAVES_FROMWIND                       &
       .or. waveforcing_method.eq.WAVES_FROMFILE) then
   !call StateAddField(exportState,trim(name_waveDir),getmGrid2D,units="rad")
   !call StateAddField(exportState,trim(name_waveH  ),getmGrid2D,units="m")
   !call StateAddField(exportState,trim(name_waveK  ),getmGrid2D,units="m-1")
   !call StateAddField(exportState,trim(name_waveT  ),getmGrid2D,units="s")
   end if

   !call ESMF_StatePrint(exportState)

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
! !ROUTINE: init_exportStateP2 -
!
! !INTERFACE:
   subroutine init_exportStateP2(getmComp,exportState)
!
! !DESCRIPTION:
!
! !USES:
   use initialise     ,only: runtype
   use domain         ,only: grid_type
   use variables_2d   ,only: D,velx,vely
#ifndef NO_3D
   use m3d            ,only: calc_salt,calc_temp
   use variables_3d   ,only: Dn,hn,velx3d,vely3d,velx2dadv,vely2dadv
   use variables_3d   ,only: nuh
#ifndef NO_BAROCLINIC
   use variables_3d   ,only: S,T
#endif
#endif
   use meteo          ,only: met_method,METEO_CONST,METEO_FROMFILE
   use meteo          ,only: swr,u10,v10
   use waves          ,only: waveforcing_method
   use waves          ,only: NO_WAVES,WAVES_FROMWIND,WAVES_FROMFILE
   use waves          ,only: waves_method,WAVES_NOSTOKES
   use variables_waves,only: waveH,waveK,waveT
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
   logical                    :: frc
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_exportStateP2() # ',Ncall
#endif

!  force allocation of new memory if Stokes drift needs to be removed
   frc = (waveforcing_method.ne.NO_WAVES .and. waves_method.ne.WAVES_NOSTOKES)

   if (runtype .eq. 1) then

!     2D-only and general 3D variables

      call StateCompleteConnectedField(exportState,trim(name_depth  ), &
                                       farray2d=D)
      p3dr => NULL()
      call StateCompleteConnectedField(exportState,trim(name_h3D    ), &
                                       p3dr=p3dr)
      if (associated(p3dr)) then
         p2dr => p3dr(:,:,1)
      else
         p2dr => NULL()
      end if
      call StateCompleteConnectedField(exportState,trim(name_hbot   ), &
                                       farray2d=D,p2dr=p2dr)
      call StateCompleteConnectedField(exportState,trim(name_U2D    ), &
                                       farray2D=velx,frc=frc)
      p3dr => NULL()
      call StateCompleteConnectedField(exportState,trim(name_U3D    ), &
                                       p3dr=p3dr)
      if (associated(p3dr)) then
         p2dr => p3dr(:,:,1)
      else
         p2dr => NULL()
      end if
      call StateCompleteConnectedField(exportState,trim(name_Ubot   ), &
                                       farray2d=velx,p2dr=p2dr,frc=frc)
      call StateCompleteConnectedField(exportState,trim(name_V2D    ), &
                                       farray2D=vely,frc=frc)
      p3dr => NULL()
      call StateCompleteConnectedField(exportState,trim(name_V3D    ), &
                                       p3dr=p3dr)
      if (associated(p3dr)) then
         p2dr => p3dr(:,:,1)
      else
         p2dr => NULL()
      end if
      call StateCompleteConnectedField(exportState,trim(name_Vbot   ), &
                                       farray2d=vely,p2dr=p2dr,frc=frc)

   else ! runtype.gt.1

#ifndef NO_3D

      call StateCompleteConnectedField(exportState,trim(name_depth  ), &
                                       farray2d=Dn)
      p3dr => NULL()
      call StateCompleteConnectedField(exportState,trim(name_h3D    ), &
                                       farray3D=hn,p3dr=p3dr)
      if (associated(p3dr)) then
         p2dr => p3dr(:,:,1)
      else
         p2dr => NULL()
      end if
      call StateCompleteConnectedField(exportState,trim(name_hbot   ), &
                                       farray2D=hn(:,:,1),p2dr=p2dr)
      call StateCompleteConnectedField(exportState,trim(name_U2D    ), &
                                       farray2D=velx2dadv,frc=frc)
      p3dr => NULL()
      call StateCompleteConnectedField(exportState,trim(name_U3D    ), &
                                       farray3D=velx3d,p3dr=p3dr,frc=frc)
      if (associated(p3dr)) then
         p2dr => p3dr(:,:,1)
      else
         p2dr => NULL()
      end if
      call StateCompleteConnectedField(exportState,trim(name_Ubot   ), &
                                       farray2D=velx3d(:,:,1),p2dr=p2dr,frc=frc)
      call StateCompleteConnectedField(exportState,trim(name_V2D    ), &
                                       farray2D=vely2dadv,frc=frc)
      p3dr => NULL()
      call StateCompleteConnectedField(exportState,trim(name_V3D    ), &
                                       farray3D=vely3d,p3dr=p3dr,frc=frc)
      if (associated(p3dr)) then
         p2dr => p3dr(:,:,1)
      else
         p2dr => NULL()
      end if
      call StateCompleteConnectedField(exportState,trim(name_Vbot   ), &
                                       farray2D=vely3d(:,:,1),p2dr=p2dr,frc=frc)

!     additional 3D-only variables

      p3dr => NULL()
      call StateCompleteConnectedField(exportState,trim(name_eps3D  ), &
                                       p3dr=p3dr)
      if (associated(p3dr)) then
         p2dr => p3dr(:,:,1)
      else
         p2dr => NULL()
      end if
      call StateCompleteConnectedField(exportState,trim(name_epsbot ), &
                                       p2dr=p2dr)
      p3dr => NULL()
      call StateCompleteConnectedField(exportState,trim(name_num3D  ), &
                                       p3dr=p3dr)
      if (associated(p3dr)) then
         p2dr => p3dr(:,:,1)
      else
         p2dr => NULL()
      end if
      call StateCompleteConnectedField(exportState,trim(name_numbot ), &
                                       p2dr=p2dr)
      call StateCompleteConnectedField(exportState,trim(name_SS3D   ))
      call StateCompleteConnectedField(exportState,trim(name_taubmax))
      p3dr => NULL()
      call StateCompleteConnectedField(exportState,trim(name_tke3D  ), &
                                       p3dr=p3dr)
      if (associated(p3dr)) then
         p2dr => p3dr(:,:,1)
      else
         p2dr => NULL()
      end if
      call StateCompleteConnectedField(exportState,trim(name_tkebot ), &
                                       p2dr=p2dr)

      if (runtype .ge. 3) then

#ifndef NO_BAROCLINIC

      call StateCompleteConnectedField(exportState,trim(name_NN3D   ))
      call StateCompleteConnectedField(exportState,trim(name_nuh3D  ), &
                                       farray3D=nuh)

      if (calc_salt) then
      call StateCompleteConnectedField(exportState,trim(name_S3D    ), &
                                       farray3D=S)
      end if ! calc_salt

      if (calc_temp) then
      p3dr => NULL()
      call StateCompleteConnectedField(exportState,trim(name_T3D    ), &
                                       farray3D=T,p3dr=p3dr)
      if (associated(p3dr)) then
         p2dr => p3dr(:,:,1)
      else
         p2dr => NULL()
      end if
      call StateCompleteConnectedField(exportState,trim(name_Tbot   ), &
                                       p2dr=p2dr)
      end if ! calc_temp

#endif

      end if ! runtype.ge.3

#endif

   end if ! runtype


   if (met_method.eq.METEO_CONST .or. met_method.eq.METEO_FROMFILE) then
      call StateCompleteConnectedField(exportState,trim(name_swr    ), &
                                       farray2D=swr)
!     force allocation of new memory if grid rotation needs to be removed
      frc = (grid_type .ne. 2)
      call StateCompleteConnectedField(exportState,trim(name_windU  ), &
                                       farray2D=u10,frc=frc)
      call StateCompleteConnectedField(exportState,trim(name_windV  ), &
                                       farray2D=v10,frc=frc)
   end if ! meteo

   if (     waveforcing_method.eq.WAVES_FROMWIND                       &
       .or. waveforcing_method.eq.WAVES_FROMFILE) then
      call StateCompleteConnectedField(exportState,trim(name_waveDir))
      call StateCompleteConnectedField(exportState,trim(name_waveH  ), &
                                       farray2D=waveH)
      call StateCompleteConnectedField(exportState,trim(name_waveK  ), &
                                       farray2D=waveK)
      call StateCompleteConnectedField(exportState,trim(name_waveT  ), &
                                       farray2D=waveT)
   end if ! waves

   !call ESMF_StatePrint(exportState)

#ifdef DEBUG
   write(debug,*) 'Leaving init_exportStateP2()'
   write(debug,*)
#endif
   return

   end subroutine init_exportStateP2
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
         call StateSetCompleteField(exportState,trim(name_swr    ),    &
                                    farray2D=swr)
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

         call StateSetCompleteField(exportState,trim(name_waveH  ),    &
                                    farray2D=waveH)
         call StateSetCompleteField(exportState,trim(name_waveK  ),    &
                                    farray2D=waveK)
         call StateSetCompleteField(exportState,trim(name_waveT  ),    &
                                    farray2D=waveT)
   end if ! waves

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
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: createCoordArrays1D - 
!
! !INTERFACE:
   subroutine createCoordArrays1D(xc1d,yc1d,xx1d,yx1d,                     &
                                  getmGrid2D,getmGrid3D,                   &
                                  xcArray2D,ycArray2D,xxArray2D,yxArray2D, &
                                  xcArray3D,ycArray3D,xxArray3D,yxArray3D)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(ESMF_KIND_R8),dimension(:),pointer,intent(in) :: xc1d,yc1d,xx1d,yx1d
   type(ESMF_Grid),intent(in)                         :: getmGrid2D,getmGrid3D
!
! !OUTPUT PARAMETERS:
   type(ESMF_Array),intent(out) :: xcArray2D,ycArray2D,xxArray2D,yxArray2D
   type(ESMF_Array),intent(out) :: xcArray3D,ycArray3D,xxArray3D,yxArray3D
!
! !REVISION HISTORY:
!  Original Author(s): Knut Klingbeil
!
! !LOCAL VARIABLES
   type(ESMF_DistGrid)      :: distgrid
   type(ESMF_StaggerLoc)    :: staggerloc
   logical                  :: abort
   integer                  :: rc
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'createCoordArrays1D() # ',Ncall
#endif

!  Note (KK): datacopyflag=ESMF_DATACOPY_REFERENCE (default)
!             DATACOPY_VALUE not needed, because source memory persists
!             when pointer goes out of scope
!
!  1D xcord is replicated automatically along 2nd DistGrid dimension
!  1D ycord is replicated along 1st DistGrid dimension as specified
!  by distgridToArrayMap.
!  total[L|U]Width are automatically determined from shape of [x|y]cord

!  internal call to ESMF_ArrayCreateAssmdShape<rank><type><kind>()
!  (because of required indexflag)
!  Note (KK): These ArrayCreate()'s only work for 1DE per PET!!!
!             Automatically determined coordDimMap for rectilinear
!             coordinates is incorrect!
!             automatically determined total[L|U]Width are not
!             consistent with specified gridAlign

!  2D grid

   staggerloc = ESMF_STAGGERLOC_CENTER
   call ESMF_GridGet(getmGrid2D, staggerloc, distgrid=distgrid, rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   xcArray2D = ESMF_ArrayCreate(distgrid,xc1D,                         &
                                indexflag=ESMF_INDEX_DELOCAL,          &
                                rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   ycArray2D = ESMF_ArrayCreate(distgrid,yc1D,                         &
                                indexflag=ESMF_INDEX_DELOCAL,          &
                                distgridToArrayMap=(/0,1/),            &
                                rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   staggerloc = ESMF_STAGGERLOC_CORNER
   call ESMF_GridGet(getmGrid2D, staggerloc, distgrid=distgrid, rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   xxArray2D = ESMF_ArrayCreate(distgrid,xx1D,                         &
                                indexflag=ESMF_INDEX_DELOCAL,          &
!                                totalLWidth=(/HALO+1/),                &
                                totalUWidth=(/HALO  /),                &
                                rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   yxArray2D = ESMF_ArrayCreate(distgrid,yx1D,                         &
                                indexflag=ESMF_INDEX_DELOCAL,          &
                                distgridToArrayMap=(/0,1/),            &
!                                totalLWidth=(/HALO+1/),                &
                                totalUWidth=(/HALO  /),                &
                                rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)


!  3D grid

   staggerloc = ESMF_STAGGERLOC_CENTER
   call ESMF_GridGet(getmGrid3D, staggerloc, distgrid=distgrid, rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   xcArray3D = ESMF_ArrayCreate(distgrid,xc1D,                         &
                                indexflag=ESMF_INDEX_DELOCAL,          &
                                rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   ycArray3D = ESMF_ArrayCreate(distgrid,yc1D,                         &
                                indexflag=ESMF_INDEX_DELOCAL,          &
                                distgridToArrayMap=(/0,1,0/),          &
                                rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   staggerloc = ESMF_STAGGERLOC_CORNER
   call ESMF_GridGet(getmGrid3D, staggerloc, distgrid=distgrid, rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   xxArray3D = ESMF_ArrayCreate(distgrid,xx1D,                         &
                                indexflag=ESMF_INDEX_DELOCAL,          &
!                                totalLWidth=(/HALO+1/),                &
                                totalUWidth=(/HALO  /),                &
                                rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   yxArray3D = ESMF_ArrayCreate(distgrid,yx1D,                         &
                                indexflag=ESMF_INDEX_DELOCAL,          &
                                distgridToArrayMap=(/0,1,0/),          &
!                                totalLWidth=(/HALO+1/),                &
                                totalUWidth=(/HALO  /),                &
                                rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

#ifdef DEBUG
   write(debug,*) 'Leaving createCoordArrays1D()'
   write(debug,*)
#endif
   return

   end subroutine createCoordArrays1D
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: createCoordArrays2D - 
!
! !INTERFACE:
   subroutine createCoordArrays2D(xc2d,yc2d,xx2d,yx2d,                     &
                                  getmGrid2D,getmGrid3D,                   &
                                  xcArray2D,ycArray2D,xxArray2D,yxArray2D, &
                                  xcArray3D,ycArray3D,xxArray3D,yxArray3D)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(ESMF_KIND_R8),dimension(:,:),pointer,intent(in) :: xc2d,yc2d,xx2d,yx2d
   type(ESMF_Grid),intent(in)                           :: getmGrid2D,getmGrid3D
!
! !OUTPUT PARAMETERS:
   type(ESMF_Array),intent(out) :: xcArray2D,ycArray2D,xxArray2D,yxArray2D
   type(ESMF_Array),intent(out) :: xcArray3D,ycArray3D,xxArray3D,yxArray3D
!
! !REVISION HISTORY:
!  Original Author(s): Knut Klingbeil
!
! !LOCAL VARIABLES
   type(ESMF_DistGrid)      :: distgrid
   type(ESMF_StaggerLoc)    :: staggerloc
   logical                  :: abort
   integer                  :: rc
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'createCoordArrays2D() # ',Ncall
#endif

!  Note (KK): datacopyflag=ESMF_DATACOPY_REFERENCE (default)
!             DATACOPY_VALUE not needed, because source memory persists
!             when pointer goes out of scope
!
!  2D farrays are replicated automatically along 3rd DistGrid dimension
!  internal call to ESMF_ArrayCreateAssmdShape<rank><type><kind>()
!  (because of required indexflag)
!  Note (KK): These ArrayCreate()'s only work for 1DE per PET!!!
!             automatically determined total[L|U]Width are not
!             consistent with specified gridAlign

!  2D grid

   staggerloc = ESMF_STAGGERLOC_CENTER
   call ESMF_GridGet(getmGrid2D, staggerloc, distgrid=distgrid, rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   xcArray2D = ESMF_ArrayCreate(distgrid,xc2D,                         &
                                indexflag=ESMF_INDEX_DELOCAL,          &
                                rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   ycArray2D = ESMF_ArrayCreate(distgrid,yc2D,                         &
                                indexflag=ESMF_INDEX_DELOCAL,          &
                                rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   staggerloc = ESMF_STAGGERLOC_CORNER
   call ESMF_GridGet(getmGrid2D, staggerloc, distgrid=distgrid, rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   xxArray2D = ESMF_ArrayCreate(distgrid,xx2D,                         &
                                indexflag=ESMF_INDEX_DELOCAL,          &
!                                totalLWidth=(/HALO+1,HALO+1/),         &
                                totalUWidth=(/HALO  ,HALO  /),         &
                                rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   yxArray2D = ESMF_ArrayCreate(distgrid,yx2D,                         &
                                indexflag=ESMF_INDEX_DELOCAL,          &
!                                totalLWidth=(/HALO+1,HALO+1/),         &
                                totalUWidth=(/HALO  ,HALO  /),         &
                                rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)


!  3D grid

   staggerloc = ESMF_STAGGERLOC_CENTER
   call ESMF_GridGet(getmGrid3D, staggerloc, distgrid=distgrid, rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   xcArray3D = ESMF_ArrayCreate(distgrid,xc2D,                         &
                                indexflag=ESMF_INDEX_DELOCAL,          &
                                rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   ycArray3D = ESMF_ArrayCreate(distgrid,yc2D,                         &
                                indexflag=ESMF_INDEX_DELOCAL,          &
                                rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   staggerloc = ESMF_STAGGERLOC_CORNER
   call ESMF_GridGet(getmGrid3D, staggerloc, distgrid=distgrid, rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   xxArray3D = ESMF_ArrayCreate(distgrid,xx2D,                         &
                                indexflag=ESMF_INDEX_DELOCAL,          &
!                                totalLWidth=(/HALO+1,HALO+1/),         &
                                totalUWidth=(/HALO  ,HALO  /),         &
                                rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   yxArray3D = ESMF_ArrayCreate(distgrid,yx2D,                         &
                                indexflag=ESMF_INDEX_DELOCAL,          &
!                                totalLWidth=(/HALO+1,HALO+1/),         &
                                totalUWidth=(/HALO  ,HALO  /),         &
                                rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

#ifdef DEBUG
   write(debug,*) 'Leaving createCoordArrays2D()'
   write(debug,*)
#endif
   return

   end subroutine createCoordArrays2D
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: setCoordUnits - 
!
! !INTERFACE:
   subroutine setCoordUnits(xunits,yunits,                             &
                            xcArray2D,ycArray2D,xxArray2D,yxArray2D,   &
                            xcArray3D,ycArray3D,xxArray3D,yxArray3D)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*),intent(in) :: xunits,yunits
!
! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_Array) :: xcArray2D,ycArray2D,xxArray2D,yxArray2D
   type(ESMF_Array) :: xcArray3D,ycArray3D,xxArray3D,yxArray3D
!
! !REVISION HISTORY:
!  Original Author(s): Knut Klingbeil
!
! !LOCAL VARIABLES
   logical                  :: abort
   integer                  :: rc
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'setCoordUnits() # ',Ncall
#endif

   call ESMF_AttributeSet(xcArray2D,'units',trim(xunits),rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
   call ESMF_AttributeSet(xcArray3D,'units',trim(xunits),rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
   call ESMF_AttributeSet(ycArray2D,'units',trim(yunits),rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
   call ESMF_AttributeSet(ycArray3D,'units',trim(yunits),rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
   call ESMF_AttributeSet(xxArray2D,'units',trim(xunits),rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
   call ESMF_AttributeSet(xxArray3D,'units',trim(xunits),rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
   call ESMF_AttributeSet(yxArray2D,'units',trim(yunits),rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
   call ESMF_AttributeSet(yxArray3D,'units',trim(yunits),rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

#ifdef DEBUG
   write(debug,*) 'Leaving setCoordUnits()'
   write(debug,*)
#endif
   return

   end subroutine setCoordUnits
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: StateAddField -
!
! !INTERFACE:
   subroutine StateAddField(state,name,grid,kwe,units,staggerloc,      &
                            farray2D,farray3D,frc)
!
! !DESCRIPTION:
!  Creates empty or complete field in state.
!
! !USES:
   use domain    ,only: imin,jmin,imax,jmax,kmax
   use initialise,only: runtype
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*)     ,intent(in)          :: name
   type(ESMF_Grid)      ,intent(in)          :: grid
   logical              ,intent(in),optional :: kwe !keyword-enforcer
   character(len=*)     ,intent(in),optional :: units
   type(ESMF_StaggerLoc),intent(in),optional :: staggerloc
   REALTYPE,target      ,intent(in),optional :: farray2D(I2DFIELD)
   REALTYPE,target      ,intent(in),optional :: farray3D(I3DFIELD)
   logical              ,intent(in),optional :: frc
!
! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_State)                          :: state
!
! !REVISION HISTORY:
!  Original Author(s): Knut Klingbeil
!
! !LOCAL VARIABLES
   type(ESMF_Field)           :: field
   real(ESMF_KIND_R8),pointer :: p2dr_(:,:),p3dr_(:,:,:)
   REALTYPE                   :: getmreal
   integer                    :: rc,dimCount,klen
   integer,target             :: elb(3),eub(3)
   logical                    :: abort,frc_,isCreated,noKindMatch
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'StateAddField() # ',Ncall
#endif

   if (present(frc)) then
      frc_ = frc
   else
      frc_ = .false.
   end if

   isCreated = ESMF_StateIsCreated(state,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   if (.not. isCreated) then
      call ESMF_LogWrite('skip adding of field '//trim(name)//'',  &
                         ESMF_LOGMSG_WARNING,line=__LINE__,file=FILENAME)
      !call ESMF_Finalize(endflag=ESMF_END_ABORT)
      return
   end if

#ifdef _GETM_NUOPC_
   if (NUOPC_FieldDictionaryHasEntry(trim(name))) then
      call NUOPC_Advertise(state, trim(name), Units=units, rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

      call ESMF_StateGet(state, trim(name), field=field, rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
   else
      call NUOPC_Advertise(state, trim(name), rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

      call ESMF_StateGet(state, trim(name), field=field, rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

      if (present(units)) then
         call ESMF_AttributeSet(field,'units',trim(units),rc=rc)
         abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
         if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
   end if
#else
   field = ESMF_FieldEmptyCreate(name=trim(name),rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   if (present(units)) then
      call ESMF_AttributeSet(field,'units',trim(units),rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
   end if

   call ESMF_AttributeAdd(field,convention="NUOPC",purpose="Instance", &
                          attrList=(/"Connected"/),rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call NUOPC_SetAttribute(field,"Connected","false",rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
#endif

   call ESMF_FieldEmptySet(field,grid,staggerloc=staggerloc,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridGet(grid,dimCount=dimCount,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   noKindMatch = ( kind(getmreal) .ne. ESMF_KIND_R8 )

   if (dimCount .eq. 2) then

      p2dr_ => NULL()
      if (frc_) then
         allocate(p2dr_(I2DFIELD))
         if (present(farray2D)) p2dr_ = farray2D
      else if (present(farray2D) .and. .not.noKindMatch) then
         p2dr_(imin-HALO:,jmin-HALO:) => farray2D
      end if
      if (associated(p2dr_)) then
         call ESMF_GridGetFieldBounds(grid,staggerloc=staggerloc,      &
                                      totalLBound=elb(1:2),totalUBound=eub(1:2),rc=rc)
         abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
         if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!        in contrast to ESMF_ArrayCreate() no automatic determination of total[L|U]Width
#if 1
!        Note (KK): in former times ESMF_FieldCreateGridDataPtr<rank><type><kind>() failed
         call ESMF_FieldEmptyComplete(field,p2dr_,                             &
#else
!        internal call to ESMF_FieldCreateGridData<rank><type><kind>()
!        forced by indexflag argument.
         call ESMF_FieldEmptyComplete(field,p2dr_,ESMF_INDEX_DELOCAL,          &
#endif
                                      totalLWidth=int(elb(1:2)-lbound(p2dr_)), &
                                      totalUWidth=int(ubound(p2dr_)-eub(1:2)), &
                                      rc=rc)
         abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
         if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

   else if (dimCount .eq. 3) then

      p3dr_ => NULL()
      if (runtype .eq. 1) then
         klen = 1
      else
         klen = kmax
      end if
      if (frc_) then
         allocate(p3dr_(I2DFIELD,0:klen))
         if (present(farray3D)) p3dr_ = farray3D
      else if (present(farray3D) .and. .not.noKindMatch) then
         p3dr_(imin-HALO:,jmin-HALO:,0:) => farray3D
      end if
      if (associated(p3dr_)) then
         call ESMF_GridGetFieldBounds(grid,staggerloc=staggerloc,      &
                                      totalLBound=elb,totalUBound=eub,rc=rc)
         abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
         if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!        in contrast to ESMF_ArrayCreate() no automatic determination of total[L|U]Width
#if 1
!        Note (KK): in former times ESMF_FieldCreateGridDataPtr<rank><type><kind>() failed
         call ESMF_FieldEmptyComplete(field,p3dr_,                        &
#else
!        internal call to ESMF_FieldCreateGridData<rank><type><kind>()
!        forced by indexflag argument.
         call ESMF_FieldEmptyComplete(field,p3dr_,ESMF_INDEX_DELOCAL,     &
#endif
                                      totalLWidth=int(elb-lbound(p3dr_)), &
                                      totalUWidth=int(ubound(p3dr_)-eub), &
                                      rc=rc)
         abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
         if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

   end if

#ifndef _GETM_NUOPC_
   call ESMF_StateAdd(state,(/field/),rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving StateAddField()'
   write(debug,*)
#endif
   return

   end subroutine StateAddField
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: StateCompleteConnectedField -
!
! !INTERFACE:
   subroutine StateCompleteConnectedField(state,name,kwe,&
                                          farray2d,farray3d,p2dr,p3dr,frc)
!
! !DESCRIPTION:
!  Links allocated arrays to empty field in state or removes unconnected
!  field. Returns pointer to array (also for an already completed field).
!  Note, that 3D fields are expected to contain k=0 layer!!!
!  priority order: 1) associated Ptr to ESMF_KIND_R8 array
!                  2) forced allocation
!                  3) Ptr to kind-matching farray
!                  4) allocation
!
! !USES:
   use domain    ,only: imin,jmin,imax,jmax,kmax
   use initialise,only: runtype
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*),intent(in)          :: name
   logical         ,intent(in),optional :: kwe !keyword-enforcer
   REALTYPE,target ,intent(in),optional :: farray2d(I2DFIELD)
   REALTYPE,target ,intent(in),optional :: farray3d(I3DFIELD)
   logical         ,intent(in),optional :: frc
!
! !INPUT/OUTPUT PARAMETERS:
   type(ESMF_State)                     :: state
   real(ESMF_KIND_R8),pointer ,optional :: p2dr(:,:),p3dr(:,:,:)
!
! !REVISION HISTORY:
!  Original Author(s): Knut Klingbeil
!
! !LOCAL VARIABLES
   type(ESMF_Field)            :: field
   type(ESMF_FieldStatus_Flag) :: status
   type(ESMF_Grid)             :: grid
   type(ESMF_StaggerLoc)       :: staggerloc
   type(ESMF_StateItem_Flag)   :: itemType
   real(ESMF_KIND_R8),pointer  :: p2dr_(:,:),p3dr_(:,:,:)
   REALTYPE                    :: getmreal
   integer                     :: rc,dimCount,klen
   integer,target              :: elb(3),eub(3)
   logical                     :: abort,isPresent,isConnected,noKindMatch
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'StateCompleteConnectedField() # ',Ncall
#endif

   call ESMF_StateGet(state,trim(name),itemType,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!  skip non-existing fields
   if (itemType .eq. ESMF_STATEITEM_NOTFOUND) return

   call ESMF_StateGet(state,trim(name),field,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_AttributeGetAttPack(field,                                 &
                                 convention="NUOPC",purpose="Instance", &
                                 isPresent=isPresent,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   if (isPresent) then
      call ESMF_AttributeGet(field,"Connected",                        &
                             convention="NUOPC",purpose="Instance",    &
                             isPresent=isPresent,rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
   end if

   if (.not. isPresent) then
      call ESMF_LogWrite('hope for non-NUOPC field '//trim(name)//'',  &
                         ESMF_LOGMSG_WARNING,line=__LINE__,file=FILENAME)
      return
   end if

   isConnected = NUOPC_IsConnected(field,rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   if (.not. isConnected) then
      call ESMF_LogWrite('remove unconnected field '//trim(name)//'',  &
                         ESMF_LOGMSG_WARNING,line=__LINE__,file=FILENAME)
      call ESMF_StateRemove(state,(/trim(name)/),rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      return
   end if

   call ESMF_FieldGet(field,status=status,grid=grid,staggerloc=staggerloc,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_GridGet(grid,dimCount=dimCount,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   noKindMatch = ( kind(getmreal) .ne. ESMF_KIND_R8 )

   if (runtype .eq. 1) then
      klen = 1
   else
      klen = kmax
   end if

   if (dimCount .eq. 2) then
      if (status.eq.ESMF_FIELDSTATUS_COMPLETE .and. present(p2dr)) then
         call ESMF_FieldGet(field,farrayPtr=p2dr_,rc=rc)
         abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
         p2dr => p2dr_
         return
      end if
      p2dr_ => NULL()
      if (present(p2dr)) then
         if (associated(p2dr)) p2dr_(imin-HALO:,jmin-HALO:) => p2dr
      end if
      if (.not. associated(p2dr_) .and. present(frc)) then
         if (frc) allocate(p2dr_(I2DFIELD))
      end if
      if (.not. associated(p2dr_) .and. present(farray2d)) then
         if (.not. noKindMatch) p2dr_ => farray2d
      end if
      if (.not. associated(p2dr_)) allocate(p2dr_(I2DFIELD))
      if (present(p2dr)) p2dr => p2dr_

      call ESMF_GridGetFieldBounds(grid,staggerloc=staggerloc,         &
                                   totalLBound=elb(1:2),totalUBound=eub(1:2),rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!     in contrast to ESMF_ArrayCreate() no automatic determination of total[L|U]Width
#if 1
!     Note (KK): in former times ESMF_FieldCreateGridDataPtr<rank><type><kind>() failed
      call ESMF_FieldEmptyComplete(field,p2dr_,                             &
#else
!     internal call to ESMF_FieldCreateGridData<rank><type><kind>()
!     forced by indexflag argument.
      call ESMF_FieldEmptyComplete(field,p2dr_,ESMF_INDEX_DELOCAL,          &
#endif
                                   totalLWidth=int(elb(1:2)-lbound(p2dr_)), &
                                   totalUWidth=int(ubound(p2dr_)-eub(1:2)), &
                                   rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   else if (dimCount .eq. 3) then
      if (status.eq.ESMF_FIELDSTATUS_COMPLETE .and. present(p3dr)) then
         call ESMF_FieldGet(field,farrayPtr=p3dr_,rc=rc)
         abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
         p3dr => p3dr_
         return
      end if
      p3dr_ => NULL()
      if (present(p3dr)) then
         if (associated(p3dr)) p3dr_(imin-HALO:,jmin-HALO:,0:) => p3dr
      end if
      if (.not. associated(p3dr_) .and. present(frc)) then
         if (frc) allocate(p3dr_(I2DFIELD,0:klen))
      end if
      if (.not. associated(p3dr_) .and. present(farray3d)) then
         if (.not. noKindMatch) p3dr_ => farray3d
      end if
      if (.not. associated(p3dr_)) allocate(p3dr_(I2DFIELD,0:klen))
      if (present(p3dr)) p3dr => p3dr_

      call ESMF_GridGetFieldBounds(grid,staggerloc=staggerloc,            &
                                   totalLBound=elb,totalUBound=eub,rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!     in contrast to ESMF_ArrayCreate() no automatic determination of total[L|U]Width
#if 1
!     Note (KK): in former times ESMF_FieldCreateGridDataPtr<rank><type><kind>() failed
      call ESMF_FieldEmptyComplete(field,p3dr_,                        &
#else
!     internal call to ESMF_FieldCreateGridData<rank><type><kind>()
!     forced by indexflag argument.
      call ESMF_FieldEmptyComplete(field,p3dr_,ESMF_INDEX_DELOCAL,     &
#endif
                                   totalLWidth=int(elb-lbound(p3dr_)), &
                                   totalUWidth=int(ubound(p3dr_)-eub), &
                                   rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   end if

   call NUOPC_Realize(state, field, rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

#ifdef DEBUG
   write(debug,*) 'Leaving StateCompleteConnectedField()'
   write(debug,*)
#endif
   return

   end subroutine StateCompleteConnectedField
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: StateReadCompleteField -
!
! !INTERFACE:
   subroutine StateReadCompleteField(state,name,kwe,                   &
                                     farray2d,farray3d,p2dr,p3dr,frc,ign)
!
! !DESCRIPTION:
!  Gets pointer or values from field in state.
!  returns if KindMatch and unforced
!  frc tries to read also if KindMatch
!  ign skips missing fields that were forced to be read
!
! !USES:
   use domain, only: imin,jmin,imax,jmax,kmax
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type(ESMF_State),intent(in)           :: state
   character(len=*),intent(in)           :: name
   logical         ,intent(in) ,optional :: kwe !keyword-enforcer
   logical         ,intent(in) ,optional :: frc,ign
!
! !OUTPUT PARAMETERS:
   REALTYPE        ,intent(out),optional :: farray2d(I2DFIELD)
   REALTYPE        ,intent(out),optional :: farray3d(I3DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
   real(ESMF_KIND_R8),pointer  ,optional :: p2dr(:,:),p3dr(:,:,:)
!
! !REVISION HISTORY:
!  Original Author(s): Knut Klingbeil
!
! !LOCAL VARIABLES
   type(ESMF_Field)            :: field
   type(ESMF_FieldStatus_Flag) :: status
   type(ESMF_Grid)             :: grid
   type(ESMF_StateItem_Flag)   :: itemType
   real(ESMF_KIND_R8),pointer  :: p2dr_(:,:),p3dr_(:,:,:)
   REALTYPE                    :: getmreal
   integer                     :: rc,dimCount
   logical                     :: abort,frc_,noKindMatch,ign_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'StateReadCompleteField() # ',Ncall
#endif

   if (present(frc)) then
      frc_ = frc
   else
      frc_ = .false.
   end if

   if (present(ign)) then
      ign_ = ign
   else
      ign_ = .false.
   end if

   noKindMatch = ( kind(getmreal) .ne. ESMF_KIND_R8 )

   if (.not.frc_ .and. .not.noKindMatch) return

   call ESMF_StateGet(state,trim(name),itemType,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!  unconnected fields might have been removed
   if (itemType .eq. ESMF_STATEITEM_NOTFOUND) then
      if (ign_) return
      call ESMF_LogWrite('missing field '//trim(name)//'',             &
                         ESMF_LOGMSG_ERROR,line=__LINE__,file=FILENAME)
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
   end if

   call ESMF_StateGet(state,trim(name),field,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_FieldGet(field,status=status,grid=grid,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   if (status .ne. ESMF_FIELDSTATUS_COMPLETE) then
      if (ign_) return
      call ESMF_LogWrite('incomplete field '//trim(name)//'',          &
                         ESMF_LOGMSG_ERROR,line=__LINE__,file=FILENAME)
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
   end if

   call ESMF_GridGet(grid,dimCount=dimCount,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   if (dimCount .eq. 2) then
      call ESMF_FieldGet(field,farrayPtr=p2dr_,rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      if (present(p2dr)) then
         p2dr => p2dr_
         if (frc_) return
      end if
      if (present(farray2d)) farray2d = p2dr_
   else if (dimCount .eq. 3) then
      call ESMF_FieldGet(field,farrayPtr=p3dr,rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      if (present(p3dr)) then
         p3dr => p3dr_
         if (frc_) return
      end if
      if (present(farray3d)) farray3d = p3dr_
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving StateReadCompleteField()'
   write(debug,*)
#endif
   return

   end subroutine StateReadCompleteField
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: StateSetCompleteField -
!
! !INTERFACE:
   subroutine StateSetCompleteField(state,name,kwe,                    &
                                    farray2d,farray3d,p2dr,p3dr,frc)
!
! !DESCRIPTION:
!  Update field in state (either its pointer or values).
!  priority order: 1) return if KindMatch and unforced
!                  2) forced update of associated Ptr
!                  3) forced linking to farray if KindMatch
!                  4) copy of farray if noKindMatch
!
! !USES:
   use domain, only: imin,jmin,imax,jmax,kmax
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type(ESMF_State),intent(in)          :: state
   character(len=*),intent(in)          :: name
   logical         ,intent(in),optional :: kwe !keyword-enforcer
   REALTYPE,target ,intent(in),optional :: farray2d(I2DFIELD)
   REALTYPE,target ,intent(in),optional :: farray3d(I3DFIELD)

   logical         ,intent(in),optional :: frc
!
! !OUTPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
   real(ESMF_KIND_R8),pointer,optional   :: p2dr(:,:),p3dr(:,:,:)
!
! !REVISION HISTORY:
!  Original Author(s): Knut Klingbeil
!
! !LOCAL VARIABLES
   type(ESMF_Field)            :: field
   type(ESMF_FieldStatus_Flag) :: status
   type(ESMF_Grid)             :: grid
   type(ESMF_StateItem_Flag)   :: itemType
   real(ESMF_KIND_R8),pointer  :: p2dr_(:,:),p3dr_(:,:,:)
   REALTYPE                    :: getmreal
   integer                     :: rc,dimCount
   logical                     :: abort,frc_,noKindMatch,ptr_associated
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'StateSetCompleteField() # ',Ncall
#endif

   if (present(frc)) then
      frc_ = frc
   else
      frc_ = .false.
   end if

   noKindMatch = ( kind(getmreal) .ne. ESMF_KIND_R8 )

   if (.not.frc_ .and. .not.noKindMatch) return

   call ESMF_StateGet(state,trim(name),itemType,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

!  unconnected fields might have been removed
   if (itemType .eq. ESMF_STATEITEM_NOTFOUND) return

   call ESMF_StateGet(state,trim(name),field,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   call ESMF_FieldGet(field,status=status,grid=grid,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   if (status .ne. ESMF_FIELDSTATUS_COMPLETE) then
      !call ESMF_LogWrite('skip incomplete field '//trim(name)//'',  &
      !                   ESMF_LOGMSG_WARNING,line=__LINE__,file=FILENAME)
      !call ESMF_Finalize(endflag=ESMF_END_ABORT)
      return
   end if

   call ESMF_GridGet(grid,dimCount=dimCount,rc=rc)
   abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
   if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)

   if (dimCount .eq. 2) then
      call ESMF_FieldGet(field,farrayPtr=p2dr_,rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!     do not overwrite values of p2dr_ in case of KindMatch !!!
!     (might point to internally swapped GETM array)
      ptr_associated = .false.
      if (present(p2dr)) ptr_associated = associated(p2dr)
      if (frc_) then
         if (ptr_associated) then
            p2dr_(imin-HALO:,jmin-HALO:) => p2dr
         else if (present(farray2d)) then
            if (noKindMatch) then
               p2dr_ = farray2d
            else
               p2dr_ => farray2d
            end if
         end if
      else if (noKindMatch) then
         if (present(farray2d)) p2dr_ = farray2d
      end if
      if (present(p2dr)) p2dr => p2dr_
   else if (dimCount .eq. 3) then
      call ESMF_FieldGet(field,farrayPtr=p3dr_,rc=rc)
      abort = ESMF_LogFoundError(rc,line=__LINE__,file=FILENAME)
      if (abort) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      ptr_associated = .false.
      if (present(p3dr)) ptr_associated = associated(p3dr)
      if (frc_) then
         if (ptr_associated) then
            p3dr_(imin-HALO:,jmin-HALO:,0:) => p3dr
         else if (present(farray3d)) then
            if (noKindMatch) then
               p3dr_ = farray3d
            else
               p3dr_ => farray3d
            end if
         end if
      else if (noKindMatch) then
         if (present(farray3d)) p3dr_ = farray3d
      end if
      if (present(p3dr)) p3dr => p3dr_
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving StateSetCompleteField()'
   write(debug,*)
#endif
   return

   end subroutine StateSetCompleteField
!EOC
!-----------------------------------------------------------------------

   end module getm_esmf

!-----------------------------------------------------------------------
! Copyright (C) 2013 - Knut Klingbeil                                  !
!-----------------------------------------------------------------------
