#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  integration - Initialise the time and do the time loop
!
! !INTERFACE:
   module integration
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   integer                             :: MinN=1,MaxN=-1
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
! !ROUTINE: time_loop - the main loop of getm
!
! !INTERFACE:
   subroutine time_loop(runtype)
!
! !DESCRIPTION:
!  A wrapper that calls meteo\_forcing, integrate\_2d, integrate\_3d,
!  do\_getm\_bio and output
!  within a time loop.
!
! !USES:
   use time,     only: update_time,timestep
   use domain,   only: kmax
   use meteo,    only: do_meteo,tausx,tausy,airp,fwf_method,evap,precip
   use m2d,      only: integrate_2d
   use variables_2d, only: fwf,fwf_int
#ifndef NO_3D
   use m3d,      only: integrate_3d,M
#ifndef NO_BAROCLINIC
   use variables_3d, only: T
#endif
   use rivers,   only: do_rivers
#ifdef _FABM_
   use gotm_fabm, only: fabm_calc
   use getm_fabm, only: do_getm_fabm
#endif
#ifdef GETM_BIO
   use bio, only: bio_calc
   use getm_bio, only: do_getm_bio
#endif
#endif
#ifdef SPM
   use suspended_matter, only: spm_calc,do_spm
#endif
   use input,    only: do_input
   use output,   only: do_output
#ifdef TEST_NESTING
   use nesting,   only: nesting_file
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: runtype
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES
   logical                   :: do_3d=.false.
   integer                   :: n
   integer                   :: progress=100
   character(8)              :: d_
   character(10)             :: t_
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'time_loop() # ',Ncall
#endif

   STDERR LINE
   LEVEL1 'integrating....'
   STDERR LINE

   do n=MinN,MaxN

      if (progress .gt. 0 .and. mod(n,progress) .eq. 0) then
         call date_and_time(date=d_,time=t_)
         LEVEL1 t_(1:2),':',t_(3:4),':',t_(5:10),' n=',n
      end if

#ifndef NO_3D
      do_3d = (runtype .ge. 2 .and. mod(n,M) .eq. 0)
#endif
      call do_input(n,do_3d)
      if(runtype .le. 2) then
         call do_meteo(n)
#ifndef NO_3D
#ifndef NO_BAROCLINIC
      else
         call do_meteo(n,T(:,:,kmax))
#endif
#endif
      end if

      if (fwf_method .ge. 1) then
         fwf = evap+precip
#ifndef NO_3D
         fwf_int = fwf_int+timestep*fwf
#endif
      end if

#ifndef NO_BAROTROPIC
      call integrate_2d(runtype,n,tausx,tausy,airp)
#endif
#ifndef NO_3D
      call do_rivers(n,do_3d)
      if (do_3d) then
         call integrate_3d(runtype,n)
#ifdef SPM
         if (spm_calc) call do_spm()
#endif
#ifdef _FABM_
         if (fabm_calc) call do_getm_fabm(M*timestep)
#endif
#ifdef GETM_BIO
         if (bio_calc) call do_getm_bio(M*timestep)
#endif
#ifndef NO_3D
         if (fwf_method .ge. 1) then
            fwf_int = _ZERO_
         end if
#endif
      end if
#endif

#ifdef TEST_NESTING
      if (mod(n,80) .eq. 0) then
         call nesting_file(WRITING)
      end if
#endif
      call update_time(n)

      call do_output(runtype,n,timestep)
#ifdef DIAGNOSE
      call diagnose(n,MaxN,runtype)
#endif
   end do


#ifdef DEBUG
   write(debug,*) 'Leaving time_loop()'
   write(debug,*)
#endif
   return
   end subroutine time_loop
!EOC

!-----------------------------------------------------------------------

   end module integration

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
