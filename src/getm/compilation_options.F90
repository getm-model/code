#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: compilation_options() - 
!
! !INTERFACE:
   subroutine compilation_options()
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
!EOP
!-----------------------------------------------------------------------
!BOC

#ifdef GETM_PARALLEL
   LEVEL1 'Compiled for parallel execution'
#else
   LEVEL1 'Compiled for serial execution'
#endif

#ifdef GETM_OMP
   LEVEL1 '   with OpenMP thread capability'
#else
   LEVEL1 '   without OpenMP thread capability'
#endif

!
#ifdef SLICE_MODEL
   LEVEL1 'SLICE_MODEL'
#endif
#ifdef NO_3D
   LEVEL1 'NO_3D'
#endif
#ifdef NO_BAROCLINIC
   LEVEL1 'NO_BAROCLINIC'
#endif
!
#ifdef FORTRAN90
   LEVEL1 'Fortran 90 compilation'
#endif
!
#ifdef FORTRAN95
   LEVEL1 'Fortran 95 compilation'
#endif
!
#ifdef PRODUCTION
   LEVEL1 'Production compilation'
#endif
!
#ifdef PROFILING
   LEVEL1 'Profiling is enabled'
#endif
!
#ifdef DEBUG
   LEVEL1 'Debugging enabled'
#endif
!
#ifdef STATIC
   LEVEL1 'Using STATIC memory allocation'
#else
   LEVEL1 'Using DYNAMIC memory allocation'
#endif
!
#ifdef SINGLE
   LEVEL1 'Using single precision'
#else
   LEVEL1 'Using double precision'
#endif
!
! Various tests
#ifdef CARTESIAN
   LEVEL1 'CARTESIAN'
#endif
#ifdef SPHERICAL
   LEVEL1 'SPHERICAL'
#endif
#ifdef CURVILINEAR
   LEVEL1 'CURVILINEAR'
#endif
#ifdef TURB_ADV
   LEVEL1 'TURB_ADV'
#endif
#ifdef NO_BOTTFRIC
   LEVEL1 'NO_BOTTFRIC'
#endif
#ifdef NO_ADVECT
   LEVEL1 'NO_ADVECT'
#endif
#ifdef NO_SLR
   LEVEL1 'NO_SLR'
#endif
#ifdef _SLR_NOCLIP_
   LEVEL1 '_SLR_NOCLIP_'
#endif
#ifdef CONSTANT_VISCOSITY
   LEVEL1 'CONSTANT_VISCOSITY'
#endif
#ifdef PARABOLIC_VISCOSITY
   LEVEL1 'PARABOLIC_VISCOSITY'
#endif
#ifdef MIN_VEL_DEPTH
   LEVEL1 'MIN_VEL_DEPTH'
#endif
#ifdef NEW_SS
   LEVEL1 'NEW_SS'
#endif
#ifdef NEW_CORI
   LEVEL1 'NEW_CORI'
#endif
#ifdef SMOOTH_BVF_HORI
   LEVEL1 'SMOOTH_BVF_HORI'
#endif
#ifdef _SMOOTH_BVF_VERT_
   LEVEL1 '_SMOOTH_BVF_VERT_'
#endif
#ifdef NONNEGSALT
   LEVEL1 'NONNEGSALT'
#endif
#ifdef USE_BREAKS
   LEVEL1 'USE_BREAKS'
#endif
#ifdef PRESS_GRAD_Z
   LEVEL1 'PRESS_GRAD_Z'
#endif
#ifdef ITERATE_VERT_ADV
   LEVEL1 'ITERATE_VERT_ADV'
#endif
#ifdef SUBSTR_INI_PRESS
   LEVEL1 'SUBSTR_INI_PRESS'
#endif
#ifdef SONG_WRIGHT
   LEVEL1 'SONG_WRIGHT'
#endif
#ifdef NO_TIMERS
   LEVEL1 'NO_TIMERS'
#endif
#ifdef OLD_WRONG_FLUXES
   LEVEL1 'OLD_WRONG_FLUXES'
#endif
#ifdef _WRITE_HALOS_
   LEVEL1 '_WRITE_HALOS_'
#endif
#ifdef _WRITE_HOT_HALOS_
   LEVEL1 '_WRITE_HOT_HALOS_'
#endif
#ifdef _READ_HOT_HALOS_
   LEVEL1 '_READ_HOT_HALOS_'
#endif
#ifdef GETM_BIO
   LEVEL1 'GETM_BIO'
#endif
#ifdef _FABM_
   LEVEL1 '_FABM_'
#endif
#ifdef _POINTER_REMAP_
   LEVEL1 '_POINTER_REMAP_'
#endif
#ifdef _NCDF_SAVE_DOUBLE_
   LEVEL1 '_NCDF_SAVE_DOUBLE_'
#endif
#ifdef _FLEXIBLE_OUTPUT_
   LEVEL1 '_FLEXIBLE_OUTPUT_'
#endif

   STDERR LINE

   return
   end subroutine compilation_options
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2018 - Knut Klingbeil (IOW)                            !
!-----------------------------------------------------------------------

