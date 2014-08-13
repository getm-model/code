#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: variables_2d - global variables for 2D model
!
! !INTERFACE:
   module variables_2d
!
! !DESCRIPTION:
!  This modules contains declarations for all variables related to 2D
!  hydrodynamical calculations. Information about the calculation domain
!  is included from the {\tt domain} module.
!  The module contains public subroutines to initialise and cleanup.
!  Depending whether the compiler option {STATIC} is set or not,
!  memory for 2D variables is statically or dynamically allocated, see
!  {\tt PUBLIC DATA MEMBERS}.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
   use waves , only: waves_method,NO_WAVES,WAVES_VF
   use domain, only: bottfric_method,rdrag
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   REALTYPE                            :: dtm
   REALTYPE,dimension(:,:),pointer     :: zo,z
   logical                             :: do_numerical_analyses_2d=.false.
   logical                             :: calc_taubmax=.false.

#ifdef STATIC
#include "static_2d.h"
#else
#include "dynamic_declarations_2d.h"
#endif

   REALTYPE,dimension(:,:),pointer     :: Uf=>NULL(),Vf=>NULL()
   REALTYPE,dimension(:,:),pointer     :: UEuler=>NULL(),VEuler=>NULL()
   REALTYPE,dimension(:,:),pointer     :: UEulerInt=>NULL(),VEulerInt=>NULL()
   REALTYPE,dimension(:,:),pointer     :: taubmax=>NULL()
!  the following fields will only be allocated if deformCX=.true.
   REALTYPE,dimension(:,:),allocatable :: dudxC,dvdyC
   REALTYPE,dimension(:,:),pointer     :: dudyX=>null()
   REALTYPE,dimension(:,:),pointer     :: dvdxX=>null()
   REALTYPE,dimension(:,:),allocatable :: shearX
!  the following fields will only be allocated if deformUV=.true.
   REALTYPE,dimension(:,:),allocatable :: dudxV,dvdyU,shearU
!  the followinf fields will only be allocated if An_method=2
   REALTYPE,dimension(:,:),allocatable :: AnC,AnX

   REALTYPE,dimension(:,:),pointer     :: numdis_2d=>null()
   REALTYPE,dimension(:,:),pointer     :: phydis_2d=>null()
   REALTYPE,dimension(:,:),allocatable :: numdis_2d_old

   integer                             :: size2d_field
   integer                             :: mem2d
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
! !IROUTINE: init_variables_2d - initialise 2D related stuff.
!
! !INTERFACE:
   subroutine init_variables_2d(runtype,no_2d)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Allocates memory (unless {\tt STATIC} is set) for 2D related fields,
!  by an include statement. Furthermore all public 2D variables are
!  initialised to zero. Those are listed in table \ref{table_2d_variables}
!  on page \pageref{table_2d_variables}.
!
!  \begin{table}[h]
!  \begin{center}
!  \begin{tabular}{lll}
! {\tt z } & sea surface elevation in T-point & [m] \\
! {\tt U } & $x$ component of transport in U-point & [m$^2$s$^{-1}$] \\
! {\tt DU } & water depth in U-point & [m] \\
! {\tt fU } & Coriolis term for $V$-equation in V-point & [m$^2$s$^{-2}$] \\
! {\tt SlUx } & slow term for $U$-equation in U-point & [m$^2$s$^{-2}$] \\
! {\tt Slru } &slow bottom friction for $U$-equation in U-point &
! [m$^2$s$^{-2}$]\\
! {\tt V } & $y$ component of transport in V-point & [m$^2$s$^{-1}$]\\
! {\tt DV } & water depth in V-point & [m] \\
! {\tt fV } &  Coriolis term for $U$-equation in U-point & [m$^2$s$^{-2}$]\\
! {\tt SlVx } & slow term for $V$-equation in V-point & [m$^2$s$^{-2}$] \\
! {\tt Slrv } &slow bottom friction for $V$-equation in V-point &
! [m$^2$s$^{-2}$]\\
! {\tt Uint } & $x$-component of mean transport in U-point & [m$^2$s$^{-1}$]\\
! {\tt Vint } & $y$-component of mean transport in V-point & [m$^2$s$^{-1}$]\\
! {\tt UEx } & sum of explicit terms for for $U$-equation in U-point & [m$^2$s$^{-2}$]\\
! {\tt VEx } &sum of explicit terms for for $V$-equation in V-point & [m$^2$s$^{-2}$]\\
! {\tt ru } & bottom friction for $U$-equation in U-point & [m$^2$s$^{-2}$]\\
! {\tt rv } &bottom friction for $V$-equation in V-point & [m$^2$s$^{-2}$]\\
! {\tt res\_du } & residual depth in U-point & [m]\\
! {\tt res\_u } & $x$-component of residual transport in U-point &
! [m$^2$s$^{-1}$]\\
! {\tt res\_dv } & residual depth in V-point & [m] \\
! {\tt res\_v } &$y$-component of residual transport in V-point &
! [m$^2$s$^{-1}$]\\
! \end{tabular}
! \caption{Public 2D variables.}
! \label{table_2d_variables}
! \end{center}
! \end{table}
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: runtype
   logical, intent(in)                 :: no_2d
!
! !LOCAL VARIABLES:
   integer                   :: rc
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_variables_2d() # ',Ncall
#endif

   LEVEL2 'init_variables_2d'
   size2d_field=((imax+HALO)-(imin+HALO)+1)*((jmax+HALO)-(jmin+HALO)+1)
   mem2d=n2d_fields*size2d_field*REAL_SIZE

!  Allocates memory for the public data members - if not static
#ifndef STATIC
#include "dynamic_allocations_2d.h"
#endif

   zo => t_zo ; z => t_z

#ifdef USE_BREAKS
   break_stat = 0
#endif

   Uf        => U    ; Vf        => V
   UEuler    => U    ; VEuler    => V
   UEulerInt => Uint ; VEulerInt => Vint


   if (waves_method .ne. NO_WAVES) then

      allocate(UEuler(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (UEuler)'
      allocate(VEuler(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (VEuler)'

      if (runtype .gt. 1) then
         allocate(UEulerInt(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_2d: Error allocating memory (UEulerInt)'
         UEulerInt = _ZERO_
         allocate(VEulerInt(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_2d: Error allocating memory (VEulerInt)'
         VEulerInt = _ZERO_
      end if

      if (waves_method .eq. WAVES_VF) then
         Uf => UEuler ; Vf => VEuler
      end if

   end if


   z  = _ZERO_; zo =_ZERO_
   D = _ZERO_ ; Dvel = _ZERO_
   U = _ZERO_; DU = _ZERO_; Uint = _ZERO_; UEx = _ZERO_
   V = _ZERO_; DV = _ZERO_; Vint = _ZERO_; VEx = _ZERO_


   if (bottfric_method .eq. 1) then
      ru = rdrag
      rv = rdrag
   else
      ru = _ZERO_
      rv = _ZERO_
   end if

   res_du = _ZERO_; res_u = _ZERO_
   res_dv = _ZERO_; res_v = _ZERO_

   SlUx=_ZERO_; Slru=_ZERO_
   SlVx=_ZERO_; Slrv=_ZERO_

   fwf     = _ZERO_
   fwf_int = _ZERO_

   EWbdy=_ZERO_
   ENbdy=_ZERO_
   EEbdy=_ZERO_
   ESbdy=_ZERO_

#ifdef DEBUG
   write(debug,*) 'Leaving init_variables_2d()'
   write(debug,*)
#endif
   return
   end subroutine init_variables_2d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clean_variables_2d - cleanup after 2D run.
!
! !INTERFACE:
   subroutine clean_variables_2d()
   IMPLICIT NONE
!
! !DESCRIPTION:
!  This routine is currently empty.
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'clean_variables_2d() # ',Ncall
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving clean_variables_2d()'
   write(debug,*)
#endif
   return
   end subroutine clean_variables_2d
!EOC

!-----------------------------------------------------------------------

   end module variables_2d

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
