!$Id: variables_2d.F90,v 1.5 2006-02-04 11:21:52 hb Exp $
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
   use domain, only: imin,imax,jmin,jmax,H,HU,HV,min_depth
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
#ifdef STATIC
#include "static_2d.h"
#else
#include "dynamic_declarations_2d.h"
#endif
   integer                             :: size2d_field
   integer                             :: mem2d
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: variables_2d.F90,v $
!  Revision 1.5  2006-02-04 11:21:52  hb
!  Source code documentation extended
!
!  Revision 1.4  2003-04-23 12:09:44  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.3  2003/04/07 15:50:20  kbk
!  initialise variables
!
!  Revision 1.1.1.1  2002/05/02 14:00:47  gotm
!  recovering after CVS crash
!
!  Revision 1.1  2001/05/03 19:30:41  bbh
!  2D variables seperated from m2d
!
! !LOCAL VARIABLES:
   integer                   :: rc
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
   subroutine init_variables_2d(runtype)
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
! {\tt zu } & sea surface elevation in U-point & [m]\\ 
! {\tt zv } & sea surface elevation in V-point & [m]\\ 
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
! {\tt surfdiv } &divergence of surface currents in T-point & [s$^{-1}$]\\ 
! \end{tabular}
! \caption{Public 2D variables.} 
! \label{table_2d_variables}
! \end{center}
! \end{table}
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: runtype
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!
!  See log for module.
!
! !LOCAL VARIABLES:
!
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

   z = _ZERO_; zu = _ZERO_; zv = _ZERO_
   U = _ZERO_; DU = _ZERO_; fU = _ZERO_; SlUx = _ZERO_; Slru = _ZERO_
   V = _ZERO_; DV = _ZERO_; fV = _ZERO_; SlVx = _ZERO_; Slrv = _ZERO_

   Uint = _ZERO_; Vint = _ZERO_
   UEx = _ZERO_; VEx = _ZERO_
   ru = _ZERO_; rv = _ZERO_
   res_du = _ZERO_; res_u = _ZERO_; res_dv = _ZERO_; res_v =  _ZERO_
   surfdiv = _ZERO_

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
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This routine is currently empty.
!
! !REVISION HISTORY:
!  See log for module.
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
