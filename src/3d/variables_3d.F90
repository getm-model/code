!$Id: variables_3d.F90,v 1.8 2006-02-10 22:41:56 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: variables_3d - global 3D related variables
!
! !INTERFACE:
   module variables_3d
!
! !DESCRIPTION:
!  This modules contains declarations for all variables related to 3D
!  hydrodynamical calculations. Information about the calculation domain
!  is included from the {\tt domain} module.
!  The variables are either statically defined in {\tt static\_3d.h} or
!  dynamically allocated in {\tt dynamic\_declarations\_3d.h}.
!  The variables which need to be declared have the following dimensions,
!  units and meanings:
!
! \vspace{0.5cm}
! \begin{supertabular}{llll}
! {\tt kmin} & 2D & [-] & lowest index in T-point \\
! {\tt kumin} & 2D &[-]  & lowest index in U-point \\
! {\tt kvmin} & 2D &[-]  & lowest index in V-point \\
! {\tt kmin\_pmz} & 2D &[-]  & lowest index in T-point (poor man's
! $z$-coordinate)\\
! {\tt kumin\_pmz} & 2D &[-]  & lowest index in U-point (poor man's
! $z$-coordinate)\\
! {\tt kvmin\_pmz} & 2D &[-]  & lowest index in V-point (poor man's
! $z$-coordinate)\\
! {\tt uu} & 3D & [m$^2$s$^{-1}$] & layer integrated $u$ transport
! $p_k$\\
! {\tt vv} & 3D & [m$^2$s$^{-1}$] & layer integrated $v$ transport
! $q_k$\\
! {\tt ww} & 3D & [m\,s$^{-1}$] & grid-related vertical velocity
! $\bar w_k$\\
! {\tt ho} & 3D & [m] & old layer height in T-point \\
! {\tt hn} & 3D & [m]& new layer height in T-point \\
! {\tt huo} & 3D &[m]& old layer height in U-point \\
! {\tt hun} & 3D & [m]& new layer height in U-point \\
! {\tt hvo} & 3D & [m]& old layer height in V-point \\
! {\tt hvn} & 3D & [m]& new layer height in V-point \\
! {\tt hcc} & 3D &[-] & hydrostatic consistency index in T-points\\
! {\tt uuEx} & 3D & [m$^2$s$^{-2}$] & sum of advection and 
! diffusion for $u$-equation\\
! {\tt vvEx} & 3D &  [m$^2$s$^{-2}$]& sum of advection and 
! diffusion for $v$-equation\\
! {\tt num} & 3D &  [m$^2$s$^{-1}$]& eddy viscosity on $w$-points
! $\nu_t$\\
! {\tt nuh} & 3D &  [m$^2$s$^{-1}$]& eddy diffusivity on $w$-points $\nu'_t$\\
! {\tt tke} & 3D &  [m$^2$s$^{-2}$]& turbulent kinetic energy $k$\\
! {\tt eps} & 3D &  [m$^2$s$^{-3}$]& turbulent dissipation rate
! $\eps$ \\
! {\tt SS} & 3D & [s$^{-2}$]& shear-frequency squared $M^2$ \\
! {\tt NN} & 3D &  [s$^{-2}$]& Brunt-V\"ais\"al\"a frequency squared$N^2$ \\
! {\tt S} & 3D & [psu] & salinity $S$ \\
! {\tt T} & 3D & [$^{\circ}$C]& potential temperature $\theta$ \\
! {\tt rho} & 3D & [m\,s$^{-2}$]& buoyancy $b$ \\
! {\tt idpdx} & 3D & [m$^2$s$^{-2}$] & $x$-component of internal
! pressure gradient \\
! {\tt idpdy} & 3D & [m$^2$s$^{-2}$]& $y$-component of internal
! pressure gradient\\
! {\tt spm} & 3D & [kg\,m$^{-3}$] & suspended matter concentration \\
! {\tt spm\_ws} & 3D & [m\,s$^{-1}$] & settling velocity of
! suspended matter \\
! {\tt spm\_pool} & 2D & [kg\,m$^{-2}$] & bottom pool of suspended
! matter\\
! {\tt uadv} & 3D & [m\,s$^{-1}$] & interpolated $x$-component of
! momentum advection velocity \\
! {\tt vadv} & 3D &  [m\,s$^{-1}$]& interpolated $y$-component of 
! momentum advection velocity \\
! {\tt wadv} & 3D &  [m\,s$^{-1}$]& interpolated  vertical component of 
! momentum advection velocity \\
! {\tt huadv} & 3D &[m] & interpolated height of advective flux
! layer ($x$-component) \\
! {\tt hvadv} & 3D &[m] & interpolated height of advective flux
! layer ($y$-component) \\
! {\tt hoadv} & 3D &[m] & old height of advective finite volume cell
! \\
! {\tt hnadv} & 3D &[m] & new height of advective finite volume
! cell\\
! {\tt sseo} & 2D & [m]& sea surface elevation before macro time
! step (T-point)\\
! {\tt ssen} & 2D & [m]& sea surface elevation after macro time
! step (T-point)\\
! {\tt ssuo} & 2D & [m]& sea surface elevation before macro time
! step (U-point)\\
! {\tt ssun} & 2D & [m]&sea surface elevation after macro time step
! (U-point)\\
! {\tt ssvo} & 2D & [m]& sea surface elevation before macro time
! step (V-point)\\
! {\tt ssvn} & 2D & [m]& sea surface elevation after macro time
! step (V-point)\\
! {\tt rru} & 2D & [m\,s$^{-1}$]&drag coefficient times curret speed
! in U-point\\
! {\tt rrv} & 2D & [m\,s$^{-1}$]&drag coefficient times curret speed
! in V-point\\
! {\tt taus} & 2D & [m$^2$s$^{-2}$]& normalised surface stress
! (T-point) \\
! {\tt taub} & 2D & [m$^2$s$^{-2}$]& normalised bottom stress
! (T-point) \\
! \end{supertabular}
!
! \vspace{0.5cm}
!
! It should be noted that depending on compiler options and runtype not
! all these variables are defined.
!
! The module contains public subroutines to initialise (see
! {\tt init\_variables\_3d}) and cleanup (see {\tt clean\_variables\_3d}).
!
! !USES:
   use domain,     only: iimin,iimax,jjmin,jjmax,kmax
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   REALTYPE                            :: dt,cnpar=0.9
   REALTYPE                            :: avmback=_ZERO_,avhback=_ZERO_
   character(len=64)                   :: adv_schemes(7)
!
#ifdef STATIC
#include "static_3d.h"
#else
#include "dynamic_declarations_3d.h"
#endif

#ifdef GETM_BIO
   REALTYPE, allocatable               :: cc3d(:,:,:,:)
   REALTYPE, allocatable               :: ws3d(:,:,:,:)
#endif
   integer                             :: size3d_field
   integer                             :: mem3d
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: variables_3d.F90,v $
!  Revision 1.8  2006-02-10 22:41:56  hb
!  Source code documentation extended
!
!  Revision 1.7  2005-09-23 11:27:10  kbk
!  support for biology via GOTMs biology modules
!
!  Revision 1.6  2004/01/06 15:04:00  kbk
!  FCT advection + split of advection_3d.F90 + extra adv. input checks
!
!  Revision 1.5  2003/12/16 15:58:54  kbk
!  back ground viscosity and diffusivity (manuel)
!
!  Revision 1.4  2003/04/23 12:16:34  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.3  2003/04/07 16:22:31  kbk
!  initialise variables
!
!  Revision 1.1.1.1  2002/05/02 14:00:58  gotm
!  recovering after CVS crash
!
!  Revision 1.6  2001/09/19 13:07:00  bbh
!  Moved advection related 3D fields to global allocation
!
!  Revision 1.5  2001/09/01 17:10:25  bbh
!  Vertical coordinate definition now specified via namelist
!
!  Revision 1.4  2001/08/27 11:51:45  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.3  2001/05/21 13:07:19  bbh
!  dt and cnpar is in variables_3d.F90
!
!  Revision 1.2  2001/05/18 08:25:52  bbh
!  Added zooming variables
!
!  Revision 1.1  2001/05/03 19:31:56  bbh
!  3D variables seperated from m3d
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_variables_3d - initialise 3D related stuff
! \label{sec-init-variables}
!
! !INTERFACE:
   subroutine init_variables_3d(runtype)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: runtype
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Dynamic allocation of memory for 3D related fields via
!  {\tt dynamic\_allocations\_3d.h} (unless the compiler option
!  {\tt STATIC} is set). Furthermore, most variables are initialised here.
!
! !REVISION HISTORY:
!
!  See log for the module.
!
! !LOCAL VARIABLES:
   integer                   :: rc
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_variables_3d() # ',Ncall
#endif

   LEVEL2 'init_variables_3d'
   size3d_field=((iimax+HALO)-(iimin+HALO)+1)*        &
                ((jjmax+HALO)-(jjmin+HALO)+1)*(kmax+1)
   mem3d=n3d_fields*size3d_field*REAL_SIZE

!  Allocates memory for the public data members - if not static
#ifndef STATIC
#include "dynamic_allocations_3d.h"
#endif

   hn = _ZERO_ ; hun = _ZERO_ ; hvn = _ZERO_
   uu = _ZERO_ ; vv = _ZERO_ ; ww = _ZERO_
   ssen = _ZERO_ ; ssun = _ZERO_ ; ssvn = _ZERO_
   rru= _ZERO_ ; rrv= _ZERO_
   uuEx= _ZERO_ ; vvEx= _ZERO_
   tke=1.e-10 ; eps=1.e-10

   light=_ONE_

#ifdef UV_TVD
   uadv = _ZERO_ ; vadv = _ZERO_ ; wadv = _ZERO_
   hnadv = _ZERO_ ; hoadv = _ZERO_
   huadv = _ZERO_ ; hvadv = _ZERO_
#endif

   adv_schemes(1) = "3D first-order upstream advection"
   adv_schemes(2) = "upstream advection (first-order, monotone)"
   adv_schemes(3) = "P2-PDM advection (third-order, non-monotone)"
   adv_schemes(4) = "TVD-Superbee advection (second-order, monotone)"
   adv_schemes(5) = "TVD-MUSCL advection (second-order, monotone)"
   adv_schemes(6) = "TVD-P2-PDM advection (third-order, monotone)"
   adv_schemes(7) = "2D-FCT advection"

#ifdef DEBUG
   write(debug,*) 'Leaving init_variables_3d()'
   write(debug,*)
#endif
   return
   end subroutine init_variables_3d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clean_variables_3d - cleanup after 3D run.
!
! !INTERFACE:
   subroutine clean_variables_3d()
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This routine cleans up after a 3D integrationby doing nothing so far.
!
! !REVISION HISTORY:
!  See log for the module.
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'clean_3d() # ',Ncall
#endif

! Deallocates memory for the public data members

#ifdef DEBUG
     write(debug,*) 'Leaving clean_variables_3d()'
     write(debug,*)
#endif
   return
   end subroutine clean_variables_3d
!EOC

!-----------------------------------------------------------------------

   end module variables_3d

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
