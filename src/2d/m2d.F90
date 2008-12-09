!$Id: m2d.F90,v 1.25 2008-12-09 00:31:57 kb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: m2d - depth integrated hydrodynamical model (2D)
!
! !INTERFACE:
   module m2d
!
! !DESCRIPTION:
!  This module contains declarations for all variables related to 2D
!  hydrodynamical calculations. Information about the calculation domain
!  is included from the {\tt domain} module.
!  The module contains public subroutines for initialisation, integration
!  and clean up of the 2D model component.
!  The actual calculation routines are called in {\tt integrate\_2d}
!  and are linked
!  in from the library {\tt lib2d.a}.
!
! !USES:
   use time, only: julianday,secondsofday
   use parameters, only: avmmol
   use domain, only: imin,imax,jmin,jmax,az,au,av,H,HU,HV,min_depth
   use domain, only: openbdy,z0_method,z0_const,z0
   use halo_zones, only : update_2d_halo,wait_halo
   use halo_zones, only : U_TAG,V_TAG,H_TAG
   use variables_2d
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   logical                   :: have_boundaries
   REALTYPE                  :: dtm,Am=-_ONE_,An=-_ONE_
   integer                   :: MM=1,residual=-1
   integer                   :: sealevel_check=0
   logical                   :: bdy2d=.false.
   integer                   :: bdyfmt_2d,bdytype,bdyramp_2d=-1
   character(len=PATH_MAX)   :: bdyfile_2d
   REAL_4B                   :: bdy_data(1500)
   REAL_4B                   :: bdy_data_u(1500)
   REAL_4B                   :: bdy_data_v(1500)
   REAL_4B, allocatable      :: bdy_times(:)
   integer, parameter        :: comm_method=-1
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_2d - initialise 2D related stuff.
!
! !INTERFACE:
   subroutine init_2d(runtype,timestep,hotstart)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: runtype
   REALTYPE, intent(in)                :: timestep
   logical, intent(in)                 :: hotstart
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Here, the {\tt m2d} namelist is read from {\tt getm.inp}, and the check
!  for the fulfilment of the CFL criterium for shallow water theory
!  {\tt cfl\_check} is called. A major part of this subroutine deals
!  then with the setting of local bathymetry values and initial surface
!  elevations in $u$- and $v$-points, also by calls to the subroutines 
!  {\tt uv\_depths} and {\tt depth\_update}.
!
! !LOCAL VARIABLES:
   integer                   :: rc
   integer                   :: i,j
   integer                   :: vel_depth_method=0
   namelist /m2d/ &
          MM,vel_depth_method,Am,An,residual,sealevel_check, &
          bdy2d,bdyfmt_2d,bdyramp_2d,bdyfile_2d
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_2d() # ',Ncall
#endif

   LEVEL1 'init_2d'

!  Read 2D-model specific things from the namelist.
   read(NAMLST,m2d)

   dtm = timestep

!  Allocates memory for the public data members - if not static
   call init_variables_2d(runtype)

#ifdef PARALLEL
!   STDERR 'Not calling cfl_check() - PARALLEL'
   call cfl_check()
#else
   call cfl_check()
#endif

   if (Am .lt. _ZERO_) then
      LEVEL2 'Am < 0 --> horizontal momentum diffusion not included'
   end if
   if (An .lt. _ZERO_) then
      LEVEL2 'An < 0 --> numerical momentum diffusion not included'
   end if

   if (sealevel_check .eq. 0) then
      LEVEL2 'sealevel_check=0 --> NaN checks disabled'
   else if (sealevel_check .gt. 0) then
      LEVEL2 'sealevel_check>0 --> NaN values will result in error conditions'
   else
      LEVEL2 'sealevel_check<0 --> NaN values will result in warnings'
   end if

   if (.not. openbdy)  bdy2d=.false.
   LEVEL2 'Open boundary=',bdy2d
   if (bdy2d) then
      if (hotstart .and. bdyramp_2d .gt. 0) then
          LEVEL2 'WARNING: hotstart is .true. AND bdyramp_2d .gt. 0'
          LEVEL2 'WARNING: .. be sure you know what you are doing ..'
      end if
      LEVEL2 TRIM(bdyfile_2d)
      LEVEL2 'Format=',bdyfmt_2d
   end if

   call uv_depths(vel_depth_method)

   where ( -H+min_depth .gt. _ZERO_ )
      z = -H+min_depth
   end where
   zo=z

   where (-HU+min_depth .gt. _ZERO_ )
      zu = -HU+min_depth
   end where

   where (-HV+min_depth .gt. _ZERO_ )
      zv = -HV+min_depth
   end where

!  bottom roughness
   if (z0_method .eq. 0) then
      zub0 = z0_const
      zvb0 = z0_const
   end if
   if (z0_method .eq. 1) then
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO-1
           if (au(i,j) .gt. 0) zub0(i,j) = 0.5*(z0(i,j)+z0(i+1,j))
         end do
      end do
      do j=jmin-HALO,jmax+HALO-1
         do i=imin-HALO,imax+HALO
           if (av(i,j) .gt. 0) zvb0(i,j) = 0.5*(z0(i,j)+z0(i,j+1))
         end do
      end do
   end if
   zub=zub0
   zvb=zvb0

   call depth_update()

#ifdef DEBUG
   write(debug,*) 'Leaving init_2d()'
   write(debug,*)
#endif
   return
   end subroutine init_2d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: integrate_2d - sequence of calls to do 2D model integration
!
! !INTERFACE:
   subroutine integrate_2d(runtype,loop,tausx,tausy,airp)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: runtype,loop
   REALTYPE, intent(in)                :: tausx(E2DFIELD)
   REALTYPE, intent(in)                :: tausy(E2DFIELD)
   REALTYPE, intent(in)                :: airp(E2DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Here, all 2D related subroutines are called. The major calls and their
!  meaning are:
!
!  \vspace{0.5cm}
!
!  \begin{tabular}{ll}
!  {\tt call update\_2d\_bdy} & read in new lateral boundary conditions \\
!  {\tt call bottom\_friction} & update bottom friction\\
!  {\tt call uv\_advect} & calculate 2D advection terms\\
!  {\tt call uv\_diffusion} & calculate 2D  diffusion terms\\
!  {\tt call momentum} & iterate 2D momemtum equations\\
!  {\tt call sealevel} & update sea surface elevation\\
!  {\tt call depth\_update} & update water depths\\
!  {\tt call do\_residual} & calculate intermdediate values for residual currents
!  \end{tabular}
!
!  \vspace{0.5cm}
!
!  It should be noted that some of these calls may be excluded for certain
!  compiler options set in the {\tt Makefile} of the application.
!
! !LOCAL VARIABLES:
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'integrate_2d() # ',Ncall
#endif

   if (mod(loop-1,MM) .eq. 0) then        ! MacroMicro time step
#ifndef NO_BOTTFRIC
      call bottom_friction(runtype)
#endif
   end if
   UEx=_ZERO_ ; VEx=_ZERO_
#ifdef NO_ADVECT
   STDERR 'NO_ADVECT 2D'
#else
#ifndef UV_ADV_DIRECT
   call uv_advect()
   if (Am .gt. _ZERO_ .or. An .gt. _ZERO_) then
      call uv_diffusion(Am,An) ! Has to be called after uv_advect.
   end if
   call mirror_bdy_2d(UEx,U_TAG)
   call mirror_bdy_2d(VEx,V_TAG)
#endif
#endif
   call momentum(loop,tausx,tausy,airp)
   if (runtype .gt. 1) then
      Uint=Uint+U
      Vint=Vint+V
   end if
   if (have_boundaries) call update_2d_bdy(loop,bdyramp_2d)
   call sealevel()
   call depth_update()

   if(residual .gt. 0 .and. loop .ge. residual) then
      call do_residual(0)
   end if

#ifdef DEBUG
     write(debug,*) 'Leaving integrate_2d()'
     write(debug,*)
#endif
   return
   end subroutine integrate_2d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clean_2d - cleanup after 2D run.
!
! !INTERFACE:
   subroutine clean_2d()
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This routine executes a final call to {\tt do\_residual} where the residual
!  current calculations are finished.
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'clean_2d() # ',Ncall
#endif

   if(residual .gt. 0) then
      call do_residual(1)
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving clean_2d()'
   write(debug,*)
#endif
   return
   end subroutine clean_2d
!EOC

!-----------------------------------------------------------------------

   end module m2d

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
