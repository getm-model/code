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
   use exceptions
   use time, only: julianday,secondsofday
   use parameters, only: avmmol
   use domain, only: imin,imax,iextr,jmin,jmax,jextr,az,au,av,ax,H,HU,HV,min_depth,z0_method
   use domain, only: ilg,ihg,jlg,jhg
   use domain, only: ill,ihl,jll,jhl
   use domain, only: rigid_lid,openbdy
   use advection, only: init_advection,print_adv_settings,NOADV
   use m2d_general, only: bottom_friction,calc_uvex
   use les, only: les_mode,LES_MOMENTUM
   use halo_zones, only : update_2d_halo,wait_halo,H_TAG
   use variables_2d
   use bdy_2d, only: bdyfmt_2d,bdyramp_2d
   IMPLICIT NONE

! Temporary interface (should be read from module):
   interface
      subroutine get_2d_field(fn,varname,il,ih,jl,jh,f)
         character(len=*),intent(in)   :: fn,varname
         integer, intent(in)           :: il,ih,jl,jh
         REALTYPE, intent(out)         :: f(:,:)
      end subroutine get_2d_field
   end interface
!
! !PUBLIC DATA MEMBERS:
   logical                   :: no_2d
   logical                   :: have_boundaries
   integer                   :: vel_adv_split2d=0
   integer                   :: vel_adv_scheme=1
   logical                   :: deformCX=.false.,deformUV=.false.
   integer,parameter         :: NO_AM=0,AM_CONSTANT=1,AM_LES=2
   integer                   :: Am_method=NO_AM
   REALTYPE                  :: Am_const=1.8d-6
!  method for specifying horizontal numerical diffusion coefficient
!     (0=const, 1=from named nc-file)
   integer                   :: An_method=0
   REALTYPE                  :: An_const=-_ONE_
   character(LEN = PATH_MAX) :: An_file
   integer                   :: MM=1,residual=-1
   integer                   :: sealevel_check=0
   logical                   :: bdy2d=.false.
   character(len=PATH_MAX)   :: bdyfile_2d
   integer, parameter        :: comm_method=-1
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
   integer                   :: elev_method=1
   REALTYPE                  :: elev_const=_ZERO_
   character(LEN = PATH_MAX) :: elev_file='elev.nc'
   namelist /m2d/ &
          elev_method,elev_const,elev_file,              &
          MM,vel_adv_split2d,vel_adv_scheme,             &
          Am_method,Am_const,An_method,An_const,An_file, &
          residual,sealevel_check,                       &
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

   dtm = timestep

!  Read 2D-model specific things from the namelist.
   read(NAMLST,m2d)

!  Allocates memory for the public data members - if not static
   call init_variables_2d(runtype)
   call init_advection()

   LEVEL2 'Advection of depth-averaged velocities'
#ifdef NO_ADVECT
   if (vel_adv_scheme .ne. NOADV) then
      LEVEL2 "reset vel_adv_scheme= ",NOADV," because of"
      LEVEL2 "obsolete NO_ADVECT macro. Note that this"
      LEVEL2 "behaviour will be removed in the future."
      vel_adv_scheme = NOADV
   end if
#endif
   call print_adv_settings(vel_adv_split2d,vel_adv_scheme,_ZERO_)

   if (.not. hotstart) then
      select case (elev_method)
         case(1)
            LEVEL2 'setting initial surface elevation to ',real(elev_const)
            z = elev_const
         case(2)
            LEVEL2 'getting initial surface elevation from ',trim(elev_file)
            call get_2d_field(trim(elev_file),"elev",ilg,ihg,jlg,jhg,z(ill:ihl,jll:jhl))
!           Note (KK): we need halo update only for periodic domains
            call update_2d_halo(z,z,az,imin,jmin,imax,jmax,H_TAG)
            call wait_halo(H_TAG)
         case default
            stop 'init_2d(): invalid elev_method'
      end select

      where ( z .lt. -H+min_depth)
         z = -H+min_depth
      end where
      zo = z
      call depth_update()
   end if

#if defined(GETM_PARALLEL) || defined(NO_BAROTROPIC)
!   STDERR 'Not calling cfl_check() - GETM_PARALLEL or NO_BAROTROPIC'
!   call cfl_check()
#else
!  KK-TODO: why is cfl_check not in terms of D?
   if (.not. rigid_lid) call cfl_check()
#endif

   select case (Am_method)
      case(NO_AM)
         LEVEL2 'Am_method=0 -> horizontal momentum diffusion not included'
      case(AM_CONSTANT)
         LEVEL2 'Am_method=1 -> Using constant horizontal momentum diffusion'
         if (Am_const .lt. _ZERO_) then
           call getm_error("init_2d()", &
                           "Constant horizontal momentum diffusion <0");
         end if
         LEVEL3 real(Am_const)
         deformCX=.true.
      case(AM_LES)
         LEVEL2 'Am_method=2 -> using LES parameterisation'
         les_mode=LES_MOMENTUM
         deformCX=.true.
         deformUV=.true.
      case default
         call getm_error("init_2d()", &
                         "A non valid Am_method has been chosen");
   end select

   select case (An_method)
      case(0)
         LEVEL2 'An_method=0 -> horizontal numerical diffusion not included'
      case(1)
         LEVEL2 'An_method=1 -> Using constant horizontal numerical diffusion'
         if (An_const .lt. _ZERO_) then
            call getm_error("init_2d()", &
                            "Constant horizontal numerical diffusion <0");
         end if
      case(2)
         LEVEL2 'An_method=2 -> Using space varying horizontal numerical diffusion'
         LEVEL2 '..  will read An from An_file ',trim(An_file)

         allocate(AnC(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_2d: Error allocating memory (AnC)'
         AnC = _ZERO_

         call get_2d_field(trim(An_file),"An",ilg,ihg,jlg,jhg,AnC(ill:ihl,jll:jhl))

         if (MINVAL(AnC(imin:imax,jmin:jmax),mask=(az(imin:imax,jmin:jmax).ge.1)) .lt. _ZERO_) then
            call getm_error("init_2d()", &
                            "negative numerical diffusivity in An field");
         end if

!        Note (KK): halo update is only needed for periodic domains
         call update_2d_halo(AnC,AnC,az,imin,jmin,imax,jmax,H_TAG)
         call wait_halo(H_TAG)

         if (MAXVAL(AnC(imin-1:imax+1,jmin-1:jmax+1),mask=(az(imin-1:imax+1,jmin-1:jmax+1).ge.1)) .eq. _ZERO_) then
!           Note (BJB): If all An values are really zero, then we should not use An-smoothing at all...
!                       Note that smoothing may be on in other subdomains.
            LEVEL2 '  All An is zero for this (sub)domain - switching to An_method=0'
            An_method=0
         else
!           Note (KK): since a HALO update of AnX is not needed,
!                      the allocation of AnX can be done locally
!                      and dependent on the test above

            allocate(AnX(E2DFIELD),stat=rc)
            if (rc /= 0) stop 'init_2d: Error allocating memory (AnX)'

            ! Compute AnX (An in X-points) based on AnC and the X- and T- masks
            ! We loop over the X-points in the present domain.
            do j=jmin-1,jmax
               do i=imin-1,imax
                  if (ax(i,j) .ge. 1) then
                     AnX(i,j) = _QUART_*( AnC(i,j) + AnC(i+1,j) + AnC(i,j+1) + AnC(i+1,j+1) )
                  end if
               end do
            end do
         end if
      case default
         call getm_error("init_2d()", &
                         "A non valid An method has been chosen");
   end select

   if (.not. openbdy)  bdy2d=.false.

   if (rigid_lid) then
      if (bdy2d) then
         LEVEL2 'Reset bdy2d=F because of rigid lid'
         bdy2d=.false.
      end if
   else
      LEVEL2 'Open boundary=',bdy2d
      if (bdy2d) then
         if (hotstart .and. bdyramp_2d .gt. 0) then
             LEVEL2 'WARNING: hotstart is .true. AND bdyramp_2d .gt. 0'
             LEVEL2 'WARNING: .. be sure you know what you are doing ..'
         end if
         LEVEL2 TRIM(bdyfile_2d)
         LEVEL2 'Format=',bdyfmt_2d
      end if
      if (sealevel_check .eq. 0) then
         LEVEL2 'sealevel_check=0 --> NaN checks disabled'
      else if (sealevel_check .gt. 0) then
         LEVEL2 'sealevel_check>0 --> NaN values will result in error conditions'
      else
         LEVEL2 'sealevel_check<0 --> NaN values will result in warnings'
      end if
   end if

   if (deformCX) then

      allocate(dudxC(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (dudxC)'
      dudxC=_ZERO_

#ifndef SLICE_MODEL
      allocate(dvdyC(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (dvdyC)'
      dvdyC=_ZERO_
#endif

      allocate(shearX(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (shearX)'
      shearX=_ZERO_

      if (deformUV) then
         allocate(dudxV(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_2d: Error allocating memory (dudxV)'
         dudxV=_ZERO_

#ifndef SLICE_MODEL
         allocate(dvdyU(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_2d: Error allocating memory (dvdyU)'
         dvdyU=_ZERO_
#endif

         allocate(shearU(E2DFIELD),stat=rc)
         if (rc /= 0) stop 'init_2d: Error allocating memory (shearU)'
         shearU=_ZERO_
      end if

   end if

#ifdef SLICE_MODEL
!  Note (KK): sse=0,U=0,dyV=0,V set in 3d
   no_2d = rigid_lid
#else
!  Note (KK): sse=0,U=V=0
   no_2d = (rigid_lid .and. (imin.eq.iextr .or. jmin.eq.jextr))
#endif

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
! !IROUTINE: postinit_2d - re-initialise some 2D after hotstart read.
!
! !INTERFACE:
   subroutine postinit_2d(runtype,timestep,hotstart)
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
!  This routine provides possibility to reset/initialize 2D variables to
!  ensure that velocities are correctly set on land cells after read
!  of a hotstart file.
!
! !LOCAL VARIABLES:
   integer                   :: i,j, ischange
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'postinit_2d() # ',Ncall
#endif

   LEVEL1 'postinit_2d'

!
! It is possible that a user changes the land mask and reads an "old" hotstart file.
! In this case the "old" velocities will need to be zeroed out.
   if (hotstart) then

      ischange = 0
!     The first two loops are pure diagnostics, logging where changes will actually take place
!     (and if there is something to do at all, to be able to skip the second part)
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO
            if ( au(i,j).eq.0 .and. U(i,j).ne._ZERO_ ) then
               LEVEL3 'postinit_2d: Reset to mask(au), U=0 for i,j=',i,j
               ischange = 1
            end if
         end do
      end do
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO
            if ( av(i,j).eq.0 .and. V(i,j).ne._ZERO_ ) then
               LEVEL3 'postinit_2d: Reset to mask(av), V=0 for i,j=',i,j
               ischange = 1
            end if
         end do
      end do

!     The actual reset is below here - independent of the above diagnostics (except for the if)
      if (ischange.ne.0) then
         where (au .eq. 0)
            U    = _ZERO_
            Uint = _ZERO_
         end where
         where (av .eq. 0)
            V    = _ZERO_
            Vint = _ZERO_
         end where
!        This is probably not absolutely necessary:
         where (az .eq. 0)
            z  = _ZERO_
            zo = _ZERO_
         end where
      end if

      call depth_update()

   end if

   return
   end subroutine postinit_2d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: integrate_2d - sequence of calls to do 2D model integration
!
! !INTERFACE:
   subroutine integrate_2d(runtype,loop,tausx,tausy,airp)
   use getm_timers, only: tic, toc, TIM_INTEGR2D

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
   call tic(TIM_INTEGR2D)

   if (mod(loop-1,MM) .eq. 0) then        ! MacroMicro time step
      if (z0_method .ne. 0) then
         call bottom_friction(U,V,DU,DV,ru,rv)
      end if
   end if

   call calc_uvex(An_method,U,V,D,DU,DV)

   call toc(TIM_INTEGR2D)

   call momentum(loop,tausx,tausy,airp)

   if (rigid_lid) then
!     Note (KK): we need to solve Poisson equation to get final transports
!                that fulfill dxU+dyV=0
      stop 'integrate_2d(): Poisson solver for rigid lid computations not implemented yet!'
   end if

   if (runtype .gt. 1) then
      call tic(TIM_INTEGR2D)
      Uint=Uint+U
      Vint=Vint+V
      call toc(TIM_INTEGR2D)
   end if

   if (.not. rigid_lid) then
      call sealevel(loop)
      call depth_update()
   end if

   if(residual .gt. 0 .and. loop .ge. residual) then
      call tic(TIM_INTEGR2D)
      call do_residual(0)
      call toc(TIM_INTEGR2D)
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
