#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: internal_pressure
!
! !INTERFACE:
   module internal_pressure
!
! !DESCRIPTION:
!
! In GETM, various methods are provided for the calculation of the
! internal pressure gradients terms in $x$- and $y$-direction.
! These terms which appear as layer-integrated terms in the
! equations for the layer-integrated momentum are for the
! eastward momentum $p_k$ (see equation (\ref{uEqvi})):
!
! \begin{equation}
! h_k\left(\frac12h_N(\partial^*_xb)_N
! +\sum_{j=k}^{N-1}\frac12(h_j+h_{j+1})(\partial^*_xb)_j
! \right)
! \end{equation}
!
! and for the northward layer-integrated momentum $q_k$
! (see equation (\ref{vEqvi})):
!
! \begin{equation}
! h_k\left(\frac12h_N(\partial^*_yb)_N
! +\sum_{j=k}^{N-1}\frac12(h_j+h_{j+1})(\partial^*_yb)_j
! \right)
! \end{equation}
!
! The major problem is how to calculate the horizontal (with respect
! to isogeopotentials) buoyancy gradients $\partial^*_xb$ and $\partial^*_yb$,
! which need to be defined at the interfaces positioned vertically
! between two velocity points.
!
! The methods for calculating the internal pressure gradient included in
! GETM are currently:
!
! \begin{enumerate}
! \item Method by \cite{MELLORea94}, see routine {\tt ip\_blumberg\_mellor}
! \item Modified \cite{MELLORea94} method, exact for linear density profiles
!       with $z$-dependence only, see routine {\tt ip\_blumberg\_mellor\_lin}
! \item Calculation by mean of linear interpolation to $z$-levels, see routine
!       {\tt ip\_z\_interpol}
! \item Method by \cite{SONG98}, see routine {\tt ip\_song\_wright}
! \item Method by \cite{CHUea03}, see routine {\tt ip\_chu\_fan}
! \item Method by \cite{SHCHEPETKINea03}, see routine {\tt
! ip\_shchepetkin\_mcwilliams}
! \item Method by \cite{STELLINGea94}, see routine {\tt
! ip\_stelling\_vankester.F90}
!  \end{enumerate}
!
! It is possible, by setting the compiler option {\tt SUBSTR\_INI\_PRESS},
! to substract the initial pressure gradient from all pressure
! gradients. This is only advisable for strong stratification
! without any initial internal pressure gradients. In this case
! any non-zero values of the resulting numerical initial pressure gradient
! are due to discretisation errors.
!
! !USES:
   use exceptions
   use domain, only: imin,imax,jmin,jmax,kmax,az,au,av,H,HU,HV
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxu,dyv
#else
   use domain, only: dx,dy
#endif
#ifndef NO_BAROCLINIC
   use variables_3d, only: buoy_eos=>buoy
#endif
   use variables_3d, only: kmin,hn,ssen,minus_bnh
   use variables_3d, only: idpdx_m3d=>idpdx,idpdy_m3d=>idpdy
   use variables_3d, only: idpdx_hs,idpdy_hs,idpdx_nh,idpdy_nh,idpdx_full,idpdy_full

   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   public init_internal_pressure, do_internal_pressure
   integer,public                           :: ip_method=1
   integer,public                           :: ip_ramp=-1
   logical,public                           :: ip_ramp_is_active=.false.
   logical,public                           :: calc_ipfull=.false.
   REALTYPE,dimension(:,:,:),pointer,public :: idpdx,idpdy,buoy
#ifdef STATIC
!  KK-TODO: this should become an automatic array in each ip routine
   REALTYPE,public                  :: zz(I3DFIELD)
#else
   REALTYPE,allocatable,public     :: zz(:,:,:)
#endif
!
! !PRIVATE DATA MEMBERS:
#ifdef SUBSTR_INI_PRESS
   REALTYPE,dimension(:,:,:),allocatable,private        :: idpdx0,idpdy0
#endif
   integer, private, parameter         :: BLUMBERG_MELLOR=1
   integer, private, parameter         :: BLUMBERG_MELLOR_LIN=2
   integer, private, parameter         :: Z_INTERPOL=3
   integer, private, parameter         :: SONG_WRIGHT=4
   integer, private, parameter         :: CHU_FAN=5
   integer, private, parameter         :: SHCHEPETKIN_MCWILLIAMS=6
   integer, private, parameter         :: STELLING_VANKESTER=7
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_internal_pressure - initialising internal pressure gradient
! \label{sec-init-internal pressure}
!
! !INTERFACE:
   subroutine init_internal_pressure(runtype,hotstart,nonhyd_method)
!
! !DESCRIPTION:
!
! Here, some necessary memory is allocated (in case of the compiler option
! {\tt STATIC}), and information is written to the log-file of
! the simulation.

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in) :: runtype,nonhyd_method
   logical, intent(in) :: hotstart
!
! !LOCAL VARIABLES
   integer         :: rc
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_internal_pressure() # ',Ncall
#endif

   LEVEL2 'init_internal_pressure()'

   select case (ip_method)
      case(BLUMBERG_MELLOR)
         LEVEL3 'Blumber-Mellor scheme'
      case(BLUMBERG_MELLOR_LIN)
         LEVEL3 'Blumber-Mellor linear scheme'
      case(Z_INTERPOL)
         LEVEL3 'Z-interpolation'
      case(SONG_WRIGHT)
         LEVEL3 'Song and Wright scheme'
      case(CHU_FAN)
         LEVEL3 'Chu and Fan scheme'
         call getm_error("init_internal_pressure()", &
             "Not working, use other internal pressure gradient scheme")
      case(SHCHEPETKIN_MCWILLIAMS)
         LEVEL3 'Shchepetkin and McWilliams (2003) scheme'
      case(STELLING_VANKESTER)
         LEVEL3 'Stelling and vanKester (1994) scheme'
      case default
         FATAL 'Not valid ip_method specified'
         stop 'init_internal_pressure()'
   end select

   if (ip_ramp .gt. 1) then
      LEVEL3 'ip_ramp=',ip_ramp
      ip_ramp_is_active = .true.
      if (hotstart) then
         LEVEL3 'WARNING: hotstart is .true. AND ip_ramp .gt. 1'
         LEVEL3 'WARNING: .. be sure you know what you are doing ..'
      end if
   end if

   if (runtype .ge. 3) then
      allocate(idpdx_hs(I3DFIELD),stat=rc) ! Internal pressure gradient - x
      if (rc /= 0) stop 'init_internal_pressure(): Error allocating memory (idpdx_hs)'
      idpdx_hs = _ZERO_
      idpdx_m3d => idpdx_hs

#ifndef SLICE_MODEL
      allocate(idpdy_hs(I3DFIELD),stat=rc) ! Internal pressure gradient - y
      if (rc /= 0) stop 'init_internal_pressure(): Error allocating memory (idpdy_hs)'
      idpdy_hs = _ZERO_
      idpdy_m3d => idpdy_hs
#endif
   end if

   if (nonhyd_method .eq. 1) then
      allocate(idpdx_nh(I3DFIELD),stat=rc) ! Internal pressure gradient - x
      if (rc /= 0) stop 'init_internal_pressure(): Error allocating memory (idpdx_nh)'
      idpdx_nh = _ZERO_
      idpdx_m3d => idpdx_nh

#ifndef SLICE_MODEL
      allocate(idpdy_nh(I3DFIELD),stat=rc) ! Internal pressure gradient - y
      if (rc /= 0) stop 'init_internal_pressure(): Error allocating memory (idpdy_nh)'
      idpdy_nh = _ZERO_
      idpdy_m3d => idpdy_nh
#endif
   end if


   if (runtype.ge.3 .and. nonhyd_method.eq.1) then

      calc_ipfull = .true.

      allocate(idpdx_full(I3DFIELD),stat=rc) ! Internal pressure gradient - x
      if (rc /= 0) stop 'init_internal_pressure(): Error allocating memory (idpdx_full)'
      idpdx_full = _ZERO_
      idpdx_m3d => idpdx_full

#ifndef SLICE_MODEL
      allocate(idpdy_full(I3DFIELD),stat=rc) ! Internal pressure gradient - y
      if (rc /= 0) stop 'init_internal_pressure(): Error allocating memory (idpdy_full)'
      idpdy_full = _ZERO_
      idpdy_m3d => idpdy_full
#endif

   end if


#ifndef STATIC
   allocate(zz(I3DFIELD),stat=rc)
   if (rc /= 0) stop 'init_internal_pressure: Error allocating memory (zz)'
   zz=_ZERO_
#endif

#ifdef SUBSTR_INI_PRESS
   allocate(idpdx0(I3DFIELD),stat=rc) ! Initial x - pressure gradient.
   if (rc /= 0) stop 'init_internal_pressure(): Error allocating memory (idpdx0)'
   idpdx0=_ZERO_

#ifndef SLICE_MODEL
   allocate(idpdy0(I3DFIELD),stat=rc) ! Initial y - pressure gradient.
   if (rc /= 0) stop 'init_internal_pressure(): Error allocating memory (idpdy0)'
   idpdy0=_ZERO_
#endif
#endif

   return
   end subroutine init_internal_pressure
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_internal_pressure - internal pressure gradient
! \label{sec-do-internal-pressure}
!
! !INTERFACE:
   subroutine do_internal_pressure(part)
!
! !DESCRIPTION:
!
! Here, the chosen internal pressure gradient method is selected
! and (in case that the compiler option {\tt SUBSTR\_INI\_PRESS} is
! set), the initial pressure is calculated and subtracted from the
! updated internal pressure gradient.
!
! If GETM is executed as slice model (compiler option {\tt SLICE\_MODEL}
! is set, the internal pressure gradient for $j=2$ is copied to
! $j=3$.
!
! !USES:
   use getm_timers, only: tic, toc, TIM_INTPRESS
   IMPLICIT NONE
!
! !INPUT VARIABLES:
   integer, intent(in) :: part
!
! !LOCAL VARIABLES:
   integer                :: i,j,k
   logical, save          :: first=.true.
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_internal_pressure() # ',Ncall
#endif
   call tic(TIM_INTPRESS)

   select case (part)
      case(1)
#ifndef NO_BAROCLINIC
         buoy  => buoy_eos
#endif
         idpdx => idpdx_hs
         idpdy => idpdy_hs
      case(2)
         buoy  => minus_bnh
         idpdx => idpdx_nh
         idpdy => idpdy_nh
   end select

! BJB-NOTE: Initialization of zz, idpdx and ipdy moved to the
!  individual ip_methods to allow speed-up based on method (by
!  local reduction of amount of initialzation) BJB 2009-09-24.

   select case (ip_method)
      case(BLUMBERG_MELLOR)
         call ip_blumberg_mellor()
      case(BLUMBERG_MELLOR_LIN)
         call ip_blumberg_mellor_lin()
      case(Z_INTERPOL)
         call ip_z_interpol()
      case(SONG_WRIGHT)
         call ip_song_wright()
      case(CHU_FAN)
         call ip_chu_fan()
      case(SHCHEPETKIN_MCWILLIAMS)
         call ip_shchepetkin_mcwilliams()
      case(STELLING_VANKESTER)
         call ip_stelling_vankester()
      case default
         FATAL 'Not valid ip_method specified'
         stop 'do_internal_pressure()'
   end select

   if (part .eq. 2) then
!     no nh pressure gradient across open bdys
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO
            if (au(i,j) .eq. 2) then
               idpdx_nh(i,j,:) = _ZERO_
            end if
#ifndef SLICE_MODEL
            if (av(i,j) .eq. 2) then
               idpdy_nh(i,j,:) = _ZERO_
            end if
#endif
         end do
      end do
   end if

#ifdef SUBSTR_INI_PRESS
   if (first) then
      first = .false.
      do k=0,kmax
         do j=jmin,jmax
            do i=imin,imax
               idpdx0(i,j,k) = idpdx(i,j,k)
               idpdx(i,j,k) = _ZERO_
#ifndef SLICE_MODEL
               idpdy0(i,j,k) = idpdy(i,j,k)
               idpdy(i,j,k) = _ZERO_
#endif
            end do
         end do
      end do
   else
      do k=0,kmax
         do j=jmin,jmax
            do i=imin,imax
               idpdx(i,j,k) = idpdx(i,j,k) - idpdx0(i,j,k)
#ifndef SLICE_MODEL
               idpdy(i,j,k) = idpdy(i,j,k) - idpdy0(i,j,k)
#endif
            end do
         end do
      end do
   end if
#endif

#ifdef SLICE_MODEL
   do i=imin,imax
      do k=kmin(i,2),kmax
         idpdx(i,3,k)=idpdx(i,2,k)
      end do
   end do
#endif

   if (calc_ipfull) then
      idpdx_full = idpdx_hs + idpdx_nh
#ifndef SLICE_MODEL
      idpdy_full = idpdy_hs + idpdy_nh
#endif
   end if

   call toc(TIM_INTPRESS)
#ifdef DEBUG
   write(debug,*) 'Leaving do_internal_pressure()'
   write(debug,*)
#endif
   return
   end subroutine do_internal_pressure
!EOC

!-----------------------------------------------------------------------

   end module internal_pressure

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
