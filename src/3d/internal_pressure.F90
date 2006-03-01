!$Id: internal_pressure.F90,v 1.15 2006-03-01 14:45:12 hb Exp $
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
   use domain, only: iimin,iimax,jjmin,jjmax,kmax,az,au,av,H,HU,HV
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxu,dyv
#else
   use domain, only: dx,dy
#endif
   use variables_3d, only: kmin,hn,hun,hvn,idpdx,idpdy,rho
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   public init_internal_pressure, do_internal_pressure
   integer, public           :: ip_method=1
#ifdef STATIC
   REALTYPE                  :: zz(I3DFIELD)
#ifdef SUBSTR_INI_PRESS
   REALTYPE                  :: idpdx0(I3DFIELD),idpdy0(I3DFIELD)
#endif
#else
   REALTYPE, allocatable     :: zz(:,:,:)
#ifdef SUBSTR_INI_PRESS
   REALTYPE, allocatable     :: idpdx0(:,:,:),idpdy0(:,:,:)
#endif
#endif
!
! !PRIVATE DATA MEMBERS:
   integer, private, parameter         :: BLUMBERG_MELLOR=1
   integer, private, parameter         :: BLUMBERG_MELLOR_LIN=2
   integer, private, parameter         :: Z_INTERPOL=3
   integer, private, parameter         :: SONG_WRIGHT=4
   integer, private, parameter         :: CHU_FAN=5
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
   subroutine init_internal_pressure()
   IMPLICIT NONE
!
! !DESCRIPTION:
!  
! Here, some necessary memory is allocated (in case of the compiler option
! {\tt STATIC}), and information is written to the log-file of
! the simulation.
!
! !LOCAL VARIABLES
   integer         :: rc
!EOP
!-------------------------------------------------------------------------
!BOC

#ifndef STATIC
   allocate(zz(I3DFIELD),stat=rc)
   if (rc /= 0) stop 'init_internal_pressure: Error allocating memory (zz)'
#ifdef SUBSTR_INI_PRESS
   allocate(idpdx0(I3DFIELD),stat=rc) ! Initial x - pressure gradient.
   if (rc /= 0) stop &
        'init_internal_pressure(): Error allocating memory (idpdx0)'
   allocate(idpdy0(I3DFIELD),stat=rc) ! Initial y - pressure gradient.
   if (rc /= 0) stop &
         'init_internal_pressure(): Error allocating memory (idpdy0)'
   idpdx0=_ZERO_ ;  idpdy0=_ZERO_
#endif
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
      case default
         FATAL 'Not valid ip_method specified'
         stop 'init_internal_pressure()'
   end select

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
   subroutine do_internal_pressure()
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
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
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

   zz = _ZERO_
   idpdx = _ZERO_
   idpdy = _ZERO_

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
      case default
         FATAL 'Not valid ip_method specified'
         stop 'do_internal_pressure()'
   end select

#ifdef SUBSTR_INI_PRESS
   if (first) then
      first = .false.
      do k=0,kmax
         do j=jjmin,jjmax
            do i=iimin,iimax
               idpdx0(i,j,k) = idpdx(i,j,k)
               idpdx(i,j,k) = _ZERO_
               idpdy0(i,j,k) = idpdy(i,j,k)
               idpdy(i,j,k) = _ZERO_
            end do
         end do
      end do
   else
      do k=0,kmax
         do j=jjmin,jjmax
            do i=iimin,iimax
               idpdx(i,j,k) = idpdx(i,j,k) - idpdx0(i,j,k)
               idpdy(i,j,k) = idpdy(i,j,k) - idpdy0(i,j,k)
            end do
         end do
      end do
   end if
#endif

#ifdef SLICE_MODEL
   do i=iimin,iimax
      do k=kmin(i,2),kmax
         idpdx(i,3,k)=idpdx(i,2,k)
      end do   
   end do
#endif


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
