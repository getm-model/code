!$Id: internal_pressure.F90,v 1.12 2006-01-28 20:07:54 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: internal_pressure()
!
! !INTERFACE:
   module internal_pressure
!
! !DESCRIPTION:
!  ip_method=1    ! old GETM version of Blumberg-Mellor, very stable
!  ip_method=2    ! new GETM version of 1, zero for linear profiles 
!  ip_method=3    ! z-coordinates, linear interpolation
!  ip_method=4    ! Song and Wright 1998 version
!  ip_method=5    ! CHU and FAN 2003 version
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
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: internal_pressure.F90,v $
!  Revision 1.12  2006-01-28 20:07:54  hb
!  Extensions to compiler option SLICE_MODEL for better representation of zero gradients in y-direction
!
!  Revision 1.11  2005-11-17 13:50:22  kbk
!  fixes to compile with gfortran
!
!  Revision 1.10  2005/04/25 07:55:50  kbk
!  use more general frame for error handling - Umlauf
!
!  Revision 1.9  2004/06/18 13:43:50  kbk
!  fixed first argument to getm_error()
!
!  Revision 1.8  2004/06/18 13:00:30  hb
!  Chu and Fan disactivated because of inconsistencies
!
!  Revision 1.7  2004/05/03 09:31:47  lars
!  bug fix: defined i,j,k and first
!
!  Revision 1.6  2004/04/21 09:21:25  lars
!  changed to default ip_method=1
!
!  Revision 1.5  2004/04/06 12:42:50  kbk
!  internal pressure calculations now uses wrapper
!
!  Revision 1.3  2003/04/23 12:16:34  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/01 15:50:13  gotm
!  removed dead print statement
!
!  Revision 1.1.1.1  2002/05/02 14:00:59  gotm
!  recovering after CVS crash
!
!  Revision 1.15  2001/10/22 08:23:53  bbh
!  Removed reference to kplus and kminus when not -DPRESS_GRAD_Z
!
!  Revision 1.14  2001/10/22 07:47:59  bbh
!  Needed to run dos2unix
!
!  Revision 1.13  2001/10/22 07:45:26  bbh
!  Fixed a serious bug - now gives the same as old version
!
!  Revision 1.10  2001/09/04 07:29:18  bbh
!  Internal pressure based on interpolation to z-coordinates
!
!  Revision 1.9  2001/09/03 13:03:37  bbh
!  Initial pressure gradient can now be subtracted
!
!  Revision 1.8  2001/08/31 15:40:37  bbh
!  initial pressure can be subtracted now
!
!  Revision 1.7  2001/08/01 08:31:22  bbh
!  CURVILINEAR now implemented
!
!  Revision 1.6  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.5  2001/05/25 19:09:48  bbh
!  Removed dead code
!
!  Revision 1.4  2001/05/20 09:19:09  bbh
!  Use who specified twice
!
!  Revision 1.3  2001/05/20 07:51:40  bbh
!  Internal pressure included
!
!  Revision 1.2  2001/05/11 13:47:00  bbh
!  Added actual code
!
!  Revision 1.1  2001/05/10 11:30:16  bbh
!  Added further support for baroclinicity
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_internal_pressure
!
! !INTERFACE:
   subroutine init_internal_pressure()
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Reads the namelist and makes calls to the init functions of the
!  various model components.
!
! !REVISION HISTORY:
!  See the log for the module
!
!  !LOCAL VARIABLES
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
! !IROUTINE:  do_internal_pressure()
!
! !INTERFACE:
   subroutine do_internal_pressure()
!
! !DESCRIPTION:
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
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer                :: i,j,k
   logical, save          :: first=.true.
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

#ifdef SLICE_MODEL
   do i=iimin,iimax
      do k=kmin(i,2),kmax
         idpdx(i,3,k)=idpdx(i,2,k)
      end do   
   end do
#endif

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
