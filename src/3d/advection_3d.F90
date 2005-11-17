!$id: advection_3d.F90,v 1.18 2001/09/19 13:53:08 bbh Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  3D advection
!
! !INTERFACE:
   module advection_3d
!
! !DESCRIPTION:
!  This module do advection of scalars.  The module follows the same
!  convention as the other modules in 'getm'. The module is initialised
!  by calling 'init\_advection\_3d()'. In the time-loop 'do\_advection\_3d' is
!  called. 'do\_advection\_3d' is a wrapper routine which - dependent on the
!  actual advection scheme chosen - makes a call to the appropriate subroutine.
!  New advection schemes are easily implemented - at least from a program
!  point of view - since only this module needs to be changed.
!  Additional work arrays can easily be added following the stencil given
!  below. To add a new advection scheme 3 things must be done: 1) define
!  a unique constant to identify the scheme (see e.g. UPSTREAM and TVD)
!  2) adopt the 'select case' in 'do\_advection\_3d()' and 3) write the actual
!  subroutine.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
   use domain, only: iimin,iimax,jjmin,jjmax,kmax
   use halo_zones, only: update_3d_halo,wait_halo,D_TAG
   IMPLICIT NONE
!
   private
!
! !PUBLIC DATA MEMBERS:
   public init_advection_3d, do_advection_3d
#ifdef STATIC
   REALTYPE, public                    :: cu(I3DFIELD)
   REALTYPE, public                    :: hi(I3DFIELD)
   REALTYPE, public                    :: hio(I3DFIELD)
#else
   REALTYPE, public, dimension(:,:,:), allocatable       :: hi,hio,cu
#endif
   integer, public, parameter          :: UPSTREAM=1,UPSTREAM_SPLIT=2,P2=3
   integer, public, parameter          :: Superbee=4,MUSCL=5,P2_PDM=6,FCT=7
   REALTYPE, public, parameter         :: one6th=1./6.
   REALTYPE, public, parameter         :: ONE=_ONE_,TWO=2.*_ONE_
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: advection_3d.F90,v $
!  Revision 1.9  2005-11-17 13:50:22  kbk
!  fixes to compile with gfortran
!
!  Revision 1.8  2005/10/06 09:54:01  hb
!  added support for vertical slice model - via -DSLICE_MODEL
!
!  Revision 1.7  2005/05/25 10:32:13  kbk
!  merged from stabe branch v1_2_1
!
!  Revision 1.6.2.1  2005/05/25 08:39:14  kbk
!  update HALO's after each fractional step
!
!  Revision 1.6  2004/01/06 15:04:00  kbk
!  FCT advection + split of advection_3d.F90 + extra adv. input checks
!
!  Revision 1.5  2003/12/16 16:50:40  kbk
!  added support for Intel/IFORT compiler - expanded TABS, same types in subroutine calls
!
!  Revision 1.4  2003/09/03 05:38:45  kbk
!  need to call update_3d_halo() for each directional split
!
!  Revision 1.3  2003/04/23 12:16:34  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 16:30:53  kbk
!  parallel support
!
!  Revision 1.1.1.1  2002/05/02 14:00:58  gotm
!  recovering after CVS crash
!
!  Revision 1.18  2001/09/19 13:53:08  bbh
!  Typo
!
!  Revision 1.17  2001/09/19 13:46:52  bbh
!  Nasty memory leak in upstream_adv() - only if -DFORTRAN90
!
!  Revision 1.16  2001/09/03 20:13:08  bbh
!  Corrected order of select (strang)
!
!  Revision 1.15  2001/09/03 20:04:21  bbh
!  Allow individual advection settings for momentum, salinity and temperature
!
!  Revision 1.14  2001/09/03 12:56:58  bbh
!  Advection can now be split into different schemes for each direction
!
!  Revision 1.13  2001/08/29 14:22:59  bbh
!  Dimensions for masks are needed
!
!  Revision 1.12  2001/08/29 12:44:39  bbh
!  masks are E2DFIELD
!
!  Revision 1.11  2001/08/27 11:51:45  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.10  2001/08/01 08:31:22  bbh
!  CURVILINEAR now implemented
!
!  Revision 1.9  2001/07/26 12:49:33  bbh
!  Added a number of new advection schemes - UPSTREAM_SPLIT, P2,
!  Superbee, MUSCL and P2_PDM.
!
!  Revision 1.8  2001/05/22 08:24:21  bbh
!  Added adv in call to advection routines
!
!  Revision 1.7  2001/05/18 09:45:16  bbh
!  Removed some LEVEL2 statements
!
!  Revision 1.6  2001/05/18 09:37:07  bbh
!  IMPLICITE NONE was commented out  at one place
!
!  Revision 1.5  2001/05/18 07:59:46  bbh
!  UPSTREAM advection is ready
!
!  Revision 1.4  2001/05/11 13:47:53  bbh
!  First do UPSTREAM
!
!  Revision 1.3  2001/05/04 13:15:05  bbh
!  Public members was commented out
!
!  Revision 1.2  2001/05/04 09:53:06  bbh
!  Added stubs for upstream_adv() and tdv_adv + description
!
!  Revision 1.1  2001/05/03 20:20:33  bbh
!  Stubs for baroclinicity
!
! !LOCAL VARIABLES:
   integer :: advection_method
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_advection_3d
!
! !INTERFACE:
   subroutine init_advection_3d(method)
!
! !DESCRIPTION:
!  Reads the namelist and makes calls to the init functions of the
!  various model components.
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: method
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer                   :: rc
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_advection_3d() # ',Ncall
#endif

   LEVEL1 'init_advection_3d()'

#ifndef STATIC
   allocate(cu(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection_3d: Error allocating memory (cu)'

   allocate(hi(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection_3d: Error allocating memory (hi)'

   allocate(hio(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection_3d: Error allocating memory (hio)'
#endif

   cu  = _ZERO_ ; hi  = _ZERO_ ; hio = _ZERO_ 

#ifdef DEBUG
   write(debug,*) 'Leaving init_advection_3d()'
   write(debug,*)
#endif
   return
   end subroutine init_advection_3d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_advection_3d()
!
! !INTERFACE:
   subroutine do_advection_3d(dt,f,uu,vv,ww,hun,hvn,ho,hn,      &
                             delxu,delxv,delyu,delyv,area_inv,  &
                             az,au,av,hor_adv,ver_adv,adv_split,AH)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)      :: uu(I3DFIELD)
   REALTYPE, intent(in)      :: vv(I3DFIELD)
   REALTYPE, intent(in)      :: ww(I3DFIELD)
   REALTYPE, intent(in)      :: ho(I3DFIELD)
   REALTYPE, intent(in)      :: hn(I3DFIELD)
   REALTYPE, intent(in)      :: hun(I3DFIELD)
   REALTYPE, intent(in)      :: hvn(I3DFIELD)
   REALTYPE, intent(in)      :: delxu(I2DFIELD)
   REALTYPE, intent(in)      :: delxv(I2DFIELD)
   REALTYPE, intent(in)      :: delyu(I2DFIELD)
   REALTYPE, intent(in)      :: delyv(I2DFIELD)
   REALTYPE, intent(in)      :: area_inv(I2DFIELD),dt,AH
   integer, intent(in)       :: az(E2DFIELD)
   integer, intent(in)       :: au(E2DFIELD)
   integer, intent(in)       :: av(E2DFIELD)
   integer, intent(in)       :: hor_adv,ver_adv,adv_split
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout)   :: f(I3DFIELD)
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   REALTYPE, parameter       :: a1=0.5*ONE,a2=ONE
   integer         :: k
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_advection_3d() # ',Ncall
#endif

   select case (hor_adv)
      case (UPSTREAM)
         call upstream_adv(dt,f,uu,vv,ww,ho,hn, &
                           delxv,delyu,delxu,delyv,area_inv,az,AH)
      case ((UPSTREAM_SPLIT),(P2),(Superbee),(MUSCL),(P2_PDM),(FCT))
         hi=ho
         select case (adv_split)
            case (0)
               call u_split_adv(dt,f,uu,hun,delxu,delyu,area_inv,au,a2,&
                                hor_adv,az,AH)
               call update_3d_halo(f,f,az,& 
                                   iimin,jjmin,iimax,jjmax,kmax,D_TAG)
               call wait_halo(D_TAG)

#ifndef SLICE_MODEL
               call v_split_adv(dt,f,vv,hvn,delxv,delyv,area_inv,av,a2,&
                                hor_adv,az,AH)
               call update_3d_halo(f,f,az,& 
                                   iimin,jjmin,iimax,jjmax,kmax,D_TAG)
               call wait_halo(D_TAG)
#endif

               if (kmax.gt.1) then
#ifdef ITERATE_VERT_ADV
                  call w_split_it_adv(dt,f,ww,az,a2,ver_adv)
#else
                  call w_split_adv(dt,f,ww,az,a2,ver_adv)
#endif
               end if
            case (1)
               call u_split_adv(dt,f,uu,hun,delxu,delyu,area_inv,au,a1,&
                                hor_adv,az,AH)
               call update_3d_halo(f,f,az, &
                                   iimin,jjmin,iimax,jjmax,kmax,D_TAG)
               call wait_halo(D_TAG)

#ifndef SLICE_MODEL
               call v_split_adv(dt,f,vv,hvn,delxv,delyv,area_inv,av,a1,&
                                hor_adv,az,AH)
               call update_3d_halo(f,f,az, &
                                   iimin,jjmin,iimax,jjmax,kmax,D_TAG)
               call wait_halo(D_TAG)
#endif

               if (kmax.gt.1) then
#ifdef ITERATE_VERT_ADV
                  call w_split_it_adv(dt,f,ww,az,a2,ver_adv)
#else
                  call w_split_adv(dt,f,ww,az,a2,ver_adv)
#endif
                  call update_3d_halo(f,f,az, &
                                      iimin,jjmin,iimax,jjmax,kmax,D_TAG)
                  call wait_halo(D_TAG)

               end if
#ifndef SLICE_MODEL
               call v_split_adv(dt,f,vv,hvn,delxv,delyv,area_inv,av,a1,&
                                hor_adv,az,AH)
               call update_3d_halo(f,f,az, &
                                   iimin,jjmin,iimax,jjmax,kmax,D_TAG)
               call wait_halo(D_TAG)
#endif

               call u_split_adv(dt,f,uu,hun,delxu,delyu,area_inv,au,a1,&
                                hor_adv,az,AH)
               call update_3d_halo(f,f,az, &
                                   iimin,jjmin,iimax,jjmax,kmax,D_TAG)
               call wait_halo(D_TAG)

            case (2)
               select case (hor_adv)
                  case (UPSTREAM_SPLIT)
                     call upstream_2dh_adv(dt,f,uu,vv,ho,hn,hun,hvn, &
                               delxv,delyu,delxu,delyv,area_inv,az,AH)
                  case (FCT)
                     call fct_2dh_adv(dt,f,uu,vv,ho,hn,hun,hvn, &
                               delxv,delyu,delxu,delyv,area_inv,az,AH)
                  case default
                     FATAL 'For adv_split=2, hor_adv must be 2 (upstream) or 7 (fct)'
               end select
               if (kmax .gt. 1) then
#ifdef ITERATE_VERT_ADV
                  call w_split_it_adv(dt,f,ww,az,a2,ver_adv)
#else
                  call w_split_adv(dt,f,ww,az,a2,ver_adv)
#endif
               end if
            case default
               FATAL 'Not valid adv_split parameter'
         end select
      case default
         FATAL 'This is not so good - do_advection_3d()'
         stop
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving do_advection_3d()'
   write(debug,*)
#endif
   return
   end subroutine do_advection_3d
!EOC

!-----------------------------------------------------------------------

   end module advection_3d

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
