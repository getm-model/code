!$Id: rivers.F90,v 1.5 2004-11-12 07:29:36 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  rivers
!
! !INTERFACE:
   module rivers
!
! !DESCRIPTION:
!  This module includes support for river input. Rivers are treated the same
!  way as meteorology, i.e. as external to the hydrodynamic model itself.
!  The model follows the same scheme as all other modules, i.e. init\_rivers
!  sets up necessary information, do\_rivers update the relevant variables.
!  do\_river should be called before any 2d and 3d routines as it only
!  updates the sea surface elevation (2d) and sea surface elevation +
!  optional salinity and temperature (3d). Any update of other parameters,
!  (D, DU, DV...) are done through the relevant subroutines.
!  At present the momentum of the river water is not include, the model
!  however has a direct response to the river water because of the
!  pressure gradient introduced.
!
! !USES:
   use domain, only: iimin,jjmin,iimax,jjmax,ioff,joff
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: H,az,kmax,arcd1
#else
   use domain, only: H,az,kmax,ard1
#endif
   use m2d, only: dtm
   use variables_2d, only: z
#ifndef NO_BAROCLINIC
   use m3d, only: calc_salt,calc_temp
   use variables_3d, only: hn,ssen,T,S
#endif
   IMPLICIT NONE
!
   private
!
! !PUBLIC DATA MEMBERS:
   public init_rivers, do_rivers, clean_rivers
   integer, public                     :: river_method=0,nriver=0
   character(len=64), public           :: river_data="rivers.nc"
   character(len=64), public, allocatable  :: river_name(:)
   integer, public, allocatable        :: ok(:)
   REALTYPE, public, allocatable       :: river_flow(:),tr(:)
   REALTYPE, public                    :: river_factor= _ONE_
   REALTYPE, allocatable               :: macro_height(:)
!
! !PRIVATE DATA MEMBERS:
   integer                   :: river_format=2
   character(len=64)         :: river_info="riverinfo.dat"
   integer, allocatable      :: ir(:),jr(:)
   REALTYPE, allocatable     :: irr(:)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: rivers.F90,v $
!  Revision 1.5  2004-11-12 07:29:36  kbk
!  initialising ok to 0 - Buchmann
!
!  Revision 1.4  2003/10/14 10:05:54  kbk
!  checks if indices are in subdomain + cleaning
!
!  Revision 1.3  2003/04/23 12:16:34  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 13:36:38  kbk
!  parallel support, cleaned code + NO_3D, NO_BAROCLINIC
!
!  Revision 1.1.1.1  2002/05/02 14:01:00  gotm
!  recovering after CVS crash
!
!  Revision 1.3  2001/10/07 14:50:22  bbh
!  Reading river data implemented - NetCFD
!
!  Revision 1.2  2001/09/19 14:21:13  bbh
!  Cleaning
!
!  Revision 1.1  2001/09/18 17:48:32  bbh
!  Added algoritm for rivers - getting river data still missing
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_rivers
!
! !INTERFACE:
   subroutine init_rivers()
!
! !DESCRIPTION:
!  Description still missing
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
   integer                   :: i,j,n,rc
   integer                   :: unit = 25 ! kbk
   logical                   :: outside
   NAMELIST /rivers/ &
            river_method,river_info,river_format,river_data,river_factor
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_rivers() # ',Ncall
#endif

   LEVEL1 'init_rivers()'
   read(NAMLST,rivers)

   select case (river_method)
      case (0)
         LEVEL3 'River runoff not included.'
      case (1,2)
         LEVEL2 'river_method= ',river_method
         LEVEL2 'river_data=   ',trim(river_data)
         LEVEL2 'river_format= ',river_format
         open(unit,file=river_info,action='read',status='old',err=90)
         read(unit,*) nriver
         allocate(ir(nriver),stat=rc) ! i index of rivers
         allocate(jr(nriver),stat=rc) ! j index of rivers
         allocate(river_name(nriver),stat=rc) ! NetCDF name of river.
         allocate(ok(nriver),stat=rc) ! valid river spec.
         allocate(river_flow(nriver),stat=rc) ! river flux
         allocate(macro_height(nriver),stat=rc) ! height over a macro tims-step
         allocate(tr(nriver),stat=rc) ! temperature of river water
         allocate(irr(nriver),stat=rc) ! integrated river runoff
         ok = 0
         do n=1,nriver
            read(unit,*) ir(n),jr(n),river_name(n)
            LEVEL3 trim(river_name(n)),':',ir(n),jr(n)
            i = ir(n)-ioff
            j = jr(n)-joff
            tr(n) = _ZERO_
            irr(n) = _ZERO_
            macro_height(n) = _ZERO_
            outside= &
                    i .lt. iimin .or. i .gt. iimax .or.  &
                    j .lt. jjmin .or. j .gt. jjmax
            if( .not. outside) then
               if(az(i,j) .eq. 0) then
                  LEVEL3 'Warning:  river# ',n,' at (',i,j,') is on land'
                  LEVEL3 '          setting ok to 0'
                  ok(n) = 0
               else
                  ok(n) = 1
               end if
            else
!              LEVEL3 'Outside: river# ',n
            end if
         end do
      case default
         FATAL 'A non valid river_method has been selected'
         stop 'init_rivers'
   end select
   return

90 LEVEL2 'could not open ',trim(river_info),' for reading info on rivers'
   stop 'init_rivers()'

#ifdef DEBUG
   write(debug,*) 'Leaving init_rivers()'
   write(debug,*)
#endif
   return
   end subroutine init_rivers
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_rivers()
!
! !INTERFACE:
   subroutine do_rivers(do_3d)
!
! !DESCRIPTION:
!  Description still missing
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical, intent(in)                 :: do_3d
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,n
   REALTYPE                  :: vol
   REALTYPE                  :: rvol,height
   REALTYPE                  :: svol,tvol
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_rivers() # ',Ncall
#endif

   select case (river_method)
      case(0)
      case(1,2)
         do n=1,nriver
            if(ok(n) .gt. 0) then
               i = ir(n)-ioff; j = jr(n)-joff
               rvol = dtm*river_flow(n)
               irr(n) = irr(n) + rvol
               height = rvol*ARCD1
               macro_height(n)=macro_height(n)+height
               z(i,j) = z(i,j) + height
#ifndef NO_BAROCLINIC
               if (do_3d) then
                  if (calc_salt) then
                     S(i,j,1:kmax) = S(i,j,1:kmax)*(H(i,j)+ssen(i,j))   &
                                    /(H(i,j)+ssen(i,j)+macro_height(n))
                  end if

! Changes of total and layer height due to river inflow:
                  hn(i,j,1:kmax) = hn(i,j,1:kmax)/(H(i,j)+ssen(i,j)) &
                                  *(H(i,j)+ssen(i,j)+macro_height(n))
                  ssen(i,j) = ssen(i,j)+macro_height(n)
                  macro_height(n) = _ZERO_
#if 0
                  if (calc_temp .and. tr(n) .gt. _ZERO_) then
                     tvol = _ZERO_
                     do k=1,kmax
                        tvol = tvol+hn(i,j,k)*T(i,j,k)
                     end do
                     tvol = tvol*ARCD1
                     tvol = tvol + tr(n)*rvol
                     T(i,j,1:kmax) = tvol/vol
                  end if
#endif
               end if
#endif
            end if
         end do
      case default
         FATAL 'Not valid rivers_method specified'
         stop 'init_rivers'
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving do_rivers()'
   write(debug,*)
#endif
   return
   end subroutine do_rivers
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  clean_rivers()
!
! !INTERFACE:
   subroutine clean_rivers()
!
! !DESCRIPTION:
!  Description still missing
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
   integer                   :: i,j,n
   REALTYPE                  :: tot=_ZERO_
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'clean_rivers() # ',Ncall
#endif

   select case (river_method)
      case(0)
      case(1,2)
         do n=1,nriver
            if(ok(n) .gt. 0) then
               i = ir(n); j = jr(n)
               LEVEL2 trim(river_name(n)),':  ' ,irr(n)/1.e6, '10^6 m3'
               tot = tot+irr(n)
            end if
         end do
      case default
         FATAL 'Not valid rivers_method specified'
         stop 'init_rivers'
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving do_rivers()'
   write(debug,*)
#endif
   return
   end subroutine clean_rivers
!EOC

!-----------------------------------------------------------------------

   end module rivers

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Karsten Bolding and Hans Burchard               !
!-----------------------------------------------------------------------
