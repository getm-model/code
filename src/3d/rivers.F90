!$Id: rivers.F90,v 1.6 2005-01-13 09:20:46 kbk Exp $
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
   integer, public                     :: river_method=0,nriver=0,rriver=0
   logical,public                      :: use_river_temp = .false.
   logical,public                      :: use_river_salt = .false.
   character(len=64), public           :: river_data="rivers.nc"
   character(len=64), public, allocatable  :: river_name(:)
   character(len=64), public, allocatable  :: real_river_name(:)
   integer, public, allocatable        :: ok(:)
   REALTYPE, public, allocatable       :: river_flow(:)
   REALTYPE, public, allocatable       :: river_salt(:)
   REALTYPE, public, allocatable       :: river_temp(:)
   REALTYPE, public                    :: river_factor= _ONE_
   REALTYPE, public,parameter          :: temp_missing=-9999.0
   REALTYPE, public,parameter          :: salt_missing=-9999.0
   integer,  public, allocatable       :: river_split(:)
!
! !PRIVATE DATA MEMBERS:
   integer                   :: river_format=2
   character(len=64)         :: river_info="riverinfo.dat"
   integer, allocatable      :: ir(:),jr(:)
   REALTYPE, allocatable     :: irr(:)
   REALTYPE, allocatable     :: macro_height(:)
   REALTYPE, allocatable     :: flow_fraction(:)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: rivers.F90,v $
!  Revision 1.6  2005-01-13 09:20:46  kbk
!  support for T and S specifications in rivers - Stips
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
   integer                   :: i,j,n,nn,ni,rc,m
   integer                   :: unit = 25 ! kbk
   logical                   :: outside
   REALTYPE                  :: area
   NAMELIST /rivers/ &
            river_method,river_info,river_format,river_data,river_factor, &
            use_river_salt,use_river_temp
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
         LEVEL2 'use_river_temp= ',use_river_temp
         LEVEL2 'use_river_salt= ',use_river_salt
         open(unit,file=river_info,action='read',status='old',err=90)
         read(unit,*) nriver
         allocate(ir(nriver),stat=rc) ! i index of rivers
         if (rc /= 0) stop 'rivers: Error allocating memory (ir)'
         allocate(jr(nriver),stat=rc) ! j index of rivers
         if (rc /= 0) stop 'rivers: Error allocating memory (jr)'
         allocate(ok(nriver),stat=rc) ! valid river spec.
         if (rc /= 0) stop 'rivers: Error allocating memory (ok)'
         allocate(river_name(nriver),stat=rc) ! NetCDF name of river.
         if (rc /= 0) stop 'rivers: Error allocating memory (river_name)'
         allocate(river_flow(nriver),stat=rc) ! river flux
         if (rc /= 0) stop 'rivers: Error allocating memory (river_flow)'
         allocate(macro_height(nriver),stat=rc) ! height over a macro tims-step
         if (rc /= 0) stop 'rivers: Error allocating memory (macro_height)'
         allocate(river_temp(nriver),stat=rc) ! temperature of river water
         if (rc /= 0) stop 'rivers: Error allocating memory (river_temp)'
         allocate(river_salt(nriver),stat=rc) ! salinity of river water
         if (rc /= 0) stop 'rivers: Error allocating memory (river_salt)'
         allocate(river_split(nriver),stat=rc) ! split factor for river water
         if (rc /= 0) stop 'rivers: Error allocating memory (river_split)'
         allocate(flow_fraction(nriver),stat=rc) ! areafactor of data for river 
         if (rc /= 0) stop 'rivers: Error allocating memory (flow_fraction)'
         allocate(irr(nriver),stat=rc) ! integrated river runoff
         if (rc /= 0) stop 'rivers: Error allocating memory (irr)'

         ok = 0
         rriver = 0 ! number of real existing rivers...
         flow_fraction = _ZERO_
         do n=1,nriver
            read(unit,*) ir(n),jr(n),river_name(n)
            river_name(n) = trim(river_name(n))
            LEVEL3 trim(river_name(n)),':',ir(n),jr(n)
            i = ir(n)-ioff
            j = jr(n)-joff
            river_temp(n) = temp_missing
            river_salt(n) = salt_missing
            river_flow(n) = _ZERO_
            irr(n) = _ZERO_
            macro_height(n) = _ZERO_
!           calculate the number of used rivers, they must be 
!           in sequence !
            rriver = rriver +1
            if ( n .gt. 1 ) then
               if (river_name(n) .eq. river_name(n-1)) then 
                  rriver = rriver-1
               end if
            end if
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
                  flow_fraction(n) = _ONE_/ARCD1
               end if
            else
!              LEVEL3 'Outside: river# ',n
            end if
         end do

!        calculate the number of used gridboxes, they must be 
!        in sequence !
         LEVEL3 'Number of unique rivers: ',rriver
         allocate(real_river_name(rriver),stat=rc) ! NetCDF name of river.
         if (rc /= 0) stop 'rivers: Error allocating memory (rivers)'
         river_split = 1    ! normal case
         do n=2,nriver
               if (river_name(n) .eq. river_name(n-1)) then 
                  river_split(n)=river_split(n-1)+1
            end if
         end do
         ni= nriver
         do n=1,nriver
            if (ni .ge. 1) then
               if ( river_split(ni) .gt. 1 ) then  
                  do m=1,river_split(ni)
                     river_split(ni-m+1) =  river_split(ni)
                  end do
               end if
               ni = ni - river_split(ni)
            end if
         end do
         LEVEL3 'split:',river_split
!        now river_split contains the number of gridboxes used 
!        for a single river
         nn = 1
         ni = 1
         do n=1,nriver
            if (ni .le. nriver) then
               real_river_name(nn) = river_name(ni)
               if ( river_split(ni) .gt. 1 ) then
                  area = _ZERO_
                  do m=1,river_split(ni) 
                     area = area +  flow_fraction(ni+m-1)
                  end do
                  do m=1,river_split(ni)
                     if ( area .gt. _ZERO_ ) then
                        flow_fraction(ni+m-1) = flow_fraction(ni+m-1)/area
                     else
                        flow_fraction(ni+m-1) = _ZERO_
                     end if
                  end do
               else
                  flow_fraction(ni) = _ONE_
               end if
               nn = nn + 1  
               ni = ni + river_split(ni)
            end if 
            if (ok(n) .eq. 0) then
               flow_fraction(n) = _ZERO_
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
   REALTYPE                  :: rvol,height
   REALTYPE                  :: svol,tvol,vol
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
               rvol = dtm * river_flow(n) * flow_fraction(n)
               irr(n) = irr(n) + rvol
               height = rvol * ARCD1
               z(i,j) = z(i,j) + height
#ifndef NO_BAROCLINIC
               macro_height(n)=macro_height(n)+height
!              on macrotime step adjust 3d fields
               if (do_3d) then
                  if (calc_salt) then
                     if ( river_salt(n) .ne. salt_missing ) then
                        S(i,j,1:kmax) = (S(i,j,1:kmax)*(H(i,j)+ssen(i,j))   &
                                      + river_salt(n)*macro_height(n))      &
                                      / (H(i,j)+ssen(i,j)+macro_height(n))
                     else
                        S(i,j,1:kmax) = S(i,j,1:kmax)*(H(i,j)+ssen(i,j))   &
                                      / (H(i,j)+ssen(i,j)+macro_height(n))
                     end if
                  end if
                  if (calc_temp .and. river_temp(n) .ne. temp_missing) then
                     T(i,j,1:kmax) = (T(i,j,1:kmax)*(H(i,j)+ssen(i,j))   &
                                      + river_temp(n)*macro_height(n))      &
                                      / (H(i,j)+ssen(i,j)+macro_height(n))
                  end if
!                 Changes of total and layer height due to river inflow:
                  hn(i,j,1:kmax) = hn(i,j,1:kmax)/(H(i,j)+ssen(i,j)) &
                                  *(H(i,j)+ssen(i,j)+macro_height(n))
                  ssen(i,j) = ssen(i,j)+macro_height(n)
                  macro_height(n) = _ZERO_
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
