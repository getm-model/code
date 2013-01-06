#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  rivers \label{sec-rivers}
!
! !INTERFACE:
   module rivers
!
! !DESCRIPTION:
!
!  This module includes support for river input. Rivers are treated the same
!  way as meteorology, i.e.\ as external module to the hydrodynamic model
!  itself.
!  The module follows the same scheme as all other modules, i.e.\
!  {\tt init\_rivers}
!  sets up necessary information, and {\tt do\_rivers} updates
!  the relevant variables.
!  {\tt do\_river} is called in {\tt getm/integration.F90}
!  between the {\tt 2d} and {\tt 3d} routines as it only
!  updates the sea surface elevation (in {\tt 2d}) and sea surface elevation,
!  and
!  optionally salinity and temperature (in {\tt 3d}).
!  At present the momentum of the river water is not include, the model
!  however has a direct response to the river water because of the
!  pressure gradient introduced.
!
! !USES:
   use domain, only: imin,jmin,imax,jmax,ioff,joff
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: H,az,kmax,arcd1
#else
   use domain, only: H,az,kmax,ard1
#endif
   use time, only: write_time_string,timestr
   use variables_2d, only: dtm,z
#ifndef NO_BAROCLINIC
   use m3d, only: calc_salt,calc_temp
   use variables_3d, only: hn,ssen,T,S
#endif
#ifdef GETM_BIO
   use bio, only: bio_calc
   use bio_var, only: numc
   use variables_3d, only: cc3d
#endif
#ifdef _FABM_
   use gotm_fabm, only: model
   use getm_fabm, only: fabm_pel
#endif
   IMPLICIT NONE
!
   private
!
! !PUBLIC DATA MEMBERS:
   public init_rivers, do_rivers, clean_rivers
#ifdef GETM_BIO
   public init_rivers_bio
#endif
#ifdef _FABM_
   public init_rivers_fabm
#endif
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
   integer, public                     :: river_ramp= -1
   REALTYPE, public                    :: river_factor= _ONE_
   REALTYPE, public,parameter          :: temp_missing=-9999.0
   REALTYPE, public,parameter          :: salt_missing=-9999.0
   integer,  public, allocatable       :: river_split(:)
#ifdef GETM_BIO
   REALTYPE, public, allocatable       :: river_bio(:,:)
   REALTYPE, public, parameter         :: bio_missing=-9999.0
#endif
#ifdef _FABM_
   REALTYPE, public, allocatable       :: river_fabm(:,:)
#endif
!
! !PRIVATE DATA MEMBERS:
   integer                   :: river_format=2
   character(len=64)         :: river_info="riverinfo.dat"
   integer, allocatable      :: ir(:),jr(:)
   REALTYPE, allocatable     :: irr(:)
   REALTYPE, allocatable     :: macro_height(:)
   REALTYPE, allocatable     :: flow_fraction(:),flow_fraction_rel(:)
   REALTYPE                  :: ramp=_ONE_
   logical                   :: ramp_is_active=.false.
   logical                   :: river_outflow_properties_follow_source_cell=.true.
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
! !IROUTINE: init_rivers
!
! !INTERFACE:
   subroutine init_rivers(hotstart)
!
! !DESCRIPTION:
!
! First of all, the namelist {\tt rivers} is read from getm.F90 and
! a number of vectors with the length of {\tt nriver} (number of
! rivers) is allocated. Then, by looping over all rivers, the
! ascii file {\tt river\_info} is read, and checked for consistency.
! The number of used rivers {\tt rriver} is calculated and it is checked
! whether they are on land (which gives a warning) or not. When a river name
! occurs more than once in {\tt river\_info}, it means that its runoff
! is split among several grid boxed (for wide river mouths).
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical,intent(in)        :: hotstart
!
! !LOCAL VARIABLES:
   integer                   :: i,j,n,nn,ni,rc,m,iriver,jriver,numcells
   integer                   :: unit = 25 ! kbk
   logical                   :: outside,outsidehalo
   REALTYPE                  :: area, total_weight
   NAMELIST /rivers/ &
            river_method,river_info,river_format,river_data,river_ramp, &
            river_factor,use_river_salt,use_river_temp,river_outflow_properties_follow_source_cell
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
         if (river_ramp .gt. 1) then
            LEVEL2 'river_ramp=',river_ramp
            ramp_is_active = .true.
            if (hotstart) then
               LEVEL3 'WARNING: hotstart is .true. AND river_ramp .gt. 1'
               LEVEL3 'WARNING: .. be sure you know what you are doing ..'
            end if
         end if
         LEVEL2 'river_factor= ',river_factor
         LEVEL2 'use_river_temp= ',use_river_temp
         LEVEL2 'use_river_salt= ',use_river_salt
         LEVEL2 'river_outflow_properties_follow_source_cell=',river_outflow_properties_follow_source_cell
         open(unit,file=river_info,action='read',status='old',err=90)
         read(unit,*) nriver
         allocate(ir(nriver),stat=rc) ! i index of rivers
         if (rc /= 0) stop 'rivers: Error allocating memory (ir)'
         allocate(jr(nriver),stat=rc) ! j index of rivers
         if (rc /= 0) stop 'rivers: Error allocating memory (jr)'
         allocate(ok(nriver),stat=rc) ! valid river spec. 1=yes, 0=other domain, -1=other domain, but need read.
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
         allocate(flow_fraction_rel(nriver),stat=rc) ! Weight factor of data for river
         if (rc /= 0) stop 'rivers: Error allocating memory (flow_fraction_rel)'
         allocate(flow_fraction(nriver),stat=rc) ! Weight factor of data for river - scaled to unity sum for river
         if (rc /= 0) stop 'rivers: Error allocating memory (flow_fraction)'
         allocate(irr(nriver),stat=rc) ! integrated river runoff
         if (rc /= 0) stop 'rivers: Error allocating memory (irr)'

! ok(:) flags the location of each river relative to present subdomain:
!   1: River is present domain
!   2: River is in HALO of present domain
!  -1: River is multi-cell, where another part is inside present domain
!   0: River is outside present domain (no association)
         ok = 0
         rriver = 0 ! number of real existing rivers...
         flow_fraction_rel = _ZERO_
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
! Other weighting schemes could be implemented here. But we can only use
! information, which is available for cells also outside the present subdomain.
!           flow_fraction(n) = _ONE_/ARCD1 ! This does not work.
            flow_fraction_rel(n) = _ONE_
            outside= &
                    i .lt. imin .or. i .gt. imax .or.  &
                    j .lt. jmin .or. j .gt. jmax
            outsidehalo= &
                    (i .lt. imin-HALO) .or. (i .gt. imax+HALO) .or.  &
                    (j .lt. jmin-HALO) .or. (j .gt. jmax+HALO)
            if( .not. outsidehalo) then
               if(az(i,j) .eq. 0) then
                  LEVEL3 'Warning:  river# ',n,' at (',i,j,') is on land'
                  LEVEL3 '          setting ok to 0'
                  ok(n) = 0
               else if (.not. outside) then
                  LEVEL3 'Inside ',trim(river_name(n)),': river# ',n
                  ok(n) = 1
               else
                  LEVEL3 'Inhalo ',trim(river_name(n)),': river# ',n
                  ok(n) = 2
               end if
            else
              LEVEL3 'Outside: river# ',n
            end if
         end do

!  Calculate the number of used gridboxes.
!  This particular section is prepared for multi-cell rivers not-in-sequence.
         LEVEL3 'Number of unique rivers: ',rriver
         allocate(real_river_name(rriver),stat=rc) ! NetCDF name of river.
         if (rc /= 0) stop 'rivers: Error allocating memory (real_river_name)'
         river_split(:) = 1    ! normal case
!  This is a brute-force implementation, N**2, but N is presumed low.
         do iriver=1,nriver
            numcells=0
            total_weight=_ZERO_
            do jriver=1,nriver
               if (river_name(iriver) .eq. river_name(jriver)) then
                  numcells     = numcells+1
                  total_weight = total_weight+flow_fraction_rel(jriver)
                  if (ok(iriver).gt.0 .and. ok(jriver).eq.0) ok(jriver)=-1
               end if
            end do
            river_split(iriver)   = numcells
            flow_fraction(iriver) = flow_fraction_rel(iriver)/total_weight
         end do
         do iriver=1,nriver
           if (numcells.ne.1 .and. ok(iriver).ne.0) then
               LEVEL3 'Multicell river (',trim(river_name(iriver)),'):# ',iriver, &
                    'w=',flow_fraction(iriver)
            end if
         end do

         LEVEL3 'split:',river_split
!  Create a list with the river names without multiple-cells, i.e. just a
!  single entry per real river name.
         nn = 1
         ni = 1
         do n=1,nriver
            if (ni .le. nriver) then
               real_river_name(nn) = river_name(ni)
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

#ifdef GETM_BIO
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_rivers_bio
!
! !INTERFACE:
   subroutine init_rivers_bio()
!
! !DESCRIPTION:
! First, memory for storing the biological loads from rivers is
! allocated.
! The variable - {\tt river\_bio} - is initialised to  - {\tt bio\_missing}.
!
! !USES:
   IMPLICIT NONE
!
! !LOCAL VARIABLES:
   integer                   :: rc
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_rivers_bio() # ',Ncall
#endif

   LEVEL1 'init_rivers_bio()'

   allocate(river_bio(nriver,numc),stat=rc)
   if (rc /= 0) stop 'rivers: Error allocating memory (river_bio)'

   river_bio = bio_missing


#ifdef DEBUG
   write(debug,*) 'Leaving init_rivers_bio()'
   write(debug,*)
#endif
   return
   end subroutine init_rivers_bio
!EOC
#endif

#ifdef _FABM_
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_rivers_fabm
!
! !INTERFACE:
   subroutine init_rivers_fabm()
!
! !DESCRIPTION:
! First, memory for storing the biological loads from rivers is
! allocated.
! The variable - {\tt river\_fabm} - is initialised to  - variable-
! specific missing values obtained provided by FABM.
!
! !USES:
   IMPLICIT NONE
!
! !LOCAL VARIABLES:
   integer                   :: rc,m
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_rivers_fabm() # ',Ncall
#endif

   if (allocated(fabm_pel)) then
      LEVEL1 'init_rivers_fabm()'

      allocate(river_fabm(nriver,size(model%info%state_variables)),stat=rc)
      if (rc /= 0) stop 'rivers: Error allocating memory (river_fabm)'

      do m=1,size(model%info%state_variables)
         if (model%info%state_variables(m)%no_river_dilution) then
            river_fabm(:,m) = model%info%state_variables(m)%missing_value
         else
            river_fabm(:,m) = _ZERO_
         end if
      end do
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving init_rivers_fabm()'
   write(debug,*)
#endif
   return
   end subroutine init_rivers_fabm
!EOC
#endif

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_rivers - updating river points \label{sec-do-rivers}
!
! !INTERFACE:
   subroutine do_rivers(loop,do_3d)
!
! !DESCRIPTION:
!
! Here, the temperature, salinity, sea surface elevation and layer heights
! are updated in the river inflow grid boxes. Temperature and salinity
! are mixed with riverine values proportional to the old volume and the
! river inflow volume at that time step, sea surface elevation is simply
! increased by the inflow volume divided by the grid box area, and
! the layer heights are increased proportionally.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: loop
   logical, intent(in)                 :: do_3d
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,m,n
   REALTYPE                  :: rvol,height
   REALTYPE                  :: svol,tvol,vol
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_rivers() # ',Ncall
#endif

!  river spin-up
   if (ramp_is_active) then
      if (loop .ge. river_ramp) then
         ramp = _ONE_
         ramp_is_active = .false.
         STDERR LINE
         call write_time_string()
         LEVEL3 timestr,': finished river_ramp=',river_ramp
         STDERR LINE
      else
         ramp = _ONE_*loop/river_ramp
      end if
   end if

   select case (river_method)
      case(0)
      case(1,2)
         do n=1,nriver
            if(ok(n) .gt. 0) then
               i = ir(n)-ioff; j = jr(n)-joff
               rvol = ramp * dtm * river_flow(n) * flow_fraction(n)
               if (ok(n) .eq. 1) then
                  irr(n) = irr(n) + rvol
               end if
               height = rvol * ARCD1
               z(i,j) = z(i,j) + height
#ifndef NO_BAROCLINIC
               macro_height(n)=macro_height(n)+height
!              on macrotime step adjust 3d fields
               if (do_3d) then
                  if (macro_height(n).gt._ZERO_ .or. .not.river_outflow_properties_follow_source_cell) then
                     if (calc_salt ) then
                        if ( river_salt(n) .ne. salt_missing) then
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
#ifdef GETM_BIO
                     if (bio_calc) then
                        do m=1,numc
                           if ( river_bio(n,m) .ne. bio_missing ) then
                              cc3d(m,i,j,1:kmax) = &
                                    (cc3d(m,i,j,1:kmax)*(H(i,j)+ssen(i,j)) &
                                    + river_bio(n,m)*macro_height(n))      &
                                    / (H(i,j)+ssen(i,j)+macro_height(n))
                           end if
                        end do
                     end if
#endif
#ifdef _FABM_
                     if (allocated(fabm_pel)) then
                        do m=1,size(model%info%state_variables)
                           if ( river_fabm(n,m) .ne. model%info%state_variables(m)%missing_value ) then
                              fabm_pel(i,j,1:kmax,m) = &
                                    (fabm_pel(i,j,1:kmax,m)*(H(i,j)+ssen(i,j)) &
                                    + river_fabm(n,m)*macro_height(n))      &
                                    / (H(i,j)+ssen(i,j)+macro_height(n))
                           end if
                        end do
                     end if
#endif
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
         stop 'do_rivers'
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
! !IROUTINE:  clean_rivers
!
! !INTERFACE:
   subroutine clean_rivers
!
! !DESCRIPTION:
!
! This routine closes the river handling by writing the integrated
! river run-off for each river to standard output.
!
! !USES:
   IMPLICIT NONE
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
#ifdef _FABM_
         if (allocated(river_fabm)) deallocate(river_fabm)
#endif
      case default
         FATAL 'Not valid rivers_method specified'
         stop 'clean_rivers'
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving clean_rivers()'
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
