!$Id: sealevel.F90,v 1.13 2007-01-10 21:28:55 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: sealevel - using the cont. eq. to get the sealevel.
!
! !INTERFACE:
   subroutine sealevel
!
! !DESCRIPTION:
!
! Here, the sea surface elevation is iterated according to the vertically
! integrated continuity equation given in (\ref{Elevation}) on page
! \pageref{Elevation}. 
!
! When working with the option {\tt SLICE\_MODEL}, the elevations
! at $j=2$ are copied to $j=3$.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,az,H
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only : arcd1,dxv,dyu
#else
   use domain, only : dx,dy,ard1
#endif
   use m2d, only: dtm
   use variables_2d, only: z,zo,U,V
   use halo_zones, only : update_2d_halo,wait_halo,z_TAG
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: i,j
   REALTYPE                  :: kk
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'sealevel() # ',Ncall
#endif

   zo = z
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1) then
            z(i,j)=z(i,j)-dtm*((U(i,j)*DYU-U(i-1,j  )*DYUIM1) &
                              +(V(i,j)*DXV-V(i  ,j-1)*DXVJM1))*ARCD1

#ifdef FRESHWATER_LENSE_TEST
       kk=1.0
       if ((((i.eq.1).or.(i.eq.imax)).and.(j.ge.1).and.(j.le.jmax)).or. &
           (((j.eq.1).or.(j.eq.jmax)).and.(i.ge.1).and.(i.le.imax)))    &
          z(i,j)=(1.-kk)*z(i,j)
       kk=0.5625
       if ((((i.eq.2).or.(i.eq.imax-1)).and.(j.ge.2).and.(j.le.jmax-1)).or. &
           (((j.eq.2).or.(j.eq.jmax-1)).and.(i.ge.2).and.(i.le.imax-1)))    &
          z(i,j)=(1.-kk)*z(i,j)
       kk=0.25
       if ((((i.eq.3).or.(i.eq.imax-2)).and.(j.ge.3).and.(j.le.jmax-2)).or. &
           (((j.eq.3).or.(j.eq.jmax-2)).and.(i.ge.3).and.(i.le.imax-2)))    &
           z(i,j)=(1.-kk)*z(i,j)
       kk=0.0625
       if ((((i.eq.4).or.(i.eq.imax-3)).and.(j.ge.4).and.(j.le.jmax-3)).or. &
           (((j.eq.4).or.(j.eq.jmax-3)).and.(i.ge.4).and.(i.le.imax-3)))    &
           z(i,j)=(1.-kk)*z(i,j)
#endif
         end if
      end do
   end do

   call sealevel_nan_check()

#ifdef SLICE_MODEL
      do i=imin,imax
         z(i,3)=z(i,2)
      end do
#endif

   call update_2d_halo(z,z,az,imin,jmin,imax,jmax,z_TAG)
   call wait_halo(z_TAG)

#ifdef DEBUG
   write(debug,*) 'Leaving sealevel()'
   write(debug,*)
#endif
   return
   end subroutine sealevel
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: sealevel_nan_check - Sweep the sealevel (z) for NaN values
!
! !INTERFACE:
   subroutine sealevel_nan_check
!
! !DESCRIPTION:
!
!  The sea surface elevation (2d) variable is sweeped scanning for 
!  not-a-number (NaN). NaN values indicate that the integration 
!  has become unstable and that it really should be stopped.
!  First time a NaN value is found, a warning is issued and possibly
!  the code is stopped. After the first encounter, the sweep is 
!  suspended.
!
!  The behaviour of this routine is controlled by the 
!  {\tt sealevel\_check} parameter in the {\tt m2d} namelist.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
   use m2d, only: sealevel_check
   use variables_2d, only: z
   use exceptions, only: getm_error
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Bjarne B\"uchmann
!
! !LOCAL VARIABLES:
   integer, save :: Ncall = 0
   integer, save :: can_check   = 0
   integer, save :: have_warned = 0
   integer       :: num_nan
   integer       :: i,j,inan,jnan, idum
   REALTYPE      :: ahuge,zdum
!
!EOP
!-----------------------------------------------------------------------
!BOC
   Ncall = Ncall+1
#ifdef DEBUG
   write(debug,*) 'sealevel_nan_check() # ',Ncall
#endif

! Caveat! Programmers beware.
!
! In general checking for NaN is not trivial.
! Some compilers allow the use of functions ISNAN or ieee_is_nan.
! The functions are not neceassily available for all compilers, so be careful.
! The approach taken here is to compare against HUGE and 
! then say that the value is *not* NaN if it can compare less than HUGE.
! However, it is important that the comparison is not removed by the 
! compiler - as it may really test out as "always true" from a compile-
! time point of view. Function inlining and other funny stuff done by smart
! compilers can make this approach not working.
!

! Fast return if we should not sweep:
  if (sealevel_check.eq.0 .or. have_warned.gt.0) then
#ifdef DEBUG
     write(debug,*) 'Leaving sealevel_nan_check() early'
     write(debug,*)
#endif
     return
  end if

!
! On first call we test if the present method works with the present
! compilation.
!
   ahuge = HUGE(ahuge)
!
! See if we can check with the present compiler+settings
   if (Ncall .eq. 1) then
! Check a regular case
      zdum = 1.0
      call sealevel_nandum(zdum,ahuge,idum)
      if (idum .eq. 1) then
         can_check = can_check + 1
      end if
! Check a failing (NaN) case:
      zdum = 0.0/0.0
      call sealevel_nandum(zdum,ahuge,idum)
      if (idum .eq. 2) then
         can_check = can_check + 1
      end if
! Now we must have two successes - otherwise we cannot check.
      if (can_check .eq. 2) then
         LEVEL1 "NaN sweeping seems OK with this compiler and settings"
      else
         can_check = 0
         LEVEL1 "WARNING NaN sweeping not OK with this compiler and/or compile settings"
      end if
   end if
!
! Main bulk of work:
   if (can_check .gt. 0 .and. mod(Ncall,abs(sealevel_check)).eq.0) then
! Count number of NaNs encountered
      num_nan = 0
      do j=jmin,jmax
         do i=imin,imax
            call sealevel_nandum(z(i,j),ahuge,idum)
            if (idum .eq. 2) then
! Increment counter for NaNs:
               num_nan = num_nan + 1
! Keep at least one point with location (may be overwritten later)
               inan = i
               jnan = j
            end if
         end do
      end do
! Do something if there were NaNs:
      if (num_nan .gt. 0) then
! Warn about this badness
         LEVEL1 "WARNING: NaN occuring in elevations at call no. ",Ncall
         LEVEL1 ": Number of NaN values: ",num_nan
         LEVEL1 ": NaNs include i j = ",inan,jnan
         have_warned = 1
! Stop execution if we have to (otherwise just warn)
         if (sealevel_check .gt. 0) then
            call getm_error("catch_nan_2d()","NaN values encountered in elevations")
         end if
      end if
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving sealevel_nan_check()'
   write(debug,*)
#endif
   return
   end subroutine sealevel_nan_check
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: sealevel_nandum. Helper routine to spot NaNs
!
! !INTERFACE:
   subroutine sealevel_nandum(a,b,idum)
!
! !USES:
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)       :: a,b
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   integer, intent(out)      :: idum
!
!
! !DESCRIPTION:
!  This routine is a kind of dummy routine primarily to provide a means
!  to spot NaN  values.
!  Output is 1 or 2, based on which is smaller (a or b, respectively).
!  The default is 2, and the idea is that "imin=2" should be returned
!  also if a is NaN. If b=HUGE(b), then this provides a means to detect 
!  if a is a denormal number.
!
!
!EOP
!-----------------------------------------------------------------------
!BOC
   idum = 2
   if ( b .gt. a) then
      idum = 1
   end if

   return
   end subroutine sealevel_nandum
!EOC


!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
