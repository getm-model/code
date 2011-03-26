#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: part_domain() - partition the domain
!
! !INTERFACE:
   subroutine part_domain()
!
! !DESCRIPTION:
!  Set various integers defining the calculation domain. The settings
!  depends on STATIC vs. DYNAMIC compilation and serial vs. parallel
!  model run.
!
! !USES:
   use domain, only: iextr,jextr
   use domain, only: imin,imax,jmin,jmax,kmax
   use domain, only: ioff,joff
#ifdef GETM_PARALLEL
   use halo_mpi, only: part_domain_mpi
#endif
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'part_domain'
#endif

#ifndef GETM_PARALLEL
#ifndef STATIC
   imin=1 ; imax=iextr ; jmin=1 ; jmax=jextr
#endif
   ioff=0 ; joff=0
#else
   call part_domain_mpi(iextr,jextr,kmax,imin,imax,jmin,jmax,ioff,joff)
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving part_domain()'
   write(debug,*)
#endif
   return
   end subroutine part_domain
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2002 - Karsten Bolding and Hans Burchard (BBH)         !
!-----------------------------------------------------------------------
