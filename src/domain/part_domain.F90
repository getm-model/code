!$Id: part_domain.F90,v 1.2 2007-06-07 10:25:19 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: part_domain - partition the calculation domain
!
! !INTERFACE:
   subroutine part_domain()
!
! !USES:
   use domain, only: iextr,jextr
   use domain, only: imin,imax,jmin,jmax,kmax
   use domain, only: ioff,joff
#ifdef PARALLEL
   use halo_mpi, only: part_domain_mpi
#endif
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Set various integers defining the calculation domain
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'part_domain'
#endif

#ifndef PARALLEL
#ifndef STATIC
   imin=1 ; imax=iextr ; jmin=1 ; jmax=jextr
#endif
   ioff=0 ; joff=0
#else
   call part_domain_mpi(iextr,jextr,kmax,imin,imax,jmin,jmax,ioff,joff)
#endif

#if 0
#ifndef STATIC
      if(il .eq. -1 .or. ih .eq. -1 .or. jl .eq. -1 .or. jh .eq. -1) then
         imin = 1 ; imax = iextr ; jmin = 1 ; jmax = jextr;
         il = imin ; il = imax ; jl = jmin ; jh = jmax
      else
         imin = 1 ; imax = ih-il+1 ; jmin = 1 ; jmax = jh-jl+1;
      end if
#endif
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
