!$Id: turb.F90,v 1.1 2002-05-02 14:01:57 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: turb - turbulens model component
!
! !INTERFACE:
   module turb
!
! !DESCRIPTION: 
!  This modules contains declarations for all variables related to turbulens
!  calculations. Information about the calculation domain
!  is included from the \emph{domain.F90} module.
!  The module contains public subroutines for initialisation, integration 
!  and clean up of the turbulens model component.
!  The actual calculation routines are called in Integrate3D from 
!  \emph{m3d.F90} and is linked in from the library libturb.a.
!
! !USES:
   use parameters, only: kappa
   use domain,     only: iimin,iimax,jjmin,jjmax,kmax
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
#ifdef KAJ_KURT
#ifndef STATIC
   REALTYPE, dimension(:,:,:), allocatable	:: tke,tkeo,eps,P,B
#else
   REALTYPE	:: tke(I3DFIELD)
   REALTYPE	:: tkeo(I3DFIELD)
   REALTYPE	:: eps(I3DFIELD)
   REALTYPE	:: P(I3DFIELD)
   REALTYPE	:: B(I3DFIELD)
#endif
#endif
   REALTYPE, parameter :: tkemin=1.0e-7
   REALTYPE, parameter :: epsmin=5.0e-10
   REALTYPE, parameter 	:: cmue=0.5625d0,sigk=1.0
   REALTYPE, parameter :: ce1=1.44,ce2=1.92,ce3=-0.9
   REALTYPE, parameter :: z0s=0.1
   REALTYPE, parameter :: cde=cmue*cmue*cmue
   REALTYPE, parameter :: sige=kappa**2*cmue/(ce2-ce1)/cde

!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: turb.F90,v $
!  Revision 1.1  2002-05-02 14:01:57  gotm
!  Initial revision
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_turb - initialise turbulens related stuff.
!
! !INTERFACE:
   subroutine init_turb(hotstart)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical, intent(in):: hotstart
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Allocates memory for turbulens related fields.
!
! !REVISION HISTORY:
!
!  22Apr99   Karsten Bolding & Hans Burchard  Initial code.
!
! !LOCAL VARIABLES:
   integer 	:: rc
!
!EOP
!-------------------------------------------------------------------------
!BOC

#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_turb() # ',Ncall
#endif

   LEVEL2 'init_turb'

!kbk   READ(NAMLST,turb)
!kbk   REWIND(NAMLST)

#ifdef KAJ_KURT
#ifndef STATIC
! Allocates memory for the public data members
   allocate(tke(I3DFIELD),stat=rc)
   if (rc /= 0) stop 'init_turb: Error allocating memory (tke)'

   allocate(tkeo(I3DFIELD),stat=rc)
   if (rc /= 0) stop 'init_turb: Error allocating memory (tkeo)'

   allocate(eps(I3DFIELD),stat=rc)
   if (rc /= 0) stop 'init_turb: Error allocating memory (eps)'

   allocate(P(I3DFIELD),stat=rc)
   if (rc /= 0) stop 'init_turb: Error allocating memory (P)'

   allocate(B(I3DFIELD),stat=rc)
   if (rc /= 0) stop 'init_turb: Error allocating memory (B)'
#endif

   if (hotstart) then
   else
      tke=tkemin
      eps=epsmin
   end if
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving init_turb()'
   write(debug,*)
#endif
   return
   end subroutine init_turb
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clean_turb - cleanup after Turb run.
!
! !INTERFACE:
   subroutine clean_turb()
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This routine cleans up after a Turb integration. Close open files etc.
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'clean_turb() # ',Ncall
#endif

! Deallocates memory for the public data members

#ifdef DEBUG
     write(debug,*) 'Leaving clean_turb()'
     write(debug,*)
#endif
   return
   end subroutine clean_turb
!EOC

!-----------------------------------------------------------------------

   end module turb

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
