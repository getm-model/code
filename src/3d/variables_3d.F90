!$Id: variables_3d.F90,v 1.1 2002-05-02 14:00:58 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: variables_3d - global 3D related variables
!
! !INTERFACE:
   module variables_3d
!
! !DESCRIPTION:
!  This modules contains declarations for all variables related to 3D
!  hydrodynamical calculations. Information about the calculation domain
!  is included from the \emph{domain.F90} module.
!  The module contains public subroutines to initialise and cleanup.
!
! !USES:
   use domain,     only: iimin,iimax,jjmin,jjmax,kmax
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   REALTYPE	:: dt,cnpar=0.9
!
#ifdef STATIC
#include "static.h"
#else
#include "dynamic_declarations.h"
#endif
   integer	:: size3d_field
   integer	:: mem3d
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: variables_3d.F90,v $
!  Revision 1.1  2002-05-02 14:00:58  gotm
!  Initial revision
!
!  Revision 1.6  2001/09/19 13:07:00  bbh
!  Moved advection related 3D fields to global allocation
!
!  Revision 1.5  2001/09/01 17:10:25  bbh
!  Vertical coordinate definition now specified via namelist
!
!  Revision 1.4  2001/08/27 11:51:45  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.3  2001/05/21 13:07:19  bbh
!  dt and cnpar is in variables_3d.F90
!
!  Revision 1.2  2001/05/18 08:25:52  bbh
!  Added zooming variables
!
!  Revision 1.1  2001/05/03 19:31:56  bbh
!  3D variables seperated from m3d
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_variables_3d - initialise 3D relatedstuff.
!
! !INTERFACE:
   subroutine init_variables_3d(runtype)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)		:: runtype
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Allocates memiory for 3D related fields.
!
! !REVISION HISTORY:
!
!  See log for the module.
!
! !LOCAL VARIABLES:
   integer	:: rc
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_variables_3d() # ',Ncall
#endif

   LEVEL2 'init_variables_3d'
   size3d_field=((iimax+HALO)-(iimin+HALO)+1)* 			&
                ((jjmax+HALO)-(jjmin+HALO)+1)*(kmax+1)
   mem3d=n3d_fields*size3d_field*REAL_SIZE

!  Allocates memory for the public data members - if not static
#ifndef STATIC
#include "dynamic_allocations.h"
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving init_variables_3d()'
   write(debug,*)
#endif
   return
   end subroutine init_variables_3d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clean_variables_3d - cleanup after 3D run.
!
! !INTERFACE:
   subroutine clean_variables_3d()
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This routine cleans up after a 3D integration. Close open files etc.
!
! !REVISION HISTORY:
!  See log for the module.
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'clean_3d() # ',Ncall
#endif

! Deallocates memory for the public data members

#ifdef DEBUG
     write(debug,*) 'Leaving clean_variables_3d()'
     write(debug,*)
#endif
   return
   end subroutine clean_variables_3d
!EOC

!-----------------------------------------------------------------------

   end module variables_3d

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
