!$Id: coriolis.F90,v 1.1 2002-05-02 14:00:42 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: coriolis - 2D coriolis terms.
!
! !INTERFACE:
   subroutine coriolis(D1,vel,avgvel,D,cor)
!
! !DESCRIPTION:
!  This routine calculates the Corilis term - in both x and y direction
!  on a Arakawa C-grid.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,CORI
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)	:: D1(E2DFIELD)
   REALTYPE, intent(in)	:: vel(E2DFIELD)
   REALTYPE, intent(in)	:: D(E2DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out):: avgvel(E2DFIELD)
   REALTYPE, intent(out):: cor(E2DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: coriolis.F90,v $
!  Revision 1.1  2002-05-02 14:00:42  gotm
!  Initial revision
!
!  Revision 1.1.1.1  2001/04/17 08:43:07  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
   integer	:: i,j
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'Coriolis() # ',Ncall
#endif

   do j=jmin,jmax
      do i=imin,imax
         avgvel(i,j)=0.25*( vel(i  ,j  )/D1(i  ,j  )  &
                           +vel(i+1,j  )/D1(i+1,j  )  &
                           +vel(i  ,j-1)/D1(i  ,j-1)  &
                           +vel(i+1,j-1)/D1(i+1,j-1) )
         cor(i,j)=Cori*D(i,j)*avgvel(i,j)
      end do
   end do
#ifdef DEBUG
   write(debug,*) 'Leaving coriolis()'
   write(debug,*)
#endif
   return
   end subroutine coriolis
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
