!$Id: pran_kol.F90,v 1.1 2002-05-02 14:01:57 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: pran_kol() - using the cont. eq. to get the sealevel. 
! 
! !INTERFACE:
   subroutine pran_kol 
!
! !DESCRIPTION: 
!
! !USES:
   use domain, only: iimin,iimax,jjmin,jjmax,kmax,min_depth,crit_depth
   use m2d,    only: D
   use m3d,    only: kmin,num,tke,eps
   use turb,   only: cde,cmue
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
!  $Log: pran_kol.F90,v $
!  Revision 1.1  2002-05-02 14:01:57  gotm
!  Initial revision
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
   REALTYPE	:: dry
   integer	:: i,j,k  
! 
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'pran_kol() # ',Ncall
#endif
!
   do i=iimin,iimax 
      do j=jjmin,jjmax 
         dry=max(_ZERO_,min(_ONE_,(D(i,j)-min_depth)/(crit_depth-min_depth))) 
         do k=kmin(i,j),kmax-1 
            num(i,j,k)=cmue*cde*tke(i,j,k)**2/eps(i,j,k)+(1-dry)*1.e-2 
         end do 
      end do 
   end do 
!
#ifdef DEBUG
   write(debug,*) 'Leaving pran_kol()'
   write(debug,*)
#endif
   return
   end subroutine pran_kol
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
