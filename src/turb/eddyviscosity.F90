!$Id: eddyviscosity.F90,v 1.1.1.1 2002-05-02 14:01:57 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: eddyviscosity() - calculates the eddy-viscosity. 
! 
! !INTERFACE:
   subroutine eddyviscosity  
!
! !DESCRIPTION: 
!
! !USES:
   use domain, only: iimin,iimax,jjmin,jjmax,H,min_depth,crit_depth
   use meteo,  only: tausx,tausy
   use m3d,    only: kmax,kmin,hn,ssen,num,nuh,taub
   use m3d,    only: kappa
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
!  $Log: eddyviscosity.F90,v $
!  Revision 1.1.1.1  2002-05-02 14:01:57  gotm
!  recovering after CVS crash
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
! 
! !LOCAL VARIABLES:
   integer	:: i,j,k   
   REALTYPE	:: zz,DD,dry,taus 
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'eddyviscosity() # ',Ncall
#endif

   do j=jjmin,jjmax 
      do i=iimin,iimax 
         zz=0 
         DD=ssen(i,j)+H(i,j) 
         dry=max(_ZERO_,min(_ONE_,(DD-min_depth)/(crit_depth-min_depth)))
         taus=sqrt(tausx(i,j)**2+tausy(i,j)**2) 
         do k=kmin(i,j),kmax-1
            zz=zz+hn(i,j,k) 
            num(i,j,k)=kappa*zz*(1-zz/DD)*sqrt(max(taus,taub(i,j)))  &
                      +(1-dry)*1.e-2 
         end do 
      end do 
   end do  

#ifdef DEBUG
   write(debug,*) 'Leaving eddyviscosity()'
   write(debug,*)
#endif
   return
   end subroutine eddyviscosity
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
