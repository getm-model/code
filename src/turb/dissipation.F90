!$Id: dissipation.F90,v 1.1 2002-05-02 14:01:56 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: dissipation() - calculates turbulent dissipation. 
! 
! !INTERFACE:
   subroutine dissipation 
!
! !DESCRIPTION: 
!
! !USES:
   use parameters, only: kappa
   use domain,     only: iimin,iimax,jjmin,jjmax,kmax,az
   use m2d,        only: zub,zvb
   use m3d,        only: kmin,hn,ho,num
   use m3d,        only: dt,taub,eps,P,B,tke,tkeo
   use turb,       only: ce1,ce2,ce3,z0s,cde,sige,epsmin
   IMPLICIT NONE
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
!  $Log: dissipation.F90,v $
!  Revision 1.1  2002-05-02 14:01:56  gotm
!  Initial revision
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
   integer	:: i,j,k
   REALTYPE	:: prod,buoy,diss,dtl,cee3 
   REALTYPE	:: pplus(1:kmax-1),pminus(1:kmax-1),dif(1:kmax),aux(1:kmax) 
   REALTYPE	:: a1(0:kmax),a2(0:kmax),a3(0:kmax),a4(0:kmax),res(0:kmax) 
! 
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'dissipation() # ',Ncall
#endif

   do i=iimin,iimax 
      do j=jjmin,jjmax 
         dtl=dt*az(i,j)   ! Applying land mask to time step 

         do k=kmin(i,j)+1,kmax 
            dif(k)=0.5*(num(i,j,k-1)+num(i,j,k))/sige 
         end do 
         dif(kmin(i,j))=taub(i,j)**2			&
		/(0.5*(eps(i,j,kmin(i,j))+eps(i,j,kmin(i,j)-1)))/sige  
         do k=kmin(i,j),kmax 
            aux(k)=dtl*dif(k)/hn(i,j,k) 
         end do

         do k=kmin(i,j),kmax-1 
            prod=ce1*eps(i,j,k)/tkeo(i,j,k)*P(i,j,k) 
            if (B(i,j,k).gt.0) then 
               cee3=1. 
            else 
               cee3=ce3 
            end if 
            buoy=cee3*eps(i,j,k)/tkeo(i,j,k)*B(i,j,k) 
            diss=ce2*eps(i,j,k)/tkeo(i,j,k)*eps(i,j,k) 
            if (prod+buoy.gt.0) then 
               pplus(k)=prod+buoy 
               pminus(k)=diss 
            else 
               pplus(k)=prod 
               pminus(k)=diss-buoy 
            end if 
         end do 

         do k=kmin(i,j),kmax-1 
            a3(k)=-2.*aux(k+1)/(hn(i,j,k+1)+hn(i,j,k)) 
            a1(k)=-2.*aux(k  )/(hn(i,j,k+1)+hn(i,j,k)) 
            a2(k)=1-a1(k)-a3(k)+dtl*pminus(k)/eps(i,j,k) 
            a4(k)=(ho(i,j,k+1)+ho(i,j,k))/(hn(i,j,k+1)+hn(i,j,k))*eps(i,j,k)  &
                 +dtl*pplus(k) 
         end do 

         a3(kmin(i,j)-1)=0 
         a2(kmin(i,j)-1)=1. 
         a4(kmin(i,j)-1)=cde*(tke(i,j,kmin(i,j)-1)**1.5)/kappa                &
                        /(0.25*(zub(i-1,j)+zub(i,j)+zvb(i,j-1)+zvb(i,j))) 

         a2(kmax)=1. 
         a1(kmax)=0 
         a4(kmax)=cde*(tke(i,j,kmax)**1.5)/kappa/z0s 

         call getm_tridiagonal(kmax,kmin(i,j)-1,kmax,a1,a2,a3,a4,res) 

         do k=kmin(i,j)-1,kmax 
            eps(i,j,k)=max(epsmin,res(k)) 
         end do 

      end do 
   end do 

#ifdef DEBUG
   write(debug,*) 'Leaving dissipation()'
   write(debug,*)
#endif
   return
   end subroutine dissipation
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
