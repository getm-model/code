!$Id: tke_calc.F90,v 1.1 2002-05-02 14:01:57 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: tke_calc() - 
! 
! !INTERFACE:
   subroutine tke_calc 
!
! !DESCRIPTION: 
!
! !USES:
   use domain, only: iimin,iimax,jjmin,jjmax,kmax,az
   use m3d,    only: dt,kmin,hn,ho,num,taub,P,B,eps,tke,tkeo 
   use turb,   only: sigk,cmue,cde,tkemin 
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
!  $Log: tke_calc.F90,v $
!  Revision 1.1  2002-05-02 14:01:57  gotm
!  Initial revision
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
   integer	:: i,j,k
   REALTYPE	:: prod,buoy,diss,dtl 
   REALTYPE	:: pplus(1:kmax-1),pminus(1:kmax-1),dif(1:kmax),aux(1:kmax) 
   REALTYPE	:: a1(0:kmax),a2(0:kmax),a3(0:kmax),a4(0:kmax),Res(0:kmax) 
! 
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'tke_calc() # ',Ncall
#endif

   do i=iimin,iimax 
      do j=jjmin,jjmax 
         do k=kmin(i,j)-1,kmax 
            tkeo(i,j,k)=tke(i,j,k) 
         end do 
      end do 
   end do 
   do i=iimin,iimax 
      do j=jjmin,jjmax 

         dtl=dt*az(i,j)   ! Applying land mask to time step 

         do k=kmin(i,j)+1,kmax  
            dif(k)=0.5*(num(i,j,k-1)+num(i,j,k))/sigk 
         end do
         dif(kmin(i,j))=taub(i,j)**2			&
	               /(0.5*(eps(i,j,kmin(i,j))+eps(i,j,kmin(i,j)-1)))/sigk  
         do k=kmin(i,j),kmax 
            aux(k)=dtl*dif(k)/hn(i,j,k) 
         end do 

         do k=kmin(i,j),kmax-1 
            prod=P(i,j,k) 
            buoy=B(i,j,k) 
            diss=eps(i,j,k) 
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
            a2(k)=1-a1(k)-a3(k)+dtl*pminus(k)/tke(i,j,k) 
            a4(k)=(ho(i,j,k+1)+ho(i,j,k))/(hn(i,j,k+1)+hn(i,j,k)) &
                 *tke(i,j,k)+dtl*pplus(k) 
         end do 

         a3(kmin(i,j)-1)=0 
         a2(kmin(i,j)-1)=1. 
         a4(kmin(i,j)-1)=taub(i,j)/dsqrt(cmue*cde) 

!hb        a2(kmax)=1. 
!hb        a1(kmax)=0 
!hb        a4(kmax)=dsqrt(tausx**2+tausy**2)/dsqrt(cmue*cde) 
!hb         a2(kmax)=1. 
!hb        a1(kmax)=-1.  
!hb        a4(kmax)=0 
         k=kmax 
         a1(k)=-aux(k  )/hn(i,j,k) 
         a2(k)=1-a1(k) 
         a4(k)=ho(i,j,k)/hn(i,j,k)*tke(i,j,k) 

         call getm_tridiagonal(kmax,kmin(i,j)-1,kmax,a1,a2,a3,a4,Res) 

         do k=kmin(i,j)-1,kmax 
            tke(i,j,k)=max(tkemin,Res(k)) 
         end do 

      end do 
   end do 

#ifdef DEBUG
   write(debug,*) 'Leaving tke_calc()'
   write(debug,*)
#endif
   return
   end subroutine tke_calc
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
