!$Id: production.F90,v 1.1.1.1 2002-05-02 14:01:57 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: production() - 
! 
! !INTERFACE:
   subroutine production 
!
! !DESCRIPTION: 
!
! !USES:
   use domain, only: iimin,iimax,jjmin,jjmax,kmax
   use m3d,    only: kmin,kumin,uu,hun,kvmin,vv,hvn,num,P,B
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
!  $Log: production.F90,v $
!  Revision 1.1.1.1  2002-05-02 14:01:57  gotm
!  recovering after CVS crash
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
   integer	:: i,j,k,rc 
   REALTYPE	:: dzu(iimin-1:iimax,jjmin  :jjmax,1:kmax-1) 
   REALTYPE	:: dzv(iimin  :iimax,jjmin-1:jjmax,1:kmax-1) 
! 
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'production() # ',Ncall
#endif

   do j=jjmin,jjmax 
      do i=iimin,iimax-1 
         do k=1,kumin(i,j)-1 
            dzu(i,j,k)=0 
         end do 
         do k=kumin(i,j),kmax-1 
            dzu(i,j,k)=uu(i,j,k+1)/hun(i,j,k+1)-uu(i,j,k)/hun(i,j,k) 
            dzu(i,j,k)=dzu(i,j,k)/(0.5*(hun(i,j,k+1)+hun(i,j,k))) 
         end do 
      end do 
      do k=kumin(iimin-1,j),kmax-1 
         dzu(iimin-1,j,k)=dzu(iimin  ,j,k) 
      end do 
      do k=kumin(iimax,j),kmax-1 
         dzu(iimax  ,j,k)=dzu(iimax-1,j,k) 
      end do 
   end do 
   do i=iimin,iimax 
      do j=jjmin,jjmax-1 
         do k=1,kvmin(i,j)-1 
            dzv(i,j,k)=0 
         end do 
         do k=kvmin(i,j),kmax-1 
            dzv(i,j,k)=vv(i,j,k+1)/hvn(i,j,k+1)-vv(i,j,k)/hvn(i,j,k) 
            dzv(i,j,k)=dzv(i,j,k)/(0.5*(hvn(i,j,k+1)+hvn(i,j,k))) 
         end do 
      end do 
      do k=kvmin(i,jjmin-1),kmax-1 
         dzv(i,jjmin-1,k)=dzv(i,jjmin  ,k) 
      end do 
      do k=kvmin(i,jjmax),kmax-1 
         dzv(i,jjmax  ,k)=dzv(i,jjmax-1,k) 
      end do 
   end do
   do i=iimin+1,iimax-1 
      do j=jjmin+1,jjmax-1 
         do k=1,kmin(i,j)-1 
            P(i,j,k)=0 
            B(i,j,k)=0 
         end do 
         do k=kmin(i,j),kmax-1
            P(i,j,k)=0.5*(						& 
                0.5*(num(i-1,j  ,k)+num(i  ,j  ,k))*dzu(i-1,j  ,k)**2+	& 
                0.5*(num(i  ,j  ,k)+num(i+1,j  ,k))*dzu(i  ,j  ,k)**2+	& 
                0.5*(num(i  ,j-1,k)+num(i  ,j  ,k))*dzv(i  ,j-1,k)**2+	& 
                0.5*(num(i  ,j  ,k)+num(i  ,j+1,k))*dzv(i  ,j  ,k)**2)
            B(i,j,k)=0 
         end do 
      end do 
   end do 
   i=iimin 
   do j=jjmin+1,jjmax-1 
      do k=kmin(i,j)+1,kmax-1 
         P(i,j,k)=0.5*(							& 
                                 num(i  ,j  ,k) *dzu(i-1,j  ,k)**2+ 	&
             0.5*(num(i  ,j  ,k)+num(i+1,j  ,k))*dzu(i  ,j  ,k)**2+ 	&
             0.5*(num(i  ,j-1,k)+num(i  ,j  ,k))*dzv(i  ,j-1,k)**2+ 	&
             0.5*(num(i  ,j  ,k)+num(i  ,j+1,k))*dzv(i  ,j  ,k)**2) 
      end do 
   end do 
   i=iimax  
   do j=jjmin+1,jjmax-1 
      do k=kmin(i,j)+1,kmax-1 
         P(i,j,k)=0.5*(							&
             0.5*(num(i-1,j  ,k)+num(i  ,j  ,k))*dzu(i-1,j  ,k)**2+ 	&
                  num(i  ,j  ,k)                *dzu(i  ,j  ,k)**2+ 	&
             0.5*(num(i  ,j-1,k)+num(i  ,j  ,k))*dzv(i  ,j-1,k)**2+ 	&
             0.5*(num(i  ,j  ,k)+num(i  ,j+1,k))*dzv(i  ,j  ,k)**2) 
      end do 
   end do 
   j=jjmin 
   do i=iimin+1,iimax-1 
      do k=kmin(i,j)+1,kmax-1 
         P(i,j,k)=0.5*(							&
             0.5*(num(i-1,j  ,k)+num(i  ,j  ,k))*dzu(i-1,j  ,k)**2+	&
             0.5*(num(i  ,j  ,k)+num(i+1,j  ,k))*dzu(i  ,j  ,k)**2+	&
                                 num(i  ,j  ,k) *dzv(i  ,j-1,k)**2+	&
             0.5*(num(i  ,j  ,k)+num(i  ,j+1,k))*dzv(i  ,j  ,k)**2)
      end do
   end do
   j=jjmax
   do i=iimin+1,iimax-1
      do k=kmin(i,j)+1,kmax-1
         P(i,j,k)=0.5*(							&
             0.5*(num(i-1,j  ,k)+num(i  ,j  ,k))*dzu(i-1,j  ,k)**2+	&
             0.5*(num(i  ,j  ,k)+num(i+1,j  ,k))*dzu(i  ,j  ,k)**2+	&
             0.5*(num(i  ,j-1,k)+num(i  ,j  ,k))*dzv(i  ,j-1,k)**2+	&
                  num(i  ,j  ,k)                *dzv(i  ,j  ,k)**2)
      end do
   end do
   i=iimin
   j=jjmin
   do k=kmin(i,j)+1,kmax-1
      P(i,j,k)=0.5*(							&
                              num(i  ,j  ,k) *dzu(i-1,j  ,k)**2+	&
          0.5*(num(i  ,j  ,k)+num(i+1,j  ,k))*dzu(i  ,j  ,k)**2+	&
                              num(i  ,j  ,k) *dzv(i  ,j-1,k)**2+	&
          0.5*(num(i  ,j  ,k)+num(i  ,j+1,k))*dzv(i  ,j  ,k)**2)
   end do
   i=iimax
   j=jjmin
   do k=kmin(i,j)+1,kmax-1
      P(i,j,k)=0.5*(							&
          0.5*(num(i-1,j  ,k)+num(i  ,j  ,k))*dzu(i-1,j  ,k)**2+	&
               num(i  ,j  ,k)                *dzu(i  ,j  ,k)**2+	&
                              num(i  ,j  ,k) *dzv(i  ,j-1,k)**2+	&
          0.5*(num(i  ,j  ,k)+num(i  ,j+1,k))*dzv(i  ,j  ,k)**2)
   end do
   i=iimin
   j=jjmax
   do k=kmin(i,j)+1,kmax-1
      P(i,j,k)=0.5*(							&
                              num(i  ,j  ,k) *dzu(i-1,j  ,k)**2+	&
          0.5*(num(i  ,j  ,k)+num(i+1,j  ,k))*dzu(i  ,j  ,k)**2+	&
          0.5*(num(i  ,j-1,k)+num(i  ,j  ,k))*dzv(i  ,j-1,k)**2+	&
               num(i  ,j  ,k)                *dzv(i  ,j  ,k)**2)
   end do
   i=iimax
   j=jjmax
   do k=kmin(i,j)+1,kmax-1
      P(i,j,k)=0.5*(							&
          0.5*(num(i-1,j  ,k)+num(i  ,j  ,k))*dzu(i-1,j  ,k)**2+	&
               num(i  ,j  ,k)                *dzu(i  ,j  ,k)**2+	&
          0.5*(num(i  ,j-1,k)+num(i  ,j  ,k))*dzv(i  ,j-1,k)**2+	&
               num(i  ,j  ,k)                *dzv(i  ,j  ,k)**2)
   end do  

#ifdef DEBUG
   write(debug,*) 'Leaving production()'
   write(debug,*)
#endif
   return
   end subroutine production
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
