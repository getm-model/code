!JMB
         subroutine ztoh(zpos,hn,depthmin)
#include "cppdefs.h"
!
! !DESCRIPTION:
!
! !USES:
   use domain,   only: iimin,iimax,jjmin,jjmax,kmax,H
   IMPLICIT NONE
   REALTYPE        :: zpos(I3DFIELD),hn(I3DFIELD),depthmin
   integer         :: i,j,k
      do k=1,kmax
      do j=jjmin,jjmax
       do i=iimin,iimax
       hn(i,j,k)= zpos(i,j,k)-zpos(i,j,k-1)
       hn(i,j,k)=max(hn(i,j,k),depthmin)
       enddo
       enddo
      enddo
! End Back to layer thickness
     return
     end
         subroutine htoz(hn,zpos)
#include "cppdefs.h"
!
! !DESCRIPTION:
!
! !USES:
   use domain,   only: iimin,iimax,jjmin,jjmax,kmax,H
   IMPLICIT NONE
   REALTYPE        :: zpos(I3DFIELD),hn(I3DFIELD),depthmin
   integer         :: i,j,k
!     write(6,*) 'htoz',iimax,hn(iimax/2,2,kmax/2),H(iimax/2,2)
      do j=jjmin,jjmax
       do i=iimin,iimax
       zpos(i,j,0)=-H(i,j)
       do k=1,kmax
       zpos(i,j,k)=zpos(i,j,k-1)+hn(i,j,k)
       enddo
       enddo
      enddo
     return
     end
         subroutine hcheck(hn,ssen,h)
#include "cppdefs.h"
!
! !DESCRIPTION:
!
! !USES:
   use domain,   only: iimin,iimax,jjmin,jjmax,kmax
   IMPLICIT NONE
   REALTYPE        :: ssen(I2DFIELD),hn(I3DFIELD),h(I2DFIELD),HH,depthmin
   integer         :: i,j,k
! Final check of layer thicnkess thoug not necessary if zpos treated correctly
!     write(6,*) 'Inside',hn(iimax/2,2,kmax/2)
      do j=jjmin,jjmax
         do i=iimin,iimax
	   HH=0.
           do k=1,kmax
	   HH=HH+hn(i,j,k)
           end do
           do k=1,kmax
           hn(i,j,k)=hn(i,j,k)* (ssen(i,j)+H(i,j))/HH
	   end do
         end do
      end do
!     write(6,*) 'Inside after',hn(iimax/2,2,kmax/2)
   
     return
     end
     
