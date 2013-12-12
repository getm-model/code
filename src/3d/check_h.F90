!JMB
         subroutine ztoh(zpos,hn,depthmin)
#include "cppdefs.h"
!
! !DESCRIPTION:
!
! !USES:
   use domain,   only: imin,imax,jmin,jmax,kmax,H
   IMPLICIT NONE
   REALTYPE        :: zpos(I3DFIELD),hn(I3DFIELD),depthmin
   integer         :: i,j,k
      do k=1,kmax
      do j=jmin-HALO,jmax+HALO
       do i=imin-HALO,imax+HALO
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
   use domain,   only: imin,imax,jmin,jmax,kmax,H
   IMPLICIT NONE
   REALTYPE        :: zpos(I3DFIELD),hn(I3DFIELD),depthmin
   integer         :: i,j,k
!     write(6,*) 'htoz',imax,hn(imax/2,2,kmax/2),H(imax/2,2)
      do j=jmin-HALO,jmax+HALO
       do i=imin-HALO,imax+HALO
       zpos(i,j,0)=-H(i,j)
       do k=1,kmax
       zpos(i,j,k)=zpos(i,j,k-1)+hn(i,j,k)
       enddo
       enddo
      enddo
     return
     end
         subroutine hcheck(hn,Dn)
#include "cppdefs.h"
!
! !DESCRIPTION:
!
! !USES:
   use domain,   only: imin,imax,jmin,jmax,kmax
   IMPLICIT NONE
   REALTYPE        :: Dn(I2DFIELD),hn(I3DFIELD)
   REALTYPE        :: HH,depthmin
   integer         :: i,j,k
! Final check of layer thicnkess thoug not necessary if zpos treated correctly
!     write(6,*) 'Inside',hn(imax/2,2,kmax/2)
      do j=jmin,jmax
         do i=imin,imax
            HH=0.
           do k=1,kmax
              HH=HH+hn(i,j,k)
           end do
           do k=1,kmax
              hn(i,j,k)=hn(i,j,k)* Dn(i,j)/HH
           end do
         end do
      end do
!     write(6,*) 'Inside after',hn(imax/2,2,kmax/2)

     return
     end

