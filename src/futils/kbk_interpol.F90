!$Id: kbk_interpol.F90,v 1.1.1.1 2002-05-02 14:01:21 gotm Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: kbk_interpol - vertical interpolation.
!
! !INTERFACE:
   subroutine kbk_interpol(nlev,zlev,salt_z,imin,jmin,imax,jmax,mask,H,	&
                           iimin,jjmin,iimax,jjmax,kmax,hn,field)
!
! !DESCRIPTION:
!
!  This is a utility subroutine in which observational data which might
!  be given on an arbitrary, but ordered grid, are interpolated and
!  extrapolated to the actual (moving) model grid.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in) 	:: nlev
   integer, intent(in)	:: iimin,iimax,jjmin,jjmax,kmax
   REALTYPE, intent(in)	:: zlev(nlev),salt_z(iimin:iimax,jjmin:jjmax,nlev)
   integer, intent(in)	:: imin,imax,jmin,jmax
   integer, intent(in)	:: mask(E2DFIELD)
   REALTYPE, intent(in)	:: H(E2DFIELD)
   REALTYPE, intent(in)	:: hn(I3DFIELD)
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out):: field(I3DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: kbk_interpol.F90,v $
!  Revision 1.1.1.1  2002-05-02 14:01:21  gotm
!  recovering after CVS crash
!
!
! !LOCAL VARIABLES:
   integer 	:: i,j,k,nn,kl
   REALTYPE 	:: rat
   REALTYPE	:: zmodel(kmax),prof(nlev)
!
!EOP
!-----------------------------------------------------------------------
!BOC
!STDERR salt_z(75,72,:)
!stop
   do j=jjmin,jjmax
      do i=iimin,iimax

         if(mask(i,j) .ne. 0) then
#if 1
            do k=1,nlev
               prof(k) = salt_z(i,j,nlev-k+1)
            end do
#else
            prof = salt_z(i,j,:)
#endif

            zmodel(1) = -H(i,j) + 0.5*hn(i,j,1)
            do k=2,kmax
               zmodel(k) = zmodel(k-1) + 0.5*(hn(i,j,k-1)+hn(i,j,k))
            end do
!STDERR zmodel
!STDERR zlev
!stop

!  Set surface values to uppermost input value
            do k=kmax,1,-1
               if (zmodel(k) .ge. zlev(nlev)) then
                  field(i,j,k) = prof(nlev)
               end if
            end do

!  Search index of lowest OK value
            do k=1,nlev
               if(prof(k) .ne. -9999.0) EXIT     
            end do
	    kl = k
!  Set bottom values to lowest input value
            do k=1,kmax
               if (zmodel(k) .le. zlev(kl)) then
                  field(i,j,k) = prof(kl)
               end if
            end do

!  Interpolate inner values linearly
            do k=1,kmax
               if ((zmodel(k) .lt. zlev(nlev)) .and. (zmodel(k) .gt. zlev(kl))) then
                  nn=0
224               nn=nn+1
                  if (zlev(nn) .le. zmodel(k)) goto 224
                  rat=(zmodel(k)-zlev(nn-1))/(zlev(nn)-zlev(nn-1))
                  field(i,j,k)=(1-rat)*prof(nn-1)+rat*prof(nn)
               end if
            end do
            field(i,j,0) = field(i,j,1)
#if 1
if( i .eq. 85 .and. j .eq. 18) then
STDERR kl
do k=nlev,1,-1
STDERR zlev(k),prof(k)
end do
do k=kmax,1,-1
STDERR zmodel(k),field(i,j,k)
end do
!stop
endif
#endif
         else
            field(i,j,:) = _ZERO_
         end if
      end do
   end do

   return
   end subroutine kbk_interpol
!EOC

!-----------------------------------------------------------------------
!Copyright (C) 2001 - Karsten Bolding and Hans Burchard
!-----------------------------------------------------------------------
