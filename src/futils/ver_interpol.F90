#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ver_interpol - vertical interpolation.
!
! !INTERFACE:
   subroutine ver_interpol(nlev,zlev,prof,imin,jmin,imax,jmax,kmax, &
                           mask,H,hn,field)
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
   integer, intent(in)                 :: nlev
   REALTYPE, intent(in)                :: zlev(nlev),prof(nlev)
   integer, intent(in)                 :: imin,jmin,imax,jmax,kmax
   integer, intent(in)                 :: mask(E2DFIELD)
   REALTYPE, intent(in)                :: H(E2DFIELD)
   REALTYPE, intent(in)                :: hn(I3DFIELD)
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: field(I3DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer :: i,j,k,nn
   REALTYPE :: rat
   REALTYPE:: zmodel(kmax)
!EOP
!-----------------------------------------------------------------------
!BOC
   do j=jmin,jmax
      do i=imin,imax
         if(mask(i,j) .ne. 0) then

            zmodel(1) = -H(i,j) + 0.5*hn(i,j,1)
            do k=2,kmax
               zmodel(k) = zmodel(k-1) + 0.5*(hn(i,j,k-1)+hn(i,j,k))
            end do

!  Set surface values to uppermost input value
            do k=kmax,1,-1
               if (zmodel(k) .ge. zlev(nlev)) then
                  field(i,j,k) = prof(nlev)
               end if
            end do

!  Set bottom values to lowest input value
            do k=1,kmax
               if (zmodel(k) .le. zlev(1)) then
                  field(i,j,k) = prof(1)
               end if
            end do

!  Interpolate inner values linearly
            do k=1,kmax
               if ((zmodel(k) .lt. zlev(nlev)) .and. (zmodel(k) .gt. zlev(1))) then
                  nn=0
224               nn=nn+1
                  if (zlev(nn) .le. zmodel(k)) goto 224
                  rat=(zmodel(k)-zlev(nn-1))/(zlev(nn)-zlev(nn-1))
                  field(i,j,k)=(1-rat)*prof(nn-1)+rat*prof(nn)
               end if
            end do
            field(i,j,0) = field(i,j,1)
         else
            field(i,j,:) = _ZERO_
         end if
      end do
   end do

   return
   end
!EOC

!-----------------------------------------------------------------------
!Copyright (C) 2001 - Karsten Bolding and Hans Burchard
!-----------------------------------------------------------------------
