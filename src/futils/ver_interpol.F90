!$Id: ver_interpol.F90,v 1.2 2003-04-23 12:02:43 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ver_interpol - vertical interpolation.
!
! !INTERFACE:
   subroutine ver_interpol(nlev,zlev,prof,imin,jmin,imax,jmax,mask,H,  &
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
   integer, intent(in)                 :: nlev
   REALTYPE, intent(in)                :: zlev(nlev),prof(nlev)
   integer, intent(in)                 :: imin,imax,jmin,jmax
   integer, intent(in)                 :: mask(E2DFIELD)
   REALTYPE, intent(in)                :: H(E2DFIELD)
   integer, intent(in)                 :: iimin,iimax,jjmin,jjmax,kmax
   REALTYPE, intent(in)                :: hn(I3DFIELD)
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: field(I3DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ver_interpol.F90,v $
!  Revision 1.2  2003-04-23 12:02:43  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.1.1.1  2002/05/02 14:01:17  gotm
!  recovering after CVS crash
!
!  Revision 1.2  2001/08/29 11:21:46  bbh
!  namelists read in salinity and temperature + initialisation
!
!  Revision 1.1  2001/08/27 11:56:12  bbh
!  New file added
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
! !LOCAL VARIABLES:
   integer :: i,j,k,nn
   REALTYPE :: rat
   REALTYPE:: zmodel(kmax)
!
!EOP
!-----------------------------------------------------------------------
!BOC
   do j=jjmin,jjmax
      do i=iimin,iimax
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
