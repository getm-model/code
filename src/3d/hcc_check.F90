!$Id: hcc_check.F90,v 1.1 2004-05-04 09:23:51 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: hcc_check() - hydrostatic consistency criteria
!
! !INTERFACE:
   subroutine hcc_check()
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,az,au,av,HU,HV
   use domain, only: iimin,iimax,jjmin,jjmax,kmax
   use variables_3d, only: hn,hun,hvn,hcc
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: hcc_check.F90,v $
!  Revision 1.1  2004-05-04 09:23:51  kbk
!  hydrostatic consistency criteria stored in .3d.nc file
!
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   REALTYPE                  :: du1,du2,dv1,dv2
   REALTYPE                  :: x,y
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'HCC check'

   do j=jjmin,jjmax
      do i=iimin,iimax
         if (az(i,j) .ge. 1) then
            du1=HU(i-1,j)
            du2=HU(i,j)
            dv1=HV(i,j-1)
            dv2=HV(i,j)

            do k=1,kmax

               if (au(i-1,j) .ge. 1 .and. au(i,j) .ge. 1) then
                  du1 = du1+0.5*hun(i-1,j,k)
                  du2 = du2+0.5*hun(i,j,k)
               else
                  du1 = _ZERO_
                  du2 = _ZERO_
               end if
               if (av(i,j-1) .ge. 1 .and. av(i,j) .ge. 1) then
                  dv1 = dv1+0.5*hvn(i,j-1,k)
                  dv2 = dv2+0.5*hvn(i,j,k)
               else
                  dv1 = _ZERO_
                  dv2 = _ZERO_
               end if

               x = (du2-du1)/hn(i,j,k) 
               y = (dv2-dv1)/hn(i,j,k)
               hcc(i,j,k) = max(abs(x),abs(y))

               if (au(i-1,j) .ge. 1 .and. au(i,j) .ge. 1) then
                  du1 = du1+0.5*hun(i-1,j,k)
                  du2 = du2+0.5*hun(i,j,k)
               else
                  du1 = _ZERO_
                  du2 = _ZERO_
               end if
               if (av(i,j-1) .ge. 1 .and. av(i,j) .ge. 1) then
                  dv1 = dv1+0.5*hvn(i,j-1,k)
                  dv2 = dv2+0.5*hvn(i,j,k)
               else
                  dv1 = _ZERO_
                  dv2 = _ZERO_
               end if
            end do
         end if
      end do
   end do

   return
   end subroutine hcc_check
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
