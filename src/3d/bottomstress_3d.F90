!$Id: bottomstress_3d.F90,v 1.1 2002-05-02 14:00:53 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: bottomstress_3d() - calculates the bottom stresses.
!
! !INTERFACE:
   subroutine bottomstress_3d
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: iimin,iimax,jjmin,jjmax
   use variables_3d,    only: kumin,kvmin,uu,vv,hun,hvn,rru,rrv,taub
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
!  $Log: bottomstress_3d.F90,v $
!  Revision 1.1  2002-05-02 14:00:53  gotm
!  Initial revision
!
!  Revision 1.3  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.2  2001/05/03 20:12:31  bbh
!  Use of variables_3d
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
   integer	:: i,j,k1,k2,k3,k4
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'bottomstress_3d() # ',Ncall
#endif

   do j=jjmin,jjmax    ! Absolute Value of Bottom Friction
      do i=iimin,iimax
         k1=kumin(i-1,j  )
         k2=kumin(i  ,j  )
         k3=kvmin(i  ,j-1)
         k4=kvmin(i  ,j  )
         taub(i,j)=0.5*(                                            &
                   (uu(i-1,j  ,k1)/hun(i-1,j  ,k1)*rru(i-1,j  ))**2 &
                  +(uu(i  ,j  ,k2)/hun(i  ,j  ,k2)*rru(i  ,j  ))**2 &
                  +(vv(i  ,j-1,k3)/hvn(i  ,j-1,k3)*rrv(i  ,j-1))**2 &
                  +(vv(i  ,j  ,k4)/hvn(i  ,j  ,k4)*rrv(i  ,j  ))**2)
         taub(i,j)=sqrt(taub(i,j))
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving bottomstress_3d()'
   write(debug,*)
#endif
   return
   end subroutine bottomstress_3d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
