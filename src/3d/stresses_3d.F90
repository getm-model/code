!$Id: stresses_3d.F90,v 1.4 2003-08-14 09:48:09 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: stresses_3d() - calculates the bottom stresses.
!
! !INTERFACE:
   subroutine stresses_3d
!
! !DESCRIPTION:
!
! !USES:
   use parameters, only: rho_0
   use domain, only: az,au,av,iimin,iimax,jjmin,jjmax
   use variables_3d, only: kumin,kvmin,uu,vv,hun,hvn,rru,rrv,taus,taub
   use meteo, only: tausx,tausy
   use halo_zones, only : update_2d_halo,wait_halo,z_TAG
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
!  $Log: stresses_3d.F90,v $
!  Revision 1.4  2003-08-14 09:48:09  kbk
!  need rru,rrv in halo zones
!
!  Revision 1.3  2003/04/23 12:16:34  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 13:36:38  kbk
!  parallel support, cleaned code + NO_3D, NO_BAROCLINIC
!
!  Revision 1.1.1.1  2002/05/02 14:00:56  gotm
!  recovering after CVS crash
!
!  Revision 1.4  2001/10/26 09:11:28  bbh
!  Stresses in meteo.F90 are in N/m2 - divide by rho_0 where necessary
!
!  Revision 1.3  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.2  2001/05/03 20:22:44  bbh
!  Also uses variables_3d
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k1,k2,k3,k4
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'stresses_3d() # ',Ncall
#endif

!  we need to know rru and rrv in the halos as well
   call update_2d_halo(rru,rru,au,iimin,jjmin,iimax,jjmax,10)
   call wait_halo(10)
   call update_2d_halo(rrv,rrv,av,iimin,jjmin,iimax,jjmax,10)
   call wait_halo(10)

   do j=jjmin,jjmax    ! Absolute Value of Bottom Friction
      do i=iimin,iimax
!KBK         if(az(i,j) .ge. 1) then
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
         taus(i,j)=0.5*(                                &
                        tausx(i,j)**2+tausx(i-1,j)**2   &
                      + tausy(i,j)**2+tausy(i,j-1)**2)
         taus(i,j)=sqrt(taus(i,j))/rho_0
!KBK         end if
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving stresses_3d()'
   write(debug,*)
#endif
   return
   end subroutine stresses_3d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
