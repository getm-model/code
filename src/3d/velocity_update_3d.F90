#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: velocity_update_3d - calculate new 3D velocities.
!
! !INTERFACE:
   subroutine velocity_update_3d(calc_2d)
!
! !DESCRIPTION:
!
! !USES:
   use domain
   use m2d,          only: velocity_update
   use variables_2d, only: Uint,Vint
   use variables_3d, only: kmin,dt
   use variables_3d, only: hn,ho,uu,vv,ww,hun,hvn,hvel,velx3d,vely3d,w
   use variables_3d, only: ssen,sseo,Dun,Dvn,Dveln,velx2dadv,vely2dadv
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical,intent(in) :: calc_2d
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
   integer            :: k
   REALTYPE,parameter :: vel_missing=-9999.0
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'velocity_update_3d() # ',Ncall
#endif

   do k=1,kmax
      call velocity_update(dt,hn(:,:,k),ho(:,:,k),hvel(:,:,k),         &
                           uu(:,:,k),hun(:,:,k),vv(:,:,k),hvn(:,:,k),  &
                           wwm=ww(:,:,k-1),wwp=ww(:,:,k),              &
                           velx=velx3d(:,:,k),vely=vely3d(:,:,k))
   end do
   call to_w(imin,jmin,imax,jmax,kmin,kmax,az,                         &
             dt,                                                       &
             dxv,dyu,arcd1,                                            &
             H,HU,HV,hn,ho,hvel,uu,hun,vv,hvn,ww,vel_missing,w)

   if (calc_2d) then
      call velocity_update(dt,ssen,sseo,Dveln,Uint,Dun,Vint,Dvn,       &
                           velx=velx2dadv,vely=vely2dadv)
   end if


#ifdef DEBUG
   write(debug,*) 'Leaving velocity_update_3d()'
   write(debug,*)
#endif
   return
   end subroutine velocity_update_3d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2016 - Knut Klingbeil                                  !
!-----------------------------------------------------------------------
