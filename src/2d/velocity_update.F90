#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: velocity_update - calculate new velocities.
!
! !INTERFACE:
   subroutine velocity_update(dt,z,zo,Dvel,U,DU,V,DV,wwm,wwp,missing,  &
                              velx,vely)

!  Note (KK): keep in sync with interface in m2d.F90
!
! !DESCRIPTION:
!
! !USES:
   use domain
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                      :: dt
   REALTYPE,dimension(E2DFIELD),intent(in)  :: z,zo,Dvel,U,DU,V,DV
   REALTYPE,dimension(E2DFIELD),target,intent(in),optional :: wwm,wwp
   REALTYPE,intent(in),target,optional      :: missing
!
! !OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(out) :: velx,vely
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
   REALTYPE,dimension(E2DFIELD)        :: scalefac
   REALTYPE,dimension(E2DFIELD),target :: zeros
   REALTYPE,dimension(:,:),pointer     :: p_wwm,p_wwp
   REALTYPE                            :: missval
   REALTYPE,parameter                  :: vel_missing=-9999.0
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'velocity_update() # ',Ncall
#endif

   zeros = _ZERO_
   if (present(wwm)) then
      p_wwm => wwm
   else
      p_wwm => zeros
   end if
   if (present(wwp)) then
      p_wwp => wwp
   else
      p_wwp => zeros
   end if
   if (present(missing)) then
      missval = missing
   else
      missval = vel_missing
   end if

   if (grid_type .eq. 3) scalefac = _ONE_
   if (grid_type .eq. 4) scalefac = deg2rad*rearth*cos(deg2rad*latc)
   call to_u(imin,jmin,imax,jmax,az,                                   &
             dt,grid_type,                                             &
             dxv,dyu,arcd1,                                            &
             p_xc,p_xu,p_xv,z,zo,Dvel,U,DU,V,DV,p_wwm,p_wwp,scalefac,missval,velx)

   if (grid_type .eq. 4) scalefac = deg2rad*rearth
   call to_v(imin,jmin,imax,jmax,az,                                   &
             dt,grid_type,                                             &
             dxv,dyu,arcd1,                                            &
             p_yc,p_yu,p_yv,z,zo,Dvel,U,DU,V,DV,p_wwm,p_wwp,scalefac,missval,vely)

#ifdef DEBUG
   write(debug,*) 'Leaving velocity_update()'
   write(debug,*)
#endif
   return
   end subroutine velocity_update
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2016 - Knut Klingbeil                                  !
!-----------------------------------------------------------------------
