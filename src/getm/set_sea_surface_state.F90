#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_sea_surface_state -
!
! !INTERFACE:
   subroutine set_sea_surface_state( runtype , ssu , ssv )
!
! !USES:
   use domain      , only: imin,imax,jmin,jmax,kmax
   use domain      , only: az, grid_type, cosconv, sinconv 
   use domain      , only: xc, xu, xv, yc, yu, yv, dxv, dyu, arcd1
   use variables_2d, only: dtm, z, zo, Dvel, DU, DV, U, V
#ifndef NO_3D
   use variables_3d, only: dt, ho, hn, hvel, hun, hvn, uu, vv, ww
#ifndef NO_BAROCLINIC
   use variables_3d, only: T
#endif
#endif
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
   integer, intent(in)                        :: runtype
!
! !OUTPUT PARAMETERS:
   REALTYPE, dimension(E2DFIELD), intent(out) :: ssu, ssv
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES: 
   REALTYPE,dimension(E2DFIELD) :: wrk, u_vel, v_vel
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'set_sea_surface_state() # ',Ncall
#endif

   if (runtype .eq. 1) then
      wrk = _ZERO_
      call to_u(imin,jmin,imax,jmax,az,                            &
                dtm,grid_type,                                     &
                dxv,dyu,arcd1,                                     &
                xc,xu,xv,z,zo,Dvel,U,DU,V,DV,wrk,wrk,_ZERO_,u_vel)
      call to_v(imin,jmin,imax,jmax,az,                            &
                dtm,grid_type,                                     &
                dxv,dyu,arcd1,                                     &
                yc,yu,yv,z,zo,Dvel,U,DU,V,DV,wrk,wrk,_ZERO_,v_vel)
   else
#ifndef NO_3D
      call to_u(imin,jmin,imax,jmax,az,                                &
                dt,grid_type,                                          &
                dxv,dyu,arcd1,                                         &
                xc,xu,xv,hn(:,:,kmax),ho(:,:,kmax),hvel(:,:,kmax),     &
                uu(:,:,kmax),hun(:,:,kmax),vv(:,:,kmax),hvn(:,:,kmax), &
                ww(:,:,kmax-1),ww(:,:,kmax),_ZERO_,u_vel)
      call to_v(imin,jmin,imax,jmax,az,                                &
                dt,grid_type,                                          &
                dxv,dyu,arcd1,                                         &
                yc,yu,yv,hn(:,:,kmax),ho(:,:,kmax),hvel(:,:,kmax),     &
                uu(:,:,kmax),hun(:,:,kmax),vv(:,:,kmax),hvn(:,:,kmax), &
                ww(:,:,kmax-1),ww(:,:,kmax),_ZERO_,v_vel)
#ifndef NO_BAROCLINIC
      if (runtype .gt. 2) then
         !sst = T(:,:,kmax)
      end if
#endif
#endif
   end if

!  rotate back to grid-related coordinates
   ssu = cosconv*u_vel - sinconv*v_vel
   ssv = sinconv*u_vel + cosconv*v_vel

#ifdef DEBUG
   write(debug,*) 'Leaving set_sea_surface_state()'
   write(debug,*)
#endif
   return
   end subroutine set_sea_surface_state
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2015 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
