#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_sea_surface_state -
!
! !INTERFACE:
   subroutine set_sea_surface_state( runtype , ssu , ssv , do_3d )
!
! !USES:
   use domain      , only: imin,imax,jmin,jmax,kmax
   use domain      , only: grid_type, cosconv, sinconv
   use domain      , only: NWB, wi, wfj, wlj
   use domain      , only: NNB, nj, nfi, nli
   use domain      , only: NEB, ei, efj, elj
   use domain      , only: NSB, sj, sfi, sli
   use variables_2d, only: Dvel, U, V, velx, vely
#ifndef NO_3D
   use variables_3d, only: hvel, uu, vv, velx3d, vely3d
#ifndef NO_BAROCLINIC
   use variables_3d, only: T
#endif
#endif
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
   integer, intent(in)                        :: runtype
   logical, intent(in)                        :: do_3d
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
   REALTYPE,dimension(:,:),pointer :: p_U, p_V, p_Dvel,p_velx,p_vely,p2d
   integer                         :: i,j,n,start
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'set_sea_surface_state() # ',Ncall
#endif

   if (runtype .eq. 1) then
      p_U    => U
      p_V    => V
      p_Dvel => Dvel
      p_velx => velx
      p_vely => vely
   else
#ifndef NO_3D
      if (.not. do_3d) return
      p2d => uu    (:,:,kmax) ; p_U   (imin-HALO:,jmin-HALO:) => p2d
      p2d => vv    (:,:,kmax) ; p_V   (imin-HALO:,jmin-HALO:) => p2d
      p2d => hvel  (:,:,kmax) ; p_Dvel(imin-HALO:,jmin-HALO:) => p2d
      p2d => velx3d(:,:,kmax) ; p_velx(imin-HALO:,jmin-HALO:) => p2d
      p2d => vely3d(:,:,kmax) ; p_vely(imin-HALO:,jmin-HALO:) => p2d
#ifndef NO_BAROCLINIC
      if (runtype .gt. 2) then
         !sst = T(:,:,kmax)
      end if
#endif
#endif
   end if


   if (grid_type.eq.1 .or. grid_type.eq.2) then

      ssu = p_velx
      ssv = p_vely

   else

!     rotate back to grid-related coordinates
      ssu = cosconv*p_velx - sinconv*p_vely
      ssv = sinconv*p_velx + cosconv*p_vely

   end if

!  dirty approximation for open bdy cells
!  KK-TODO: Replace extra handling of open bdy cells
!           by valid setting of ssu,ssv(az=2)!
!           Calculation of velx,vely(az=2) seems to require ww(az=2) and
!           _MIRROR_BDY_EXTRA_. But continuity is illposed with mirrored
!           transports!!!
   do n = 1,NWB
      i = wi(n)
      start = max(jmin-HALO+1,wfj(n))
      do j = start,wlj(n)
         ssu(i,j) = p_U(i,j)/p_Dvel(i,j)
         ssv(i,j) = _HALF_*( p_V(i  ,j-1) + p_V(i,j) )/p_Dvel(i,j)
      end do
   end do
   do n = 1,NNB
      j = nj(n)
      start = max(imin-HALO+1,nfi(n))
      do i = start,nli(n)
         ssu(i,j) = _HALF_*( p_U(i-1,j  ) + p_U(i,j) )/p_Dvel(i,j)
         ssv(i,j) = p_V(i,j-1)/p_Dvel(i,j)
      end do
   end do
   do n = 1,NEB
      i = ei(n)
      start = max(jmin-HALO+1,efj(n))
      do j = start,elj(n)
         ssu(i,j) = p_U(i-1,j)/p_Dvel(i,j)
         ssv(i,j) = _HALF_*( p_V(i  ,j-1) + p_V(i,j) )/p_Dvel(i,j)
      end do
   end do
   do n = 1,NSB
      j = sj(n)
      start = max(imin-HALO+1,sfi(n))
      do i = start,sli(n)
         ssu(i,j) = _HALF_*( p_U(i-1,j  ) + p_U(i,j) )/p_Dvel(i,j)
         ssv(i,j) = p_V(i,j)/p_Dvel(i,j)
      end do
   end do

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
