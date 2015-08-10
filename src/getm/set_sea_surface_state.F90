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
   use domain      , only: az, grid_type, cosconv, sinconv 
   use domain      , only: xc, xu, xv, yc, yu, yv, dxv, dyu, arcd1
   use domain      , only: NWB, wi, wfj, wlj
   use domain      , only: NNB, nj, nfi, nli
   use domain      , only: NEB, ei, efj, elj
   use domain      , only: NSB, sj, sfi, sli
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
   REALTYPE,dimension(E2DFIELD)        :: wrk, u_vel, v_vel
   REALTYPE,dimension(:,:),pointer     :: p_U, p_V, p_Dvel, p2d
   integer,dimension(:,:),allocatable,target,save  :: mask
   integer,dimension(:,:),pointer,save :: p_mask
   integer                             :: i,j,n,start, rc
   logical,save                        :: first=.true.
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'set_sea_surface_state() # ',Ncall
#endif

   if (first) then
      if (grid_type.eq.1 .or. grid_type.eq.2) then
         allocate(mask(E2DFIELD))
         if (rc /= 0) stop 'set_sea_surface_state: Error allocating memory (mask)'
         where(az .eq. 0)
            mask = 0
         elsewhere
            mask = 1
         end where
         p_mask => mask
      else
         p_mask => az
      end if
      first = .false.
   end if

   if (runtype .eq. 1) then
      wrk = _ZERO_
      call to_u(imin,jmin,imax,jmax,p_mask,                        &
                dtm,grid_type,                                     &
                dxv,dyu,arcd1,                                     &
                xc,xu,xv,z,zo,Dvel,U,DU,V,DV,wrk,wrk,_ZERO_,u_vel)
      call to_v(imin,jmin,imax,jmax,p_mask,                        &
                dtm,grid_type,                                     &
                dxv,dyu,arcd1,                                     &
                yc,yu,yv,z,zo,Dvel,U,DU,V,DV,wrk,wrk,_ZERO_,v_vel)
      p_U    => U
      p_V    => V
      p_Dvel => Dvel
   else
#ifndef NO_3D
      if (.not. do_3d) return
      call to_u(imin,jmin,imax,jmax,p_mask,                            &
                dt,grid_type,                                          &
                dxv,dyu,arcd1,                                         &
                xc,xu,xv,hn(:,:,kmax),ho(:,:,kmax),hvel(:,:,kmax),     &
                uu(:,:,kmax),hun(:,:,kmax),vv(:,:,kmax),hvn(:,:,kmax), &
                ww(:,:,kmax-1),ww(:,:,kmax),_ZERO_,u_vel)
      call to_v(imin,jmin,imax,jmax,p_mask,                            &
                dt,grid_type,                                          &
                dxv,dyu,arcd1,                                         &
                yc,yu,yv,hn(:,:,kmax),ho(:,:,kmax),hvel(:,:,kmax),     &
                uu(:,:,kmax),hun(:,:,kmax),vv(:,:,kmax),hvn(:,:,kmax), &
                ww(:,:,kmax-1),ww(:,:,kmax),_ZERO_,v_vel)
      p2d => uu  (:,:,kmax) ; p_U   (imin-HALO:,jmin-HALO:) => p2d
      p2d => vv  (:,:,kmax) ; p_V   (imin-HALO:,jmin-HALO:) => p2d
      p2d => hvel(:,:,kmax) ; p_Dvel(imin-HALO:,jmin-HALO:) => p2d
#ifndef NO_BAROCLINIC
      if (runtype .gt. 2) then
         !sst = T(:,:,kmax)
      end if
#endif
#endif
   end if


   if (grid_type.eq.1 .or. grid_type.eq.2) then

      ssu = u_vel
      ssv = v_vel

   else

!     rotate back to grid-related coordinates
      ssu = cosconv*u_vel - sinconv*v_vel
      ssv = sinconv*u_vel + cosconv*v_vel

!     dirty approximation for open bdy cells
      do n = 1,NWB
         i = wi(n)
         start = max(jmin-HALO+1,wfj(n))
         do j = start,wlj(n)
            ssu(i,j) = _HALF_*( p_U(i-1,j  ) + p_U(i,j) )/p_Dvel(i,j)
            ssv(i,j) = _HALF_*( p_V(i  ,j-1) + p_V(i,j) )/p_Dvel(i,j)
         end do
      end do
      do n = 1,NNB
         j = nj(n)
         start = max(imin-HALO+1,nfi(n))
         do i = start,nli(n)
            ssu(i,j) = _HALF_*( p_U(i-1,j  ) + p_U(i,j) )/p_Dvel(i,j)
            ssv(i,j) = _HALF_*( p_V(i  ,j-1) + p_V(i,j) )/p_Dvel(i,j)
         end do
      end do
      do n = 1,NEB
         i = ei(n)
         start = max(jmin-HALO+1,efj(n))
         do j = start,elj(n)
            ssu(i,j) = _HALF_*( p_U(i-1,j  ) + p_U(i,j) )/p_Dvel(i,j)
            ssv(i,j) = _HALF_*( p_V(i  ,j-1) + p_V(i,j) )/p_Dvel(i,j)
         end do
      end do
      do n = 1,NSB
         j = sj(n)
         start = max(imin-HALO+1,sfi(n))
         do i = start,sli(n)
            ssu(i,j) = _HALF_*( p_U(i-1,j  ) + p_U(i,j) )/p_Dvel(i,j)
            ssv(i,j) = _HALF_*( p_V(i  ,j-1) + p_V(i,j) )/p_Dvel(i,j)
         end do
      end do

   end if

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
