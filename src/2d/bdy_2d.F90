#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  bdy_2d - 2D boundary conditions \label{bdy-2d}
!
! !INTERFACE:
   module bdy_2d
!
! !DESCRIPTION:
!
! Here, the two-dimensional boundary
! conditions for sea surface elevation and transports are handled.
!
! !USES:
   use halo_zones, only : z_TAG,H_TAG,U_TAG,V_TAG
   use domain, only: imin,jmin,imax,jmax,kmax,H,az,au,av
   use domain, only: nsbv,NWB,NNB,NEB,NSB,bdy_index,bdy_2d_type
   use domain, only: wi,wfj,wlj,nj,nfi,nli,ei,efj,elj,sj,sfi,sli
   use domain, only: min_depth
   use variables_2d, only: dtm,z,D,U,DU,V,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxc,dyc
#else
   use domain, only: dx,dy
#endif

   IMPLICIT NONE
!
   private
!
! !PUBLIC DATA MEMBERS:
   public init_bdy_2d, do_bdy_2d
   integer,public                   :: bdyfmt_2d,bdyramp_2d=-1
!  KK-TODO: static REAL_4B => allocatable REALTYPE
   REAL_4B,dimension(1500),public :: bdy_data,bdy_data_u,bdy_data_v
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_bdy_2d - initialising 2D boundary conditions
! \label{sec-init-bdy-2d}
!
! !INTERFACE:
   subroutine init_bdy_2d()
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !LOCAL VARIABLES:
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_bdy_2d() # ',Ncall
#endif

   LEVEL2 'init_bdy_2d()'

#ifdef DEBUG
   write(debug,*) 'Leaving init_bdy_2d()'
   write(debug,*)
#endif
   return
   end subroutine init_bdy_2d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_bdy_2d  - updating 2D boundary conditions
! \label{sec-do-bdy-2d}
!
! !INTERFACE:
   subroutine do_bdy_2d(loop,tag)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: loop,tag
!
! !INPUT/OUTPUT PARAMETERS:
!
! !LOCAL VARIABLES:
   REALTYPE                  :: a,fac
   integer                   :: i,j,k,l,n
   REALTYPE, parameter       :: FOUR=4.*_ONE_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_bdy_2d() # ',Ncall
#endif

   select case (bdyfmt_2d)
      case (NO_DATA)
      case (ANALYTICAL)
      case (ASCII)
      case (NETCDF)
!        Read in get_2d_bdy() via get_2d_bdy_ncdf()
      case default
         FATAL 'A non valid communication method has been chosen'
         stop 'update2dbdy'
   end select

!  Data read - do time interpolation

   fac = _ONE_
   if(bdyramp_2d .gt. 1) fac=min( _ONE_ ,FOUR*loop/float(bdyramp_2d))

   l = 0
   do n = 1,NWB
      l = l+1
      k = bdy_index(l)
      i = wi(n)
      do j = wfj(n),wlj(n)
         select case (bdy_2d_type(l))
            case (ZERO_GRADIENT)
               z(i,j) = z(i+1,j)
            case (SOMMERFELD)
!              KK-TODO: change DXC to DXU ?!
!                       change D(i,j) to _HALF_*(D(i,j)+D(i+1,j)) ?
               z(i,j) = z(i,j) + dtm/DXC*sqrt(9.81*D(i,j))*(z(i+1,j)-z(i,j))
            case (CLAMPED)
               z(i,j) = max(fac*bdy_data(k),-H(i,j)+min_depth)
            case (FLATHER_ELEV)
               a = sqrt(DU(i,j)/9.81)*(U(i,j)/DU(i,j)-bdy_data_u(k))
               z(i,j) = max(fac*(bdy_data(k) - a),-H(i,j)+min_depth)
            case default
               FATAL 'Illegal NWB 2D boundary type selection'
               stop 'update_2d_bdy()'
         end select
         k = k+1
      end do
   end do

   do n = 1,NNB
      l = l+1
      k = bdy_index(l)
      j = nj(n)
      do i = nfi(n),nli(n)
         select case (bdy_2d_type(l))
            case (ZERO_GRADIENT)
               z(i,j) = z(i,j-1)
            case (SOMMERFELD)
!              KK-TODO: change DYC to DYVJM1 ?! (not yet in cppdefs.h!)
!                       change D(i,j) to _HALF_*(D(i,j-1)+D(i,j)) ?
               z(i,j) = z(i,j) - dtm/DYC*sqrt(9.81*D(i,j))*(z(i,j)-z(i,j-1))
            case (CLAMPED)
               z(i,j) = max(fac*bdy_data(k),-H(i,j)+min_depth)
            case (FLATHER_ELEV)
               a = sqrt(DV(i,j)/9.81)*(V(i,j-1)/DV(i,j-1)-bdy_data_v(k))
               z(i,j) = max(fac*(bdy_data(k) + a),-H(i,j)+min_depth)
            case default
               FATAL 'Illegal NNB 2D boundary type selection'
               stop 'update_2d_bdy()'
         end select
         k = k+1
      end do
   end do

   do n = 1,NEB
      l = l+1
      k = bdy_index(l)
      i = ei(n)
      do j = efj(n),elj(n)
         select case (bdy_2d_type(l))
            case (ZERO_GRADIENT)
               z(i,j) = z(i-1,j)
            case (SOMMERFELD)
!              KK-TODO: change DXC to DXUIM1 ?! (not yet in cppdefs.h!)
!                       change D(i,j) to _HALF_*(D(i-1,j)+D(i,j)) ?
               z(i,j) = z(i,j) - dtm/DXC*sqrt(9.81*D(i,j))*(z(i,j)-z(i-1,j))
            case (CLAMPED)
               z(i,j) = max(fac*bdy_data(k),-H(i,j)+min_depth)
            case (FLATHER_ELEV)
               a = sqrt(DU(i,j)/9.81)*(U(i-1,j)/DU(i-1,j)-bdy_data_u(k))
               z(i,j) = max(fac*(bdy_data(k) + a),-H(i,j)+min_depth)
            case default
               FATAL 'Illegal NEB 2D boundary type selection'
               stop 'update_2d_bdy()'
         end select
         k = k+1
      end do
   end do

   do n = 1,NSB
      l = l+1
      k = bdy_index(l)
      j = sj(n)
      do i = sfi(n),sli(n)
         select case (bdy_2d_type(l))
            case (ZERO_GRADIENT)
               z(i,j) = z(i,j+1)
            case (SOMMERFELD)
!              KK-TODO: change DYC to DYV ?!
!                       change D(i,j) to _HALF_*(D(i,j)+D(i,j+1)) ?
               z(i,j) = z(i,j) + dtm/DYC*sqrt(9.81*D(i,j))*(z(i,j+1)-z(i,j))
            case (CLAMPED)
               z(i,j) = max(fac*bdy_data(k),-H(i,j)+min_depth)
            case (FLATHER_ELEV)
               a = sqrt(DV(i,j)/9.81)*(V(i,j)/DV(i,j)-bdy_data_v(k))
               z(i,j) = max(fac*(bdy_data(k) - a),-H(i,j)+min_depth)
            case default
               FATAL 'Illegal NSB 2D boundary type selection'
               stop 'update_2d_bdy()'
         end select
         k = k+1
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'leaving do_bdy_2d()'
   write(debug,*)
#endif
   return
   end subroutine do_bdy_2d
!EOC
!-----------------------------------------------------------------------

   end module bdy_2d

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Karsten Bolding and Hans Burchard               !
!-----------------------------------------------------------------------
