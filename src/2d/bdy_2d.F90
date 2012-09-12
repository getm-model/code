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
   use parameters,only: g
   use halo_zones, only : z_TAG,H_TAG,U_TAG,V_TAG
   use domain, only: imin,jmin,imax,jmax,kmax,H,az,au,av
   use domain, only: nsbv,NWB,NNB,NEB,NSB,bdy_index,bdy_2d_type
   use domain, only: wi,wfj,wlj,nj,nfi,nli,ei,efj,elj,sj,sfi,sli
   use domain, only: min_depth
   use variables_2d, only: dtm,z,zo,D,U,DU,V,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxu,dyv
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
   REALTYPE                  :: cfl,depth,a,fac
   integer                   :: i,j,k,l,n
   REALTYPE, parameter       :: theta = _HALF_
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

#if 0
   select case (bdyfmt_2d)
      case (NO_DATA)
      case (ANALYTICAL)
      case (ASCII)
      case (NETCDF)
!        Read in get_2d_bdy() via get_2d_bdy_ncdf()
      case default
         stop 'do_bdy_2d(): invalid bdyfmt_2d'
   end select
#endif

!  Data read - do time interpolation

   fac = _ONE_
   if(bdyramp_2d .gt. 1) fac=min( _ONE_ ,FOUR*loop/float(bdyramp_2d))


   select case (tag)

      case (z_TAG,H_TAG)

         l = 0
         do n = 1,NWB
            l = l+1
            k = bdy_index(l)
            i = wi(n)
            select case (bdy_2d_type(l))
               case (ZERO_GRADIENT,CLAMPED_VEL,FLATHER_VEL)
                  do j = wfj(n),wlj(n)
                     z(i,j) = z(i+1,j)
                  end do
               case (SOMMERFELD)
                  do j = wfj(n),wlj(n)
                     cfl = sqrt(g*_HALF_*(D(i,j)+D(i+1,j)))*dtm/DXU
                     z(i,j) = (                                             &
                                (_ONE_ - _TWO_*cfl*(_ONE_-theta))*z (i  ,j) &
                               +(_ONE_ + _TWO_*cfl*(_ONE_-theta))*zo(i+1,j) &
                               -(_ONE_ - _TWO_*cfl*theta        )*z (i+1,j) &
                              )/(_ONE_ + _TWO_*cfl*theta        )
                  end do
               case (CLAMPED_ELEV,CLAMPED)
                  do j = wfj(n),wlj(n)
                     z(i,j) = max(fac*bdy_data(k),-H(i,j)+min_depth)
                     k = k+1
                  end do
               case (FLATHER_ELEV)
                  do j = wfj(n),wlj(n)
!                    Note (KK): approximate interface depths at vel-time stage
!                               by spatial mean at last sse-time stage
                     depth = _HALF_*(D(i,j)+D(i+1,j))
!                    Note (KK): note approximation of sse at vel-time stage
                     a = fac*bdy_data(k) &
                         - _TWO_/sqrt(g*depth)*(U(i,j)-fac*bdy_data_u(k)*depth)
                     z(i,j) = max(a,-H(i,j)+min_depth)
                     k = k+1
                  end do
            end select
         end do
         do n = 1,NNB
            l = l+1
            k = bdy_index(l)
            j = nj(n)
            select case (bdy_2d_type(l))
               case (ZERO_GRADIENT,CLAMPED_VEL,FLATHER_VEL)
                  do i = nfi(n),nli(n)
                     z(i,j) = z(i,j-1)
                  end do
               case (SOMMERFELD)
                  do i = nfi(n),nli(n)
                     cfl = sqrt(g*_HALF_*(D(i,j-1)+D(i,j)))*dtm/DYVJM1
                     z(i,j) = (                                             &
                                (_ONE_ - _TWO_*cfl*(_ONE_-theta))*z (i,j  ) &
                               +(_ONE_ + _TWO_*cfl*(_ONE_-theta))*zo(i,j-1) &
                               -(_ONE_ - _TWO_*cfl*theta        )*z (i,j-1) &
                              )/(_ONE_ + _TWO_*cfl*theta        )
                  end do
               case (CLAMPED_ELEV,CLAMPED)
                  do i = nfi(n),nli(n)
                     z(i,j) = max(fac*bdy_data(k),-H(i,j)+min_depth)
                     k = k+1
                  end do
               case (FLATHER_ELEV)
                  do i = nfi(n),nli(n)
!                    Note (KK): approximate interface depths at vel-time stage
!                               by spatial mean at last sse-time stage
                     depth = _HALF_*(D(i,j-1)+D(i,j))
!                    Note (KK): note approximation of sse at vel-time stage
                     a = fac*bdy_data(k) &
                         + _TWO_/sqrt(g*depth)*(V(i,j-1)-fac*bdy_data_v(k)*depth)
                     z(i,j) = max(a,-H(i,j)+min_depth)
                     k = k+1
                  end do
            end select
         end do
         do n = 1,NEB
            l = l+1
            k = bdy_index(l)
            i = ei(n)
            select case (bdy_2d_type(l))
               case (ZERO_GRADIENT,CLAMPED_VEL,FLATHER_VEL)
                  do j = efj(n),elj(n)
                     z(i,j) = z(i-1,j)
                  end do
               case (SOMMERFELD)
                  do j = efj(n),elj(n)
                     cfl = sqrt(g*_HALF_*(D(i-1,j)+D(i,j)))*dtm/DXUIM1
                     z(i,j) = (                                             &
                                (_ONE_ - _TWO_*cfl*(_ONE_-theta))*z (i  ,j) &
                               +(_ONE_ + _TWO_*cfl*(_ONE_-theta))*zo(i-1,j) &
                               -(_ONE_ - _TWO_*cfl*theta        )*z (i-1,j) &
                              )/(_ONE_ + _TWO_*cfl*theta        )
                  end do
               case (CLAMPED_ELEV,CLAMPED)
                  do j = efj(n),elj(n)
                     z(i,j) = max(fac*bdy_data(k),-H(i,j)+min_depth)
                     k = k+1
                  end do
               case (FLATHER_ELEV)
                  do j = efj(n),elj(n)
!                    Note (KK): approximate interface depths at vel-time stage
!                               by spatial mean at last sse-time stage
                     depth = _HALF_*(D(i-1,j)+D(i,j))
!                    Note (KK): note approximation of sse at vel-time stage
                     a = fac*bdy_data(k) &
                         + _TWO_/sqrt(g*depth)*(U(i-1,j)-fac*bdy_data_u(k)*depth)
                     z(i,j) = max(a,-H(i,j)+min_depth)
                     k = k+1
                  end do
            end select
         end do
         do n = 1,NSB
            l = l+1
            k = bdy_index(l)
            j = sj(n)
            select case (bdy_2d_type(l))
               case (ZERO_GRADIENT,CLAMPED_VEL,FLATHER_VEL)
                  do i = sfi(n),sli(n)
                     z(i,j) = z(i,j+1)
                  end do
               case (SOMMERFELD)
                  do i = sfi(n),sli(n)
                     cfl = sqrt(g*_HALF_*(D(i,j)+D(i,j+1)))*dtm/DYV
                     z(i,j) = (                                             &
                                (_ONE_ - _TWO_*cfl*(_ONE_-theta))*z (i,j  ) &
                               +(_ONE_ + _TWO_*cfl*(_ONE_-theta))*zo(i,j+1) &
                               -(_ONE_ - _TWO_*cfl*theta        )*z (i,j+1) &
                              )/(_ONE_ + _TWO_*cfl*theta        )
                  end do
               case (CLAMPED_ELEV,CLAMPED)
                  do i = sfi(n),sli(n)
                     z(i,j) = max(fac*bdy_data(k),-H(i,j)+min_depth)
                     k = k+1
                  end do
               case (FLATHER_ELEV)
                  do i = sfi(n),sli(n)
!                    Note (KK): approximate interface depths at vel-time stage
!                               by spatial mean at last sse-time stage
                     depth = _HALF_*(D(i,j)+D(i,j+1))
!                    Note (KK): note approximation of sse at vel-time stage
                     a = fac*bdy_data(k) &
                         - _TWO_/sqrt(g*depth)*(V(i,j)-fac*bdy_data_v(k)*depth)
                     z(i,j) = max(a,-H(i,j)+min_depth)
                     k = k+1
                  end do
            end select
         end do

      case (U_TAG)

         l = 0
         do n = 1,NWB
            l = l+1
            k = bdy_index(l)
            i = wi(n)
            select case (bdy_2d_type(l))
               case (FLATHER_VEL)
                  do j = wfj(n),wlj(n)
!                    Note (KK): approximate interface depths at vel-time stage
!                               by spatial mean at last sse-time stage
                     depth = _HALF_*(D(i,j)+D(i+1,j))
!                    Note (KK): note approximation of sse at vel-time stage
                     U(i,j) = fac*bdy_data_u(k)*depth &
                              - _HALF_*sqrt(g*depth)*(z(i,j)-fac*bdy_data(k))
                     k = k+1
                  end do
               case (CLAMPED_VEL,CLAMPED)
                  do j = wfj(n),wlj(n)
!                    Note (KK): approximate interface depths at vel-time stage
!                               by spatial mean at last sse-time stage
                     depth = _HALF_*(D(i,j)+D(i+1,j))
                     U(i,j) = fac*bdy_data_u(k)*depth
                     k = k+1
                  end do
            end select
         end do
         l = l + NNB
         do n = 1,NEB
            l = l+1
            k = bdy_index(l)
            i = ei(n)
            select case (bdy_2d_type(l))
               case (FLATHER_VEL)
                  do j = efj(n),elj(n)
!                    Note (KK): approximate interface depths at vel-time stage
!                               by spatial mean at last sse-time stage
                     depth = _HALF_*(D(i-1,j)+D(i,j))
!                    Note (KK): note approximation of sse at vel-time stage
                     U(i-1,j) = fac*bdy_data_u(k)*depth &
                                + _HALF_*sqrt(g*depth)*(z(i,j)-fac*bdy_data(k))
                     k = k+1
                  end do
               case (CLAMPED_VEL,CLAMPED)
                  do j = efj(n),elj(n)
!                    Note (KK): approximate interface depths at vel-time stage
!                               by spatial mean at last sse-time stage
                     depth = _HALF_*(D(i-1,j)+D(i,j))
                     U(i-1,j) = fac*bdy_data_u(k)*depth
                     k = k+1
                  end do
            end select
         end do

      case (V_TAG)

         l = NWB
         do n = 1,NNB
            l = l+1
            k = bdy_index(l)
            j = nj(n)
            select case (bdy_2d_type(l))
               case (FLATHER_VEL)
                  do i = nfi(n),nli(n)
!                    Note (KK): approximate interface depths at vel-time stage
!                               by spatial mean at last sse-time stage
                     depth = _HALF_*(D(i,j-1)+D(i,j))
!                    Note (KK): note approximation of sse at vel-time stage
                     V(i,j-1) = fac*bdy_data_v(k)*depth &
                                + _HALF_*sqrt(g*depth)*(z(i,j)-fac*bdy_data(k))
                     k = k+1
                  end do
               case (CLAMPED_VEL,CLAMPED)
                  do i = nfi(n),nli(n)
!                    Note (KK): approximate interface depths at vel-time stage
!                               by spatial mean at last sse-time stage
                     depth = _HALF_*(D(i,j-1)+D(i,j))
                     V(i,j-1) = fac*bdy_data_v(k)*depth
                     k = k+1
                  end do
            end select
         end do
         l = l + NEB
         do n = 1,NSB
            l = l+1
            k = bdy_index(l)
            j = sj(n)
            select case (bdy_2d_type(l))
               case (FLATHER_VEL)
                  do i = sfi(n),sli(n)
!                    Note (KK): approximate interface depths at vel-time stage
!                               by spatial mean at last sse-time stage
                     depth = _HALF_*(D(i,j)+D(i,j+1))
!                    Note (KK): note approximation of sse at vel-time stage
                     V(i,j) = fac*bdy_data_v(k)*depth &
                              - _HALF_*sqrt(g*depth)*(z(i,j)-fac*bdy_data(k))
                     k = k+1
                  end do
               case (CLAMPED_VEL,CLAMPED)
                  do i = sfi(n),sli(n)
!                    Note (KK): approximate interface depths at vel-time stage
!                               by spatial mean at last sse-time stage
                     depth = _HALF_*(D(i,j)+D(i,j+1))
                     V(i,j) = fac*bdy_data_v(k)*depth
                     k = k+1
                  end do
            end select
         end do

   end select

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
