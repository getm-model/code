#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:  adv_upstream_2dh - 2DH upstream advection of 2D quantities
! \label{sec-upstream-2dh-adv}
!
! !INTERFACE:
   subroutine adv_upstream_2dh(dt,f,fi,Di,adv,U,V,Dn,DU,DV, &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                               dxv,dyu,dxu,dyv,arcd1,       &
#endif
                               action,AH,az)
!  Note (KK): Keep in sync with interface in advection.F90
!
! !DESCRIPTION:
! In this routine, the first-order upstream advection scheme
! is applied for the two horizontal directions in
! one step. The scheme should be positive definite and of
! high resolution. In order to remove truncation errors which might in Wadden
! Sea applications cause non-monotonicity, a truncation of over- and
! undershoots is carried out at the end of this subroutine. Such
! two-dimensional schemes are advantageous in Wadden Sea applications, since
! one-dimensional directional-split schemes might compute negative intermediate
! depths.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
#if !( defined(SPHERICAL) || defined(CURVILINEAR) )
   use domain, only: dx,dy,ard1
#endif
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                        :: dt,AH
   REALTYPE,dimension(E2DFIELD),intent(in)    :: f,U,V,Dn,DU,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
   REALTYPE,dimension(:,:),pointer,intent(in) :: dxu,dyu
   REALTYPE,dimension(_IRANGE_HALO_,_JRANGE_HALO_-1),intent(in) :: dxv,dyv
   REALTYPE,dimension(E2DFIELD),intent(in)    :: arcd1
#endif
   integer,intent(in)                         :: action
   integer,dimension(E2DFIELD),intent(in)     :: az
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(inout) :: fi,Di,adv
!
! !LOCAL VARIABLES:
   integer                      :: i,j,ii,jj
   REALTYPE                     :: Dio,advn
   REALTYPE,dimension(E2DFIELD) :: uflux,vflux
   REALTYPE,dimension(E2DFIELD) :: cmin,cmax
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'adv_upstream_2dh() # ',Ncall
#endif

#ifdef SLICE_MODEL
 stop 'adv_upstream_2dh(): Do not use adv_upstream_2dh in SLICE_MODEL mode'
#endif

! NOTE: With the present implementation it is not necessary
!   to initialize flx and fly. Even if they
!   are set to garbage here, then it does not change the result.
!      BJB 2009-09-25.


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ii,jj,Dio,advn)

!  Calculating u-interface low-order fluxes !
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin-1,imax
         if (U(i,j) .gt. _ZERO_) then
            uflux(i,j)=U(i,j)*f(i,j)
         else
            uflux(i,j)=U(i,j)*f(i+1,j)
         end if
         if ((AH.gt._ZERO_).and.(az(i,j).gt.0).and.(az(i+1,j).gt.0)) then
            uflux(i,j) = uflux(i,j) - AH*( f(i+1,j)                                 &
                                          -f(i  ,j))/DXU*_HALF_*(Dn(i+1,j)+Dn(i,j))
         end if
      end do
   end do
!$OMP END DO NOWAIT

#ifndef SLICE_MODEL
!  Calculating v-interface low-order fluxes !
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-1,jmax
      do i=imin,imax
         if (V(i,j) .gt. _ZERO_) then
            vflux(i,j)=V(i,j)*f(i,j)
         else
            vflux(i,j)=V(i,j)*f(i,j+1)
         end if
         if ((AH.gt._ZERO_).and.(az(i,j).gt.0).and.(az(i,j+1).gt.0)) then
            vflux(i,j) = vflux(i,j) - AH*( f(i,j+1)                                 &
                                          -f(i,j  ))/DYV*_HALF_*(Dn(i,j+1)+Dn(i,j))
         end if
      end do
   end do
!$OMP END DO NOWAIT
#endif

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1)  then
            cmin(i,j) =  100000*_ONE_
            cmax(i,j) = -100000*_ONE_
            do ii=-1,1
               if (az(i+ii,j).eq.1) then
                  if (f(i+ii,j).gt.cmax(i,j)) cmax(i,j)=f(i+ii,j)
                  if (f(i+ii,j).lt.cmin(i,j)) cmin(i,j)=f(i+ii,j)
               end if
            end do
            do jj=-1,1,2
               if (az(i,j+jj).eq.1) then
                  if (f(i,j+jj).gt.cmax(i,j)) cmax(i,j)=f(i,j+jj)
                  if (f(i,j+jj).lt.cmin(i,j)) cmin(i,j)=f(i,j+jj)
               end if
            end do
         end if
      end do
   end do
!$OMP END DO

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1)  then
            Dio = Di(i,j)
            Di(i,j) =  Dio - dt*(                                 &
                                   U(i,j)*DYU - U(i-1,j  )*DYUIM1 &
#ifndef SLICE_MODEL
                                 + V(i,j)*DXV - V(i  ,j-1)*DXVJM1 &
#endif
                                )*ARCD1
            advn = (                                     &
                      uflux(i,j)*DYU - uflux(i-1,j  )*DYUIM1 &
#ifndef SLICE_MODEL
                    + vflux(i,j)*DXV - vflux(i  ,j-1)*DXVJM1 &
#endif
                   )*ARCD1
            fi(i,j) = ( Dio*fi(i,j) - dt*advn ) / Di(i,j)
            adv(i,j) = adv(i,j) + advn
!           Force monotonicity, this is needed here for correcting truncations errors:
            if (fi(i,j).gt.cmax(i,j)) fi(i,j)=cmax(i,j)
            if (fi(i,j).lt.cmin(i,j)) fi(i,j)=cmin(i,j)
         end if
      end do
   end do
!$OMP END DO

!$OMP END PARALLEL

#ifdef DEBUG
   write(debug,*) 'Leaving adv_upstream_2dh()'
   write(debug,*)
#endif
   return
   end subroutine adv_upstream_2dh
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
