#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: tow() - calculates cell centered physical vertical velocity
!
! !INTERFACE:
   subroutine tow(imin,jmin,imax,jmax,kmin,kmax,az,                    &
                  dt,                                                  &
#if defined(CURVILINEAR) || defined(SPHERICAL)
                  dxv,dyu,arcd1,                                       &
#else
                  dx,dy,ard1,                                          &
#endif
                  H,HU,HV,hn,ho,uu,hun,vv,hvn,ww,missing,wc)
!
! !DESCRIPTION:
!
! Here the physical vertical velocity at cell centres (but velocity time
! stages) is calculated.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)                       :: imin,jmin,imax,jmax,kmax
   integer,dimension(I2DFIELD),intent(in)   :: kmin
   integer,dimension(E2DFIELD),intent(in)   :: az
   REALTYPE,intent(in)                      :: dt
#if defined(CURVILINEAR) || defined(SPHERICAL)
   REALTYPE,dimension(E2DFIELD),intent(in)  :: dxv,dyu,arcd1
#else
   REALTYPE,intent(in)                      :: dx,dy,ard1
#endif
   REALTYPE,dimension(E2DFIELD),intent(in)  :: H,HU,HV
   REALTYPE,dimension(I3DFIELD),intent(inout) :: hn,ho,hun,hvn
   REALTYPE,dimension(I3DFIELD),intent(in)  :: uu,vv,ww
   REALTYPE,intent(in)                      :: missing
!
! !OUTPUT PARAMETERS:
   REALTYPE,dimension(I3DFIELD),intent(out) :: wc
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
   REALTYPE,dimension(I2DFIELD)             :: hwc,zco,zcn,zu
#ifndef SLICE_MODEL
   REALTYPE,dimension(I2DFIELD)             :: zv
#endif
   REALTYPE,dimension(I2DFIELD,0:1)         :: zw
   REALTYPE                                 :: dtm1
   integer                                  :: i,j,k,kp,km
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'tow() # ',Ncall
#endif
#ifdef SLICE_MODEL
   Note (KK): this value MUST NOT be changed !!!
   j = jmax/2
#endif

!  KK-TODO: this should be part of do_coordinates
!           ( 0 <= k < k[ |u|v]min-1 )
   ho(:,:,0) = _ZERO_
   hn(:,:,0) = _ZERO_
   hun(:,:,0) = _ZERO_
   hvn(:,:,0) = _ZERO_

   dtm1 = _ONE_/dt

!  initialise
   zco = -H
   zcn = -H
   zu = -HU
#ifndef SLICE_MODEL
   zv = -HV
#endif
   zw(:,:,0) = -H
   kp = 0

   do k=1,kmax

      hwc = _HALF_*( ho(:,:,k) + hn(:,:,k) )

!     update z-levels
      zco = zco + _HALF_*(  ho(:,:,k-1) +  ho(:,:,k) )
      zcn = zcn + _HALF_*(  hn(:,:,k-1) +  hn(:,:,k) )
      zu  = zu  + _HALF_*( hun(:,:,k-1) + hun(:,:,k) )
#ifndef SLICE_MODEL
      zv  = zv  + _HALF_*( hvn(:,:,k-1) + hvn(:,:,k) )
#endif
      km = kp
      kp = 1-kp
      zw(:,:,kp) = zw(:,:,km) + hwc

!     calculate wc
#ifndef SLICE_MODEL
      do j=jmin,jmax
#endif
         do i=imin,imax
            if (az(i,j) .eq. 1) then
               wc(i,j,k) =                               &
                  (                                      &
                     (                                   &
                        hn(i,j,k)*zcn(i,j)               &
                      - ho(i,j,k)*zco(i,j)               &
                     )                                   &
                     *dtm1                               &
                   + (                                   &
                        uu(i  ,j  ,k)*zu(i  ,j  )*DYU    &
                      - uu(i-1,j  ,k)*zu(i-1,j  )*DYUIM1 &
#ifndef SLICE_MODEL
                      + vv(i  ,j  ,k)*zv(i  ,j  )*DXV    &
                      - vv(i  ,j-1,k)*zv(i  ,j-1)*DXVJM1 &
#endif
                     )                                   &
                     *ARCD1                              &
                   + (                                   &
                        ww(i,j,k  )*zw(i,j,kp)           &
                      - ww(i,j,k-1)*zw(i,j,km)           &
                     )                                   &
                  )                                      &
                  /hwc(i,j)
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif

   end do

#ifndef SLICE_MODEL
   do j=jmin,jmax
#endif
      do i=imin,imax
         if (az(i,j) .eq. 2) then
            wc(i,j,:) = _ZERO_
         end if
         if (az(i,j) .eq. 0) then
!           Note (KK): can be skipped if wc is initialised with missing
            wc(i,j,:) = missing
         else
            do k=0,kmin(i,j)-1
               wc(i,j,k) = missing
            end do
         end if
      end do
#ifndef SLICE_MODEL
   end do
#endif

#ifdef SLICE_MODEL
   wc(:,j+1,:) = wc(:,j,:)
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving tow()'
   write(debug,*)
#endif
   return
   end subroutine tow
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2012 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
