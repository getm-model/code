!$Id: upstream_adv.F90,v 1.1 2004-01-06 15:04:00 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:  upstream_adv()
!
! !INTERFACE:
   subroutine upstream_adv(dt,f,uu,vv,ww,ho,hn, &
                           delxv,delyu,delxu,delyv,area_inv,az,AH)
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
   use domain, only: iimin,iimax,jjmin,jjmax,kmax
   use advection_3d, only: cu
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: uu(I3DFIELD),vv(I3DFIELD)
   REALTYPE, intent(in)                :: ww(I3DFIELD)
   REALTYPE, intent(in)                :: ho(I3DFIELD),hn(I3DFIELD)
   REALTYPE, intent(in)                :: delxv(I2DFIELD),delyu(I2DFIELD)
   REALTYPE, intent(in)                :: delxu(I2DFIELD),delyv(I2DFIELD)
   REALTYPE, intent(in)                :: area_inv(I2DFIELD),dt,AH
   integer, intent(in)                 :: az(E2DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout)             :: f(I3DFIELD)
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: upstream_adv.F90,v $
!  Revision 1.1  2004-01-06 15:04:00  kbk
!  FCT advection + split of advection_3d.F90 + extra adv. input checks
!
!
! !LOCAL VARIABLES:
   integer         :: rc,i,ii,j,jj,k,kk
#ifdef USE_ALLOCATED_ARRAYS
   REALTYPE, dimension(:,:,:), allocatable       :: adv
#else
   REALTYPE        :: adv(I3DFIELD)
#endif
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'upstream_adv() # ',Ncall
#endif

#ifdef USE_ALLOCATED_ARRAYS
   allocate(adv(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'upstream_adv: Error allocating memory (adv)'
#endif

   adv = _ZERO_

   cu = _ZERO_
   do k=1,kmax   ! Calculating u-interface fluxes !
      do j=jjmin,jjmax
         do i=iimin-1,iimax
            if (uu(i,j,k) .gt. _ZERO_) then
               cu(i,j,k)=uu(i,j,k)*f(i,j,k)
            else
               cu(i,j,k)=uu(i,j,k)*f(i+1,j,k)
            end if
            if ((AH.gt.0.).and.(az(i,j).gt.0).and.(az(i+1,j).gt.0))    &
               cu(i,j,k)=cu(i,j,k)-AH*(f(i+1,j,k)-f(i,j,k))/delxu(i,j) &
                         *0.5*(hn(i+1,j,k)+hn(i,j,k))
         end do
      end do
   end do
   do k=1,kmax   ! Updating the advection term for u-advection !
      do j=jjmin,jjmax
         do i=iimin,iimax
            adv(i,j,k)=(cu(i  ,j,k)*delyu(i  ,j)    &
                       -cu(i-1,j,k)*delyu(i-1,j))*area_inv(i,j)
         end do
      end do
   end do

   cu = _ZERO_
   do k=1,kmax   ! Calculating v-interface fluxes !
      do j=jjmin-1,jjmax
         do i=iimin,iimax
            if (vv(i,j,k) .gt. _ZERO_) then
               cu(i,j,k)=vv(i,j,k)*f(i,j,k)
            else
               cu(i,j,k)=vv(i,j,k)*f(i,j+1,k)
            end if
            if ((AH.gt.0.).and.(az(i,j).gt.0).and.(az(i,j+1).gt.0))   &
               cu(i,j,k)=cu(i,j,k)-AH*(f(i,j+1,k)-f(i,j,k))/delyv(i,j)   &
                         *0.5*(hn(i,j+1,k)+hn(i,j,k))
         end do
      end do
   end do
   do k=1,kmax   ! Updating the advection term for v-advection !
      do j=jjmin,jjmax
         do i=iimin,iimax
            adv(i,j,k)=adv(i,j,k)+(cu(i,j  ,k)*delxv(i,j  )   &
                                  -cu(i,j-1,k)*delxv(i,j-1))*area_inv(i,j)
         end do
      end do
   end do

   cu = _ZERO_
   if (kmax.gt.1) then
      do k=1,kmax-1   ! Calculating w-interface fluxes !
         do j=jjmin,jjmax
            do i=iimin,iimax
               if (ww(i,j,k) .gt. _ZERO_) then
                  cu(i,j,k)=ww(i,j,k)*f(i,j,k)
               else
                  cu(i,j,k)=ww(i,j,k)*f(i,j,k+1)
               end if
            end do
         end do
      end do
      do k=1,kmax   ! Updating the advection term for w-advection !
         do j=jjmin,jjmax
            do i=iimin,iimax
               adv(i,j,k)=adv(i,j,k)+(cu(i,j,k)-cu(i,j,k-1))
            end do
         end do
      end do
   end if

   do k=1,kmax   ! Doing the full advection in one step
      do j=jjmin,jjmax
         do i=iimin,iimax
            if (az(i,j) .eq. 1)                                        &
               f(i,j,k)=(f(i,j,k)*ho(i,j,k)-dt*adv(i,j,k))/hn(i,j,k)
         end do
      end do
   end do

#ifdef USE_ALLOCATED_ARRAYS
#ifdef FORTRAN90
   deallocate(adv,stat=rc)    ! work array
   if (rc /= 0) stop 'upstream_adv: Error de-allocating memory (adv)'
#endif
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving upstream_adv()'
   write(debug,*)
#endif
   return
   end subroutine upstream_adv
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
