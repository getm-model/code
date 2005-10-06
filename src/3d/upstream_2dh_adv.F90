!$Id: upstream_2dh_adv.F90,v 1.2 2005-10-06 09:54:01 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:  upstream_2dh_adv()
!
! !INTERFACE:
   subroutine upstream_2dh_adv(dt,f,uu,vv,ho,hn,hun,hvn, &
                           delxv,delyu,delxu,delyv,area_inv,az,AH)
!
! !DESCRIPTION:
! In this routine, the first-order upstream advection scheme 
! is applied for the two horizontal directions in
! one step. The scheme should be positive definite and of 
! high resolution. In order to remove truncation errors which might in Wadden
! Sea applications cause non-monotonicity, a truncation of over- and 
! undershoots is carried out at the end of this subroutine. Such 
! two-dimensional schemes are advantageous in Wadden Sea applications, since
! one-dimensional directioal-split schemes might compute negative intermediate
! depths.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
   use domain, only: iimin,iimax,jjmin,jjmax,kmax
   use advection_3d, only: hi,hio
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: uu(I3DFIELD),vv(I3DFIELD)
   REALTYPE, intent(in)                :: ho(I3DFIELD),hn(I3DFIELD)
   REALTYPE, intent(in)                :: hun(I3DFIELD),hvn(I3DFIELD)
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
!  $Log: upstream_2dh_adv.F90,v $
!  Revision 1.2  2005-10-06 09:54:01  hb
!  added support for vertical slice model - via -DSLICE_MODEL
!
!  Revision 1.1  2004/01/06 15:04:00  kbk
!  FCT advection + split of advection_3d.F90 + extra adv. input checks
!
!
! !LOCAL VARIABLES:
   integer         :: rc,i,j,k,ii,jj
#ifdef USE_ALLOCATED_ARRAYS
   REALTYPE, dimension(:,:,:), allocatable       :: flx,fly
   REALTYPE, dimension(:,:,:), allocatable       :: cmin,cmax
#else
   REALTYPE        :: flx(I3DFIELD),fly(I3DFIELD)
   REALTYPE        :: cmin(I3DFIELD),cmax(I3DFIELD)
#endif
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'upstream_2dh_adv() # ',Ncall
#endif

#ifdef USE_ALLOCATED_ARRAYS
   allocate(flx(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'upstream_2dh_adv: Error allocating memory (flx)'

   allocate(fly(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'upstream_2dh_adv: Error allocating memory (fly)'

   allocate(cmax(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'fct_2dh: Error allocating memory (cmax)'

   allocate(cmin(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'fct_2dh: Error allocating memory (cmin)'
#endif
!#ifdef SLICE_MODEL
 FATAL 'upstream_2dh_adv(): Do not use upstream_2dh_adv in SLICE_MODEL mode'
 stop
!#endif


   flx = _ZERO_
   do k=1,kmax   ! Calculating u-interface low-order fluxes !
      do j=jjmin,jjmax
         do i=iimin-1,iimax
            if (uu(i,j,k) .gt. _ZERO_) then
               flx(i,j,k)=uu(i,j,k)*f(i,j,k)
            else
               flx(i,j,k)=uu(i,j,k)*f(i+1,j,k)
            end if
            if ((AH.gt.0.).and.(az(i,j).gt.0).and.(az(i+1,j).gt.0))    &
               flx(i,j,k)=flx(i,j,k)-AH*(f(i+1,j,k)-f(i,j,k))/delxu(i,j) &
                         *0.5*(hn(i+1,j,k)+hn(i,j,k))
         end do
      end do
   end do

   fly = _ZERO_
   do k=1,kmax   ! Calculating v-interface low-order fluxes !
      do j=jjmin-1,jjmax
         do i=iimin,iimax
            if (vv(i,j,k) .gt. _ZERO_) then
               fly(i,j,k)=vv(i,j,k)*f(i,j,k)
            else
               fly(i,j,k)=vv(i,j,k)*f(i,j+1,k)
            end if
            if ((AH.gt.0.).and.(az(i,j).gt.0).and.(az(i,j+1).gt.0))   &
               fly(i,j,k)=fly(i,j,k)-AH*(f(i,j+1,k)-f(i,j,k))/delyv(i,j)   &
                         *0.5*(hn(i,j+1,k)+hn(i,j,k))
         end do
      end do
   end do

   cmax = -100000.
   cmin =  100000.
   do k=1,kmax 
      do j=jjmin,jjmax
         do i=iimin,iimax
            if (az(i,j) .eq. 1)  then                                      
               do ii=-1,1
                  if (az(i+ii,j).eq.1) then
                     if (f(i+ii,j,k).gt.cmax(i,j,k)) cmax(i,j,k)=f(i+ii,j,k)
                     if (f(i+ii,j,k).lt.cmin(i,j,k)) cmin(i,j,k)=f(i+ii,j,k)
                  end if 
               end do
               do jj=-1,1,2
                  if (az(i,j+jj).eq.1) then
                     if (f(i,j+jj,k).gt.cmax(i,j,k)) cmax(i,j,k)=f(i,j+jj,k)
                     if (f(i,j+jj,k).lt.cmin(i,j,k)) cmin(i,j,k)=f(i,j+jj,k)
                  end if 
               end do
            end if 
         end do
      end do
   end do

   do k=1,kmax 
      do j=jjmin,jjmax
         do i=iimin,iimax
            if (az(i,j) .eq. 1)  then                                      
               hio(i,j,k)=hi(i,j,k)
               hi(i,j,k)=hio(i,j,k)                               &
               -dt*((uu(i  ,j,k)*delyu(i  ,j)-uu(i-1,j,k)*delyu(i-1,j)  &
                    +vv(i,j  ,k)*delxv(i,j  )-vv(i,j-1,k)*delxv(i,j-1)  &
                    )*area_inv(i,j))
               f(i,j,k)=(f(i,j,k)*hio(i,j,k)                               &
               -dt*((flx(i  ,j,k)*delyu(i  ,j)-flx(i-1,j,k)*delyu(i-1,j)  &
                    +fly(i,j  ,k)*delxv(i,j  )-fly(i,j-1,k)*delxv(i,j-1)  &
                    )*area_inv(i,j)))/hi(i,j,k)
! Force monotonicity, this is needed here for correcting truncations errors:
               if (f(i,j,k).gt.cmax(i,j,k)) f(i,j,k)=cmax(i,j,k)            
               if (f(i,j,k).lt.cmin(i,j,k)) f(i,j,k)=cmin(i,j,k)            
            end if 
         end do
      end do
   end do

#ifdef USE_ALLOCATED_ARRAYS
#ifdef FORTRAN90
   deallocate(flx,stat=rc)    ! work array
   if (rc /= 0) stop 'upstream_2dh_adv: Error de-allocating memory (flx)'
   deallocate(fly,stat=rc)    ! work array
   if (rc /= 0) stop 'upstream_2dh_adv: Error de-allocating memory (fly)'
#endif
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving upstream_2dh_adv()'
   write(debug,*)
#endif
   return
   end subroutine upstream_2dh_adv
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
