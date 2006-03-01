!$Id: fct_2dh_adv.F90,v 1.5 2006-03-01 16:03:32 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:  fct_2dh_adv - 2D flux-corrected transport \label{sec-fct-2dh-adv}
!
! !INTERFACE:
   subroutine fct_2dh_adv(dt,f,uu,vv,ho,hn,hun,hvn, &
                           delxv,delyu,delxu,delyv,area_inv,az,AH)
! !DESCRIPTION:
! In this routine, the flux corrected transport advection scheme by
! \cite{ZALEZAK79} is applied for the two horizontal directions in
! one step. For details of this type of operator splitting, see the
! routine {\tt upstream\_2dh\_adv} (see section \ref{sec-upstream-2dh-adv} 
! on page \pageref{sec-upstream-2dh-adv}). 
! 
! The monotone low-order flux is the first-order upstream 
! scheme, the
! high-order flux is the third-order ULTIMATE QUICKEST scheme by 
! \cite{LEONARDea95}. The scheme should thus be positive definite and of 
! high resolution. In order to remove truncation errors which might in Wadden
! Sea applications cause non-monotonicity, a truncation of over- and 
! undershoots is carried out at the end of this subroutine. Such 
! two-dimensional schemes are advantageous in Wadden Sea applications, since
! one-dimensional directioal-split schemes might compute negative intermediate
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
   use domain, only: iimin,iimax,jjmin,jjmax,kmax
   use advection_3d, only: hi,hio
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in) :: uu(I3DFIELD),vv(I3DFIELD)
   REALTYPE, intent(in) :: ho(I3DFIELD),hn(I3DFIELD)
   REALTYPE, intent(in) :: hun(I3DFIELD),hvn(I3DFIELD)
   REALTYPE, intent(in) :: delxv(I2DFIELD),delyu(I2DFIELD)
   REALTYPE, intent(in) :: delxu(I2DFIELD),delyv(I2DFIELD)
   REALTYPE, intent(in) :: area_inv(I2DFIELD),dt,AH
   integer, intent(in)  :: az(E2DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout)             :: f(I3DFIELD)
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer         :: rc,i,ii,j,jj,k,kk
#ifdef USE_ALLOCATED_ARRAYS
   REALTYPE, dimension(:,:,:), allocatable       :: flx,fly
   REALTYPE, dimension(:,:,:), allocatable       :: fhx,fhy
   REALTYPE, dimension(:,:,:), allocatable       :: fi
   REALTYPE, dimension(:,:,:), allocatable       :: rp,rm
   REALTYPE, dimension(:,:,:), allocatable       :: cmin,cmax
#else
   REALTYPE        :: flx(I3DFIELD),fly(I3DFIELD)
   REALTYPE        :: fhx(I3DFIELD),fhy(I3DFIELD)
   REALTYPE        :: fi(I3DFIELD)
   REALTYPE        :: rp(I3DFIELD),rm(I3DFIELD)
   REALTYPE        :: cmin(I3DFIELD),cmax(I3DFIELD)
#endif
   REALTYPE        :: CNW,CW,CSW,CSSW,CWW,CSWW,CC,CS
   REALTYPE        :: uuu,vvv,x,CExx,Cl,Cu,fac
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'fct_2dh_adv() # ',Ncall
#endif

#ifdef USE_ALLOCATED_ARRAYS
   allocate(flx(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'fct_2dh_adv: Error allocating memory (flx)'

   allocate(fly(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'fct_2dh_adv: Error allocating memory (fly)'

   allocate(fhx(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'fct_2dh_adv: Error allocating memory (fhx)'

   allocate(fhy(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'fct_2dh_adv: Error allocating memory (fhy)'

   allocate(fi(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'fct_2dh_adv: Error allocating memory (fi)'

   allocate(rp(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'fct_2dh_adv: Error allocating memory (rp)'

   allocate(rm(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'fct_2dh_adv: Error allocating memory (rm)'

   allocate(cmax(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'fct_2dh_adv: Error allocating memory (cmax)'

   allocate(cmin(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'fct_2dh_adv: Error allocating memory (cmin)'
#endif

   flx = _ZERO_
   do k=1,kmax   ! Calculating u-interface low-order fluxes !
      do j=jjmin,jjmax
         do i=iimin-1,iimax
            if (uu(i,j,k) .gt. _ZERO_) then
               flx(i,j,k)=uu(i,j,k)*f(i,j,k)
            else
               flx(i,j,k)=uu(i,j,k)*f(i+1,j,k)
            end if
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
         end do
      end do
   end do

   fhx = _ZERO_
   do k=1,kmax   ! Calculating u-interface high-order fluxes !
      do j=jjmin,jjmax
         do i=iimin-1,iimax
            uuu=uu(i,j,k)/hun(i,j,k)*dt/delxu(i,j)
            vvv=0.25*(vv(i  ,j-1,k)/hvn(i  ,j-1,k)           &
                     +vv(i  ,j  ,k)/hvn(i  ,j  ,k)           &
                     +vv(i+1,j-1,k)/hvn(i+1,j-1,k)           &
                     +vv(i+1,j  ,k)/hvn(i+1,j  ,k))*dt/delyu(i,j)
            if (uuu .gt. _ZERO_) then
               if (vvv .gt. _ZERO_) then
                  CNW =f(i  ,j+1,k) 
                  CW  =f(i  ,j  ,k)
                  CSW =f(i  ,j-1,k) 
                  CSSW=f(i  ,j-2,k)
                  CWW =f(i-1,j  ,k)
                  CSWW=f(i-1,j-1,k)
                  CC  =f(i+1,j  ,k)
                  CS  =f(i+1,j-1,k) 
               else 
                  CNW =f(i  ,j-1,k) 
                  CW  =f(i  ,j  ,k)
                  CSW =f(i  ,j+1,k) 
                  CSSW=f(i  ,j+2,k)
                  CWW =f(i-1,j  ,k)
                  CSWW=f(i-1,j+1,k)
                  CC  =f(i+1,j  ,k)
                  CS  =f(i+1,j+1,k) 
               end if  
            else
               if (vvv.gt._ZERO_) then
                  CNW =f(i+1,j+1,k) 
                  CW  =f(i+1,j  ,k)
                  CSW =f(i+1,j-1,k) 
                  CSSW=f(i+1,j-2,k)
                  CWW =f(i+2,j  ,k)
                  CSWW=f(i+2,j-1,k)
                  CC  =f(i  ,j  ,k)
                  CS  =f(i  ,j-1,k) 
               else 
                  CNW =f(i+1,j-1,k) 
                  CW  =f(i+1,j  ,k)
                  CSW =f(i+1,j+1,k) 
                  CSSW=f(i+1,j+2,k)
                  CWW =f(i+2,j  ,k)
                  CSWW=f(i+2,j+1,k)
                  CC  =f(i  ,j  ,k)
                  CS  =f(i  ,j+1,k) 
               end if  
            end if
            uuu=abs(uuu)
            vvv=abs(vvv)
            fhx(i,j,k)=                                                     &
            0.5*(CC+CW)-uuu/2.*(CC-CW)-(1.-uuu*uuu)/6.*(CC-2.*CW+CWW)
            fhx(i,j,k)=fhx(i,j,k)                                           &
            -vvv/2.*(CW-CSW)-vvv*(0.25-uuu/3.)*(CC-CW-CS+CSW)  
            fhx(i,j,k)=fhx(i,j,k)                                           &
            -vvv*(0.25-vvv/6.)*(CNW-2.*CW+CSW)  
            fhx(i,j,k)=fhx(i,j,k)                                           &
            +vvv*(1./12.-uuu*uuu/8.)*((CC-2.*CW+CWW)-(CS-2.*CSW+CSWW))  
            fhx(i,j,k)=fhx(i,j,k)                                           &
            +vvv*(1./12.-vvv*vvv/24.)*(CNW-3.*CW+3.*CSW-CSSW)  
            fhx(i,j,k)=fhx(i,j,k)*uu(i,j,k) 
         end do
      end do
   end do

   fhy = _ZERO_
   do k=1,kmax   ! Calculating v-interface high-order fluxes !
      do j=jjmin-1,jjmax
         do i=iimin,iimax
            uuu=vv(i,j,k)*dt/delyv(i,j) 
            vvv=0.25*(                             &
                    uu(i-1,j,k)/hun(i-1,j,k)       &
                   +uu(i-1,j+1,k)/hun(i-1,j+1,k)   &
                   +uu(i,j,k)/hun(i,j,k)           &
                   +uu(i,j+1,k)/hun(i,j+1,k)       &
                   )*dt/delxv(i,j)  
            if (uuu .gt. _ZERO_) then
               if (vvv .gt. _ZERO_) then
                  CNW =f(i+1,j  ,k)
                  CW  =f(i  ,j  ,k)
                  CSW =f(i-1,j  ,k)
                  CSSW=f(i-2,j  ,k)
                  CWW =f(i  ,j-1,k)
                  CSWW=f(i-1,j-1,k)
                  CC  =f(i  ,j+1,k)
                  CS  =f(i-1,j+1,k)
               else 
                  CNW =f(i-1,j  ,k)
                  CW  =f(i  ,j  ,k)
                  CSW =f(i+1,j  ,k)
                  CSSW=f(i+2,j  ,k)
                  CWW =f(i  ,j-1,k)
                  CSWW=f(i+1,j-1,k)
                  CC  =f(i  ,j+1,k)
                  CS  =f(i+1,j+1,k)
               end if  
            else
               if (vvv .gt. _ZERO_) then
                  CNW =f(i+1,j+1,k)
                  CW  =f(i  ,j+1,k)
                  CSW =f(i-1,j+1,k)
                  CSSW=f(i-2,j+1,k)
                  CWW =f(i  ,j+2,k)
                  CSWW=f(i-1,j+2,k)
                  CC  =f(i  ,j  ,k)
                  CS  =f(i-1,j  ,k)
               else 
                  CNW =f(i-1,j+1,k)
                  CW  =f(i  ,j+1,k)
                  CSW =f(i+1,j+1,k)
                  CSSW=f(i+2,j+1,k)
                  CWW =f(i  ,j+2,k)
                  CSWW=f(i+1,j+2,k)
                  CC  =f(i  ,j  ,k)
                  CS  =f(i+1,j  ,k)
               end if  
            end if
            uuu=abs(uuu)
            vvv=abs(vvv)
            fhy(i,j,k)=                                                     &
            0.5*(CC+CW)-uuu/2.*(CC-CW)-(1.-uuu*uuu)/6.*(CC-2.*CW+CWW)
            fhy(i,j,k)=fhy(i,j,k)                                           &
            -vvv/2.*(CW-CSW)-vvv*(0.25-uuu/3.)*(CC-CW-CS+CSW)  
            fhy(i,j,k)=fhy(i,j,k)                                           &
            -vvv*(0.25-vvv/6.)*(CNW-2.*CW+CSW)  
            fhy(i,j,k)=fhy(i,j,k)                                           &
            +vvv*(1./12.-uuu*uuu/8.)*((CC-2.*CW+CWW)-(CS-2.*CSW+CSWW))  
            fhy(i,j,k)=fhy(i,j,k)                                           &
            +vvv*(1./12.-vvv*vvv/24.)*(CNW-3.*CW+3.*CSW-CSSW)  
            fhy(i,j,k)=fhy(i,j,k)*vv(i,j,k) 
         end do
      end do
   end do

! Calculate intermediate low resolution solution fi 
   do k=1,kmax 
      do j=jjmin,jjmax
         do i=iimin,iimax
            if (az(i,j) .eq. 1)  then                                      
               hio(i,j,k)=hi(i,j,k)
               hi(i,j,k)=hio(i,j,k)                               &
               -dt*((uu(i  ,j,k)*delyu(i  ,j)-uu(i-1,j,k)*delyu(i-1,j)  &
                    +vv(i,j  ,k)*delxv(i,j  )-vv(i,j-1,k)*delxv(i,j-1)  &
                    )*area_inv(i,j))
               fi(i,j,k)=(f(i,j,k)*hio(i,j,k)                               &
               -dt*((flx(i  ,j,k)*delyu(i  ,j)-flx(i-1,j,k)*delyu(i-1,j)  &
                    +fly(i,j  ,k)*delxv(i,j  )-fly(i,j-1,k)*delxv(i,j-1)  &
                    )*area_inv(i,j)))/hi(i,j,k)
            end if 
         end do
      end do
   end do

! Calculating and applying the flux limiter
   do k=1,kmax  
      do j=jjmin,jjmax
         do i=iimin,iimax
            if (az(i,j) .eq. 1) then
               cmin(i,j,k)= 10000.
               cmax(i,j,k)=-10000.
! Calculate min and max of all values around one point
               do ii=-1,1 
                  do jj=-1,1
                     if (az(i+ii,j+jj).ge.1) then
                        x=min(f(i+ii,j+jj,k),fi(i+ii,j+jj,k))
                        if (x.lt.cmin(i,j,k)) cmin(i,j,k)=x
                        x=max(f(i+ii,j+jj,k),fi(i+ii,j+jj,k))
                        if (x .gt. cmax(i,j,k)) cmax(i,j,k)=x
                            end if
                  end do
               end do
! max (Cu) and min (Cl) possible concentration after a time step
               CExx=(min(fhx(i  ,j  ,k)-flx(i  ,j  ,k),_ZERO_) &
                    -max(fhx(i-1,j  ,k)-flx(i-1,j  ,k),_ZERO_))/delxu(i,j)     &
                   +(min(fhy(i  ,j  ,k)-fly(i  ,j  ,k),_ZERO_)                 &
                    -max(fhy(i  ,j-1,k)-fly(i  ,j-1,k),_ZERO_))/delyv(i,j)
               Cu=(fi(i,j,k)*hi(i,j,k)-dt*CExx)/hi(i,j,k)
               CExx=(max(fhx(i  ,j  ,k)-flx(i  ,j  ,k),_ZERO_)                 &
                    -min(fhx(i-1,j  ,k)-flx(i-1,j  ,k),_ZERO_))/delxu(i,j)     &
                   +(max(fhy(i  ,j  ,k)-fly(i  ,j  ,k),_ZERO_)                 &
                    -min(fhy(i  ,j-1,k)-fly(i  ,j-1,k),_ZERO_))/delyv(i,j)
               Cl=(fi(i,j,k)*hi(i,j,k)-dt*CExx)/hi(i,j,k)
! calculating the maximum limiters rp and rm for each conc. cell
               if (Cu .eq. fi(i,j,k)) then
                  rp(i,j,k)=_ZERO_
               else
                  rp(i,j,k)=min((cmax(i,j,k)-fi(i,j,k))/(Cu-fi(i,j,k)),_ONE_)
               end if
               if (Cl .eq. fi(i,j,k)) then
                  rm(i,j,k)=_ZERO_
               else
                  rm(i,j,k)=min((fi(i,j,k)-cmin(i,j,k))/(fi(i,j,k)-Cl),_ONE_)
               end if
            end if
         end do
      end do
   end do

!  Limiters for the u-fluxes (fac)
   do k=1,kmax   
      do j=jjmin,jjmax
         do i=iimin-1,iimax
            if (fhx(i,j,k)-flx(i,j,k).ge.0.) then
               fac=min(rm(i,j,k),rp(i+1,j,k))
            else
               fac=min(rm(i+1,j,k),rp(i,j,k))
            end if
            fhx(i,j,k)=(1.-fac)*flx(i,j,k)+fac*fhx(i,j,k)
            if ((AH .gt. _ZERO_) .and. (az(i,j) .gt. 0) .and. (az(i+1,j) .gt. 0)) &
               fhx(i,j,k)=fhx(i,j,k)-AH*(f(i+1,j,k)-f(i,j,k))/delxu(i,j) &
                         *0.5*(hn(i+1,j,k)+hn(i,j,k))
         end do
      end do
   end do

!  Limiters for the v-fluxes (fac)
   do k=1,kmax   
      do j=jjmin-1,jjmax
         do i=iimin,iimax
            if (fhy(i,j,k)-fly(i,j,k).ge.0.) then
               fac=min(rm(i,j,k),rp(i,j+1,k))
            else
               fac=min(rm(i,j+1,k),rp(i,j,k))
            end if
            fhy(i,j,k)=(1.-fac)*fly(i,j,k)+fac*fhy(i,j,k)
            if ((AH .gt. 0.) .and. (az(i,j) .gt. 0) .and. (az(i,j+1) .gt. 0))   &
               fhy(i,j,k)=fhy(i,j,k)-AH*(f(i,j+1,k)-f(i,j,k))/delyv(i,j)   &
                         *0.5*(hn(i,j+1,k)+hn(i,j,k))
         end do
      end do
   end do


! Doing the full advection in one step
   do k=1,kmax   
      do j=jjmin,jjmax
         do i=iimin,iimax
            if (az(i,j) .eq. 1)  then                                      
! CAUTION: hi(i,j,k) already calculated above
               f(i,j,k)=(f(i,j,k)*hio(i,j,k)                              &
               -dt*((fhx(i  ,j,k)*delyu(i  ,j)-fhx(i-1,j,k)*delyu(i-1,j)  &
                    +fhy(i,j  ,k)*delxv(i,j  )-fhy(i,j-1,k)*delxv(i,j-1)  &
                    )*area_inv(i,j)))/hi(i,j,k)
! Force monotonicity, this is needed here for correcting truncations errors:
               if (f(i,j,k) .gt. cmax(i,j,k)) f(i,j,k)=cmax(i,j,k)            
               if (f(i,j,k) .lt. cmin(i,j,k)) f(i,j,k)=cmin(i,j,k)            
            end if 
         end do
      end do
   end do

#ifdef USE_ALLOCATED_ARRAYS
#ifdef FORTRAN90
   deallocate(flx,stat=rc)    ! work array
   if (rc /= 0) stop 'fct_2dh_adv: Error de-allocating memory (flx)'
   deallocate(fly,stat=rc)    ! work array
   if (rc /= 0) stop 'fct_2dh_adv: Error de-allocating memory (fly)'
   deallocate(fhx,stat=rc)    ! work array
   if (rc /= 0) stop 'fct_2dh_adv: Error de-allocating memory (fhx)'
   deallocate(fhy,stat=rc)    ! work array
   if (rc /= 0) stop 'fct_2dh_adv: Error de-allocating memory (fhy)'
   deallocate(fi,stat=rc)    ! work array
   if (rc /= 0) stop 'fct_2dh_adv: Error de-allocating memory (fi)'
   deallocate(rp,stat=rc)    ! work array
   if (rc /= 0) stop 'fct_2dh_adv: Error de-allocating memory (rp)'
   deallocate(rm,stat=rc)    ! work array
   if (rc /= 0) stop 'fct_2dh_adv: Error de-allocating memory (rm)'
#endif
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving fct_2dh_adv()'
   write(debug,*)
#endif
   return
   end subroutine fct_2dh_adv
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
