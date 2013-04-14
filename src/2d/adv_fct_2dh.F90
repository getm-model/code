#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:  adv_fct_2dh - 2DH FCT advection of 2D quantities \label{sec-fct-2dh-adv}
!
! !INTERFACE:
   subroutine adv_fct_2dh(fct,dt,f,fi,Di,adv,U,V,Dn,DU,DV, &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                          dxv,dyu,dxu,dyv,arcd1,           &
#endif
                          AH,az,                           &
                          mask_uflux,mask_vflux)
!  Note (KK): keep in sync with interface in advection.F90
!
! !DESCRIPTION:
! In this routine, the flux corrected transport advection scheme by
! \cite{ZALEZAK79} is applied for the two horizontal directions in
! one step. For details of this type of operator splitting, see  section
! \ref{sec-upstream-2dh-adv} on page \pageref{sec-upstream-2dh-adv}).
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
! solutions.
! Extra checks for boundaries including mirroring out of
! the transported quantities are performed in order to account for
! the third-order large stencils.
!
!  If GETM is executed as slice model (compiler option {\tt
!  SLICE\_MODEL}) the advection step for the $y$ direction is not
!  executed.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
#if !( defined(SPHERICAL) || defined(CURVILINEAR) )
   use domain, only: dx,dy,ard1
#endif
   use halo_zones, only : update_2d_halo,wait_halo,z_TAG
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical,intent(in)                         :: fct
   REALTYPE,intent(in)                        :: dt,AH
   REALTYPE,dimension(E2DFIELD),intent(in)    :: f,U,V,Dn,DU,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
   REALTYPE,dimension(:,:),pointer,intent(in) :: dxu,dyu
   REALTYPE,dimension(_IRANGE_HALO_,_JRANGE_HALO_-1),intent(in) :: dxv,dyv
   REALTYPE,dimension(E2DFIELD),intent(in)    :: arcd1
#endif
   integer,dimension(E2DFIELD),intent(in)     :: az
   logical,dimension(:,:),pointer,intent(in)  :: mask_uflux
   logical,dimension(_IRANGE_HALO_,_JRANGE_HALO_-1),intent(in) :: mask_vflux
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(inout) :: fi,Di,adv
!
! !LOCAL VARIABLES:
   integer         :: i,j
   REALTYPE,dimension(E2DFIELD) :: Dio
   REALTYPE,dimension(E2DFIELD) :: uflux,flx
#ifndef SLICE_MODEL
   REALTYPE,dimension(E2DFIELD) :: vflux,fly
#endif
   REALTYPE,dimension(E2DFIELD) :: faux,rp,rm,cmin,cmax
   REALTYPE        :: CNW,CW,CSW,CSSW,CWW,CSWW,CC,CS
   REALTYPE        :: advn,uuu,vvv,CExx,Cl,Cu,fac
   REALTYPE,parameter :: one12th=_ONE_/12,one6th=_ONE_/6,one3rd=_ONE_/3
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'adv_fct_2dh() # ',Ncall
#endif

! NOTE: With the present implementation it is not necessary
!   to initialize flx, fly, fhx, fhy, cmin and cmax. Even if they
!   are set to garbage here, then it does not change the result.
!      BJB 2009-09-25.

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP   PRIVATE(i,j,CNW,CW,CSW,CSSW,CWW,CSWW,CC,CS,advn,uuu,vvv,CExx,Cl,Cu,fac)

!  Calculate min and max of all values around one point
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO+1,jmax+HALO-1
      do i=imin-HALO+1,imax+HALO-1
         if (az(i,j) .eq. 1) then
            cmin(i,j) = minval(f(i-1:i+1,j-1:j+1),mask=(az(i-1:i+1,j-1:j+1).eq.1))
         end if
      end do
   end do
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO+1,jmax+HALO-1
      do i=imin-HALO+1,imax+HALO-1
         if (az(i,j) .eq. 1) then
            cmax(i,j) = maxval(f(i-1:i+1,j-1:j+1),mask=(az(i-1:i+1,j-1:j+1).eq.1))
         end if
      end do
   end do
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO+1,jmax+HALO-1
         do i=imin-HALO+1,imax+HALO-1
            if (az(i,j) .eq. 1)  then
               Dio(i,j) = Di(i,j)
               Di(i,j) = Dio(i,j) - dt*(                                &
                                         U(i  ,j)*DYU - U(i-1,j)*DYUIM1 &
#ifndef SLICE_MODEL
                                        +V(i,j  )*DXV - V(i,j-1)*DXVJM1 &
#endif
                                       )*ARCD1
            end if
         end do
      end do
!$OMP END DO NOWAIT

! OMP-TODO: The following loop gives errornous results when
!   threaded. - tried both at k,j, adn i. And I cant see what
!   should be wrong. Further, the original code needs an update
!   in any case, so for now I'll leave it in serial.
!   BJB 2009-09-23.
!$OMP MASTER
!  Calculating u-interface high-order fluxes !
   do j=jmin,jmax
      do i=imin-1,imax
         if (mask_uflux(i,j)) then
!           KK-TODO: calculation of uvel and vvel only once at start
            uuu=U(i,j)/DU(i,j)*dt/DXU
            vvv=_QUART_*(  V(i  ,j-1)/DV(i  ,j-1)         &
                         + V(i  ,j  )/DV(i  ,j  )         &
                         + V(i+1,j-1)/DV(i+1,j-1)         &
                         + V(i+1,j  )/DV(i+1,j  ))*dt/DYU
            if (uuu .gt. _ZERO_) then
               if (vvv .gt. _ZERO_) then
                  CNW =f(i  ,j+1)
                  CW  =f(i  ,j  )
                  CSW =f(i  ,j-1)
                  CSSW=f(i  ,j-2)
                  CWW =f(i-1,j  )
                  CSWW=f(i-1,j-1)
                  CC  =f(i+1,j  )
                  CS  =f(i+1,j-1)
                  if (az(i-1,j-1) .eq. 0) then
                     CSWW=CW
                  end if
                  if (az(i-1,j  ) .eq. 0) then
                     CWW=CW
                  end if
                  if (az(i  ,j+1) .eq. 0) then
                     CNW=CW
                  end if
                  if (az(i  ,j-1) .eq. 0) then
                     CSW=CW
                  end if
                  if (az(i  ,j-2) .eq. 0) then
                     CSSW=CNW
                  end if
                  if (az(i+1,j-1) .eq. 0) then
                     CS=CC
                  end if
               else
                  CNW =f(i  ,j-1)
                  CW  =f(i  ,j  )
                  CSW =f(i  ,j+1)
                  CSSW=f(i  ,j+2)
                  CWW =f(i-1,j  )
                  CSWW=f(i-1,j+1)
                  CC  =f(i+1,j  )
                  CS  =f(i+1,j+1)
                  if (az(i-1,j+1) .eq. 0) then
                     CSWW=CW
                  end if
                  if (az(i-1,j  ) .eq. 0) then
                     CWW=CW
                  end if
                  if (az(i  ,j-1) .eq. 0) then
                     CNW=CW
                  end if
                  if (az(i  ,j+1) .eq. 0) then
                     CSW=CW
                  end if
                  if (az(i  ,j+2) .eq. 0) then
                     CSSW=CNW
                  end if
                  if (az(i+1,j+1) .eq. 0) then
                     CS=CC
                  end if
               end if
            else
               if (vvv.gt._ZERO_) then
                  CNW =f(i+1,j+1)
                  CW  =f(i+1,j  )
                  CSW =f(i+1,j-1)
                  CSSW=f(i+1,j-2)
                  CWW =f(i+2,j  )
                  CSWW=f(i+2,j-1)
                  CC  =f(i  ,j  )
                  CS  =f(i  ,j-1)
                  if (az(i+2,j-1) .eq. 0) then
                     CSWW=CW
                  end if
                  if (az(i+2,j  ) .eq. 0) then
                     CWW=CW
                  end if
                  if (az(i+1,j+1) .eq. 0) then
                     CNW=CW
                  end if
                  if (az(i+1,j-1) .eq. 0) then
                     CSW=CW
                  end if
                  if (az(i+1,j-2) .eq. 0) then
                     CSSW=CNW
                  end if
                  if (az(i  ,j-1) .eq. 0) then
                     CS=CC
                  end if
               else
                  CNW =f(i+1,j-1)
                  CW  =f(i+1,j  )
                  CSW =f(i+1,j+1)
                  CSSW=f(i+1,j+2)
                  CWW =f(i+2,j  )
                  CSWW=f(i+2,j+1)
                  CC  =f(i  ,j  )
                  CS  =f(i  ,j+1)
                  if (az(i+2,j+1) .eq. 0) then
                     CSWW=CW
                  end if
                  if (az(i+2,j  ) .eq. 0) then
                     CWW=CW
                  end if
                  if (az(i+1,j-1) .eq. 0) then
                     CNW=CW
                  end if
                  if (az(i+1,j+1) .eq. 0) then
                     CSW=CW
                  end if
                  if (az(i+1,j+2) .eq. 0) then
                     CSSW=CNW
                  end if
                  if (az(i  ,j+1) .eq. 0) then
                     CS=CC
                  end if
               end if
            end if
            uuu=abs(uuu)
            vvv=abs(vvv)
            uflux(i,j) = (                                                    &
                            _HALF_                                            &
                            * ( CC + CW )                                     &
                          - _HALF_*uuu                                        &
                            * ( CC - CW )                                     &
                          - one6th*(_ONE_-uuu*uuu)                            &
                            * ( CC - _TWO_*CW + CWW )                         &
                          - _HALF_*vvv                                        &
                            * ( CW - CSW )                                    &
                          - vvv*(_QUART_-one3rd*uuu)                          &
                            * ( CC - CW - CS + CSW )                          &
                          - _HALF_*vvv*(_HALF_ -one3rd*vvv)                   &
                            * ( CNW - _TWO_*CW + CSW )                        &
                          + _QUART_*vvv*(one3rd-_HALF_*uuu*uuu)               &
                            * ( CC - _TWO_*CW + CWW - CS - _TWO_*CSW + CSWW ) &
                          + one12th*vvv*(_ONE_-_HALF_*vvv*vvv)                &
                            * ( CNW - _THREE_*CW + _THREE_*CSW - CSSW )       &
                         )                                                    &
                         *U(i,j)
         else
            uflux(i,j) = _ZERO_
         end if
      end do
   end do

#ifndef SLICE_MODEL
!  Calculating v-interface high-order fluxes !
   do j=jmin-1,jmax
      do i=imin,imax
         if (mask_vflux(i,j)) then
            uuu=V(i,j)/DV(i,j)*dt/DYV
            vvv=_QUART_*(  U(i-1,j  )/DU(i-1,j  )         &
                         + U(i-1,j+1)/DU(i-1,j+1)         &
                         + U(i  ,j  )/DU(i  ,j  )         &
                         + U(i  ,j+1)/DU(i  ,j+1))*dt/DXV
            if (uuu .gt. _ZERO_) then
               if (vvv .gt. _ZERO_) then
                  CNW =f(i+1,j  )
                  CW  =f(i  ,j  )
                  CSW =f(i-1,j  )
                  CSSW=f(i-2,j  )
                  CWW =f(i  ,j-1)
                  CSWW=f(i-1,j-1)
                  CC  =f(i  ,j+1)
                  CS  =f(i-1,j+1)
                  if (az(i-1,j-1) .eq. 0) then
                     CSWW=CW
                  end if
                  if (az(i  ,j-1) .eq. 0) then
                     CWW=CW
                  end if
                  if (az(i+1,j  ) .eq. 0) then
                     CNW=CW
                  end if
                  if (az(i-1,j  ) .eq. 0) then
                     CSW=CW
                  end if
                  if (az(i-2,j  ) .eq. 0) then
                     CSSW=CNW
                  end if
                  if (az(i-1,j+1) .eq. 0) then
                     CS=CC
                  end if
               else
                  CNW =f(i-1,j  )
                  CW  =f(i  ,j  )
                  CSW =f(i+1,j  )
                  CSSW=f(i+2,j  )
                  CWW =f(i  ,j-1)
                  CSWW=f(i+1,j-1)
                  CC  =f(i  ,j+1)
                  CS  =f(i+1,j+1)
                  if (az(i+1,j-1) .eq. 0) then
                     CSWW=CW
                  end if
                  if (az(i  ,j-1) .eq. 0) then
                     CWW=CW
                  end if
                  if (az(i-1,j  ) .eq. 0) then
                     CNW=CW
                  end if
                  if (az(i+1,j  ) .eq. 0) then
                     CSW=CW
                  end if
                  if (az(i+2,j  ) .eq. 0) then
                     CSSW=CNW
                  end if
                  if (az(i+1,j+1) .eq. 0) then
                     CS=CC
                  end if
               end if
            else
               if (vvv .gt. _ZERO_) then
                  CNW =f(i+1,j+1)
                  CW  =f(i  ,j+1)
                  CSW =f(i-1,j+1)
                  CSSW=f(i-2,j+1)
                  CWW =f(i  ,j+2)
                  CSWW=f(i-1,j+2)
                  CC  =f(i  ,j  )
                  CS  =f(i-1,j  )
                  if (az(i-1,j+2) .eq. 0) then
                     CSWW=CW
                  end if
                  if (az(i  ,j+2) .eq. 0) then
                     CWW=CW
                  end if
                  if (az(i+1,j+1) .eq. 0) then
                     CNW=CW
                  end if
                  if (az(i-1,j+1) .eq. 0) then
                     CSW=CW
                  end if
                  if (az(i-2,j+1) .eq. 0) then
                     CSSW=CNW
                  end if
                  if (az(i-1,j  ) .eq. 0) then
                     CS=CC
                  end if
               else
                  CNW =f(i-1,j+1)
                  CW  =f(i  ,j+1)
                  CSW =f(i+1,j+1)
                  CSSW=f(i+2,j+1)
                  CWW =f(i  ,j+2)
                  CSWW=f(i+1,j+2)
                  CC  =f(i  ,j  )
                  CS  =f(i+1,j  )
                  if (az(i+1,j+2) .eq. 0) then
                     CSWW=CW
                  end if
                  if (az(i  ,j+2) .eq. 0) then
                     CWW=CW
                  end if
                  if (az(i-1,j+1) .eq. 0) then
                     CNW=CW
                  end if
                  if (az(i+1,j+1) .eq. 0) then
                     CSW=CW
                  end if
                  if (az(i+2,j+1) .eq. 0) then
                     CSSW=CNW
                  end if
                  if (az(i+1,j  ) .eq. 0) then
                     CS=CC
                  end if
               end if
            end if
            uuu=abs(uuu)
            vvv=abs(vvv)
            vflux(i,j) = (                                                    &
                            _HALF_                                            &
                            * ( CC + CW )                                     &
                          - _HALF_*uuu                                        &
                            * ( CC - CW )                                     &
                          - one6th*(_ONE_-uuu*uuu)                            &
                            * ( CC - _TWO_*CW + CWW )                         &
                          - _HALF_*vvv                                        &
                            * ( CW - CSW )                                    &
                          - vvv*(_QUART_-one3rd*uuu)                          &
                            * ( CC - CW - CS + CSW )                          &
                          - _HALF_*vvv*(_HALF_-one3rd*vvv)                    &
                            * ( CNW - _TWO_*CW + CSW )                        &
                          + _QUART_*vvv*(one3rd-_HALF_*uuu*uuu)               &
                            * ( CC - _TWO_*CW + CWW - CS - _TWO_*CSW + CSWW ) &
                          + one12th*vvv*(_ONE_-_HALF_*vvv*vvv)                &
                            * ( CNW - _THREE_*CW + _THREE_*CSW - CSSW )       &
                         )                                                    &
                         *V(i,j)
         else
            vflux(i,j) = _ZERO_
         end if
      end do
   end do
#endif
!$OMP END MASTER


   if (fct) then

!     Calculate intermediate low resolution solution faux
!     Calculating u-interface low-order fluxes !
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO-1
            if (mask_uflux(i,j)) then
               if (U(i,j) .gt. _ZERO_) then
                  flx(i,j) = U(i,j)*f(i  ,j)
               else
                  flx(i,j) = U(i,j)*f(i+1,j)
               end if
            else
               flx(i,j) = _ZERO_
            end if
         end do
      end do
!$OMP END DO NOWAIT

#ifndef SLICE_MODEL
!     Calculating v-interface low-order fluxes !
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO,jmax+HALO-1
         do i=imin-HALO,imax+HALO
            if (mask_vflux(i,j)) then
               if (V(i,j) .gt. _ZERO_) then
                  fly(i,j) = V(i,j)*f(i,j  )
               else
                  fly(i,j) = V(i,j)*f(i,j+1)
               end if
            else
               fly(i,j) = _ZERO_
            end if
         end do
      end do
!$OMP END DO NOWAIT
#endif

!$OMP BARRIER
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO+1,jmax+HALO-1
         do i=imin-HALO+1,imax+HALO-1
            if (az(i,j) .eq. 1)  then
               faux(i,j) = ( f(i,j)*Dio(i,j) - dt*(                                  &
                                                  flx(i  ,j)*DYU-flx(i-1,j)*DYUIM1 &
#ifndef SLICE_MODEL
                                                 +fly(i,j  )*DXV-fly(i,j-1)*DXVJM1 &
#endif
                                                )*ARCD1 )/Di(i,j)
            end if
         end do
      end do
!$OMP END DO

!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin,imax
            if (az(i,j) .eq. 1) then
!              Calculate max and min of all values around one point
               cmax(i,j) = max( cmax(i,j) , maxval(faux(i-1:i+1,j-1:j+1),mask=(az(i-1:i+1,j-1:j+1).eq.1)) )
               cmin(i,j) = min( cmin(i,j) , minval(faux(i-1:i+1,j-1:j+1),mask=(az(i-1:i+1,j-1:j+1).eq.1)) )

!              max (Cu) possible concentration after a time step
               CExx=(                                                   &
                      ( min(uflux(i  ,j  )-flx(i  ,j  ),_ZERO_)*DYU     &
                       -max(uflux(i-1,j  )-flx(i-1,j  ),_ZERO_)*DYUIM1) &
#ifndef SLICE_MODEL
                     +( min(vflux(i  ,j  )-fly(i  ,j  ),_ZERO_)*DXV     &
                       -max(vflux(i  ,j-1)-fly(i  ,j-1),_ZERO_)*DXVJM1) &
#endif
                    )
               Cu=(faux(i,j)*Di(i,j)-dt*CExx*ARCD1)/Di(i,j)
!              calculating the maximum limiter rp for each conc. cell
               if (Cu .eq. faux(i,j)) then
                  rp(i,j)=_ZERO_
               else
                  rp(i,j)=min((cmax(i,j)-faux(i,j))/(Cu-faux(i,j)),_ONE_)
               end if

!              min (Cl) possible concentration after a time step
               CExx=(                                                   &
                      ( max(uflux(i  ,j  )-flx(i  ,j  ),_ZERO_)*DYU     &
                       -min(uflux(i-1,j  )-flx(i-1,j  ),_ZERO_)*DYUIM1) &
#ifndef SLICE_MODEL
                     +( max(vflux(i  ,j  )-fly(i  ,j  ),_ZERO_)*DXV     &
                       -min(vflux(i  ,j-1)-fly(i  ,j-1),_ZERO_)*DXVJM1) &
#endif
                    )
               Cl=(faux(i,j)*Di(i,j)-dt*CExx*ARCD1)/Di(i,j)
!              calculating the minimum limiter rm for each conc. cell
               if (Cl .eq. faux(i,j)) then
                  rm(i,j)=_ZERO_
               else
                  rm(i,j)=min((faux(i,j)-cmin(i,j))/(faux(i,j)-Cl),_ONE_)
               end if
            end if
         end do
      end do
!$OMP END DO

!     Note (KK): we need r[m|p] within [imin-1:imax+1,jmin-1:jmax+1]
!                however, only available within [imin:imax,jmin:jmax]
      call update_2d_halo(rm,rm,az,imin,jmin,imax,jmax,z_TAG)
      call wait_halo(z_TAG)
      call update_2d_halo(rp,rp,az,imin,jmin,imax,jmax,z_TAG)
      call wait_halo(z_TAG)

!     Limiters for the u-fluxes (fac)
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin-1,imax
            if (mask_uflux(i,j)) then
               if (uflux(i,j)-flx(i,j) .ge. _ZERO_) then
                  fac=min(rm(i,j),rp(i+1,j))
               else
                  fac=min(rm(i+1,j),rp(i,j))
               end if
               uflux(i,j) = (_ONE_-fac)*flx(i,j) + fac*uflux(i,j)
            end if
         end do
      end do
!$OMP END DO

#ifndef SLICE_MODEL
!     Limiters for the v-fluxes (fac)
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-1,jmax
         do i=imin,imax
            if (mask_vflux(i,j)) then
               if (vflux(i,j)-fly(i,j) .ge. _ZERO_) then
                  fac=min(rm(i,j),rp(i,j+1))
               else
                  fac=min(rm(i,j+1),rp(i,j))
               end if
               vflux(i,j) = (_ONE_-fac)*fly(i,j) + fac*vflux(i,j)
            end if
         end do
      end do
!$OMP END DO
#endif

   end if ! if (fct)

   if (AH .gt. _ZERO_) then
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin-1,imax
            if (mask_uflux(i,j)) then
               uflux(i,j) = uflux(i,j) - AH*( f(i+1,j)                                 &
                                             -f(i  ,j))/DXU*_HALF_*(Dn(i+1,j)+Dn(i,j))
            end if
         end do
      end do
!$OMP END DO NOWAIT

#ifndef SLICE_MODEL
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-1,jmax
         do i=imin,imax
            if (mask_vflux(i,j)) then
               vflux(i,j) = vflux(i,j) - AH*( f(i,j+1)                                 &
                                             -f(i,j  ))/DYV*_HALF_*(Dn(i,j+1)+Dn(i,j))
            end if
         end do
      end do
!$OMP END DO
#endif
   end if

! Doing the full advection in one step
!$OMP BARRIER
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1)  then
!           CAUTION: Di(i,j) already calculated above
            advn = (                                         &
                      uflux(i,j)*DYU - uflux(i-1,j  )*DYUIM1 &
#ifndef SLICE_MODEL
                    + vflux(i,j)*DXV - vflux(i  ,j-1)*DXVJM1 &
#endif
                   )*ARCD1
            fi(i,j) = ( Dio(i,j)*fi(i,j) - dt*advn ) / Di(i,j)
!           Force monotonicity, this is needed here for correcting truncations errors:
            if (fi(i,j).gt.cmax(i,j)) fi(i,j)=cmax(i,j)
            if (fi(i,j).lt.cmin(i,j)) fi(i,j)=cmin(i,j)
            adv(i,j) = adv(i,j) + advn
         end if
      end do
   end do
!$OMP END DO

!$OMP END PARALLEL

#ifdef DEBUG
   write(debug,*) 'Leaving adv_fct_2dh()'
   write(debug,*)
#endif
   return
   end subroutine adv_fct_2dh
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
