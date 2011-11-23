#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:  adv_fct_2dh - 2D flux-corrected transport \label{sec-fct-2dh-adv}
!
! !INTERFACE:
   subroutine adv_fct_2dh(dt,f,Di,adv,U,V,Do,Dn,DU,DV, &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                          dxv,dyu,dxu,dyv,arcd1,       &
#endif
                          az,AH,nosplit_finalise)
!  Note (KK): keep in sync with interface in advection.F90
!
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
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                        :: dt,AH
   REALTYPE,dimension(E2DFIELD),intent(in)    :: U,V,Do,Dn,DU,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
   REALTYPE,dimension(:,:),pointer,intent(in) :: dxu,dyu
   REALTYPE,dimension(_IRANGE_HALO_,_JRANGE_HALO_-1),intent(in) :: dxv,dyv
   REALTYPE,dimension(E2DFIELD),intent(in)    :: arcd1
#endif
   integer,dimension(E2DFIELD),intent(in)     :: az
   logical,intent(in),optional                :: nosplit_finalise
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(inout) :: f,Di,adv
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   logical,save    :: first=.true.
   logical         :: update_f
   integer         :: rc,i,ii,j,jj
   REALTYPE,dimension(:,:),allocatable,save    :: Dio
#ifdef USE_ALLOCATED_ARRAYS
   REALTYPE, dimension(:,:), allocatable       :: flx,fhx
#ifndef SLICE_MODEL
   REALTYPE, dimension(:,:), allocatable       :: fly,fhy
#endif
   REALTYPE, dimension(:,:), allocatable       :: fi
   REALTYPE, dimension(:,:), allocatable       :: rp,rm
   REALTYPE, dimension(:,:), allocatable       :: cmin,cmax
#else
   REALTYPE        :: flx(E2DFIELD),fhx(E2DFIELD)
#ifndef SLICE_MODEL
   REALTYPE        :: fly(E2DFIELD),fhy(E2DFIELD)
#endif
   REALTYPE        :: fi(E2DFIELD)
   REALTYPE        :: rp(E2DFIELD),rm(E2DFIELD)
   REALTYPE        :: cmin(E2DFIELD),cmax(E2DFIELD)
#endif
   REALTYPE        :: CNW,CW,CSW,CSSW,CWW,CSWW,CC,CS
   REALTYPE        :: advn,uuu,vvv,x,CExx,Cl,Cu,fac
   REALTYPE,parameter :: one12th=_ONE_/12,one6th=_ONE_/6,one3rd=_ONE_/3
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'adv_fct_2dh() # ',Ncall
#endif

   stop 'adv_fct_2dh: This routine is buggy (KK)'
   if (first) then
      allocate(Dio(E2DFIELD),stat=rc)    ! work array
      if (rc /= 0) stop 'adv_fct_2dh: Error allocating memory (Dio)'
      Dio = _ZERO_

      first = .false.
   end if

   update_f = .true.
   if (present(nosplit_finalise)) then
      if (.not. nosplit_finalise) update_f = .false.
   end if

#ifdef USE_ALLOCATED_ARRAYS
   allocate(flx(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'adv_fct_2dh: Error allocating memory (flx)'

   allocate(fhx(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'adv_fct_2dh: Error allocating memory (fhx)'

#ifndef SLICE_MODEL
   allocate(fly(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'adv_fct_2dh: Error allocating memory (fly)'

   allocate(fhy(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'adv_fct_2dh: Error allocating memory (fhy)'
#endif

   allocate(fi(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'adv_fct_2dh: Error allocating memory (fi)'

   allocate(rp(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'adv_fct_2dh: Error allocating memory (rp)'

   allocate(rm(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'adv_fct_2dh: Error allocating memory (rm)'

   allocate(cmax(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'adv_fct_2dh: Error allocating memory (cmax)'

   allocate(cmin(E2DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'adv_fct_2dh: Error allocating memory (cmin)'
#endif

! NOTE: With the present implementation it is not necessary
!   to initialize flx, fly, fhx, fhy, cmin and cmax. Even if they
!   are set to garbage here, then it does not change the result.
!      BJB 2009-09-25.

! BJB-TODO/BUG: There is an error in the use of fi, rp and rm.
!  They are used in a larger region than they are computed.
!  As a result, the result of the rpesent scheme changes if
!  they are initialized. Presently, tehre is no initialization
!  of these arrays even though perhaps there should be.
!      BJB 2009-09-25.

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP   PRIVATE(i,j,CNW,CW,CSW,CSSW,CWW,CSWW,CC,CS,advn,uuu,vvv,x,CExx,Cl,Cu,fac)

!  Calculating u-interface low-order fluxes !
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin-1,imax
         if (U(i,j) .gt. _ZERO_) then
            flx(i,j)=U(i,j)*f(i  ,j)
         else
            flx(i,j)=U(i,j)*f(i+1,j)
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
            fly(i,j)=V(i,j)*f(i,j  )
         else
            fly(i,j)=V(i,j)*f(i,j+1)
         end if
      end do
   end do
!$OMP END DO NOWAIT
#endif

! OMP-TODO: The following loop gives errornous results when
!   threaded. - tried both at k,j, adn i. And I cant see what
!   should be wrong. Further, the original code needs an update
!   in any case, so for now I'll leave it in serial.
!   BJB 2009-09-23.
!$OMP MASTER
!  Calculating u-interface high-order fluxes !
   do j=jmin,jmax
      do i=imin-1,imax
!        KK-TODO: calculation of uvel and vvel only once at start
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
         fhx(i,j) = (                                                    &
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
      end do
   end do

#ifndef SLICE_MODEL
!  Calculating v-interface high-order fluxes !
   do j=jmin-1,jmax
      do i=imin,imax
!        Note (KK): added division by DV
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
         fhy(i,j) = (                                                    &
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
      end do
   end do
#endif
!$OMP END MASTER
!$OMP BARRIER

!  Calculate intermediate low resolution solution fi
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1)  then
            Dio(i,j) = Di(i,j)
            Di(i,j) = Dio(i,j) - dt*(                                &
                                      U(i  ,j)*DYU - U(i-1,j)*DYUIM1 &
#ifndef SLICE_MODEL
                                     +V(i,j  )*DXV - V(i,j-1)*DXVJM1 &
#endif
                                    )*ARCD1
            fi(i,j) = ( f(i,j)*Dio(i,j) - dt*(                                  &
                                               flx(i  ,j)*DYU-flx(i-1,j)*DYUIM1 &
#ifndef SLICE_MODEL
                                              +fly(i,j  )*DXV-fly(i,j-1)*DXVJM1 &
#endif
                                             )*ARCD1 )/Di(i,j)
         end if
      end do
   end do
!$OMP END DO

!  Calculating and applying the flux limiter
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1) then
            cmin(i,j)= 10000*_ONE_
            cmax(i,j)=-10000*_ONE_
! Calculate min and max of all values around one point
            do ii=-1,1
               do jj=-1,1
                  if (az(i+ii,j+jj).ge.1) then
                     x=min(f(i+ii,j+jj),fi(i+ii,j+jj))
                     if (x .lt. cmin(i,j)) cmin(i,j)=x
                     x=max(f(i+ii,j+jj),fi(i+ii,j+jj))
                     if (x .gt. cmax(i,j)) cmax(i,j)=x
                  end if
               end do
            end do

! max (Cu) and min (Cl) possible concentration after a time step
            CExx=(                                              &
                   ( min(fhx(i  ,j  )-flx(i  ,j  ),_ZERO_)      &
                    -max(fhx(i-1,j  )-flx(i-1,j  ),_ZERO_))/DXU &
#ifndef SLICE_MODEL
                  +( min(fhy(i  ,j  )-fly(i  ,j  ),_ZERO_)      &
                    -max(fhy(i  ,j-1)-fly(i  ,j-1),_ZERO_))/DYV &
#endif
                 )
            Cu=(fi(i,j)*Di(i,j)-dt*CExx)/Di(i,j)

            CExx=(                                              &
                   ( max(fhx(i  ,j  )-flx(i  ,j  ),_ZERO_)      &
                    -min(fhx(i-1,j  )-flx(i-1,j  ),_ZERO_))/DXU &
#ifndef SLICE_MODEL
                  +( max(fhy(i  ,j  )-fly(i  ,j  ),_ZERO_)      &
                    -min(fhy(i  ,j-1)-fly(i  ,j-1),_ZERO_))/DYV &
#endif
                 )
            Cl=(fi(i,j)*Di(i,j)-dt*CExx)/Di(i,j)

! calculating the maximum limiters rp and rm for each conc. cell
            if (Cu .eq. fi(i,j)) then
               rp(i,j)=_ZERO_
            else
               rp(i,j)=min((cmax(i,j)-fi(i,j))/(Cu-fi(i,j)),_ONE_)
            end if
            if (Cl .eq. fi(i,j)) then
               rm(i,j)=_ZERO_
            else
               rm(i,j)=min((fi(i,j)-cmin(i,j))/(fi(i,j)-Cl),_ONE_)
            end if
         end if
      end do
   end do
!$OMP END DO

!  Limiters for the u-fluxes (fac)
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin-1,imax
         if (fhx(i,j)-flx(i,j) .ge. _ZERO_) then
            fac=min(rm(i,j),rp(i+1,j))
         else
            fac=min(rm(i+1,j),rp(i,j))
         end if
         fhx(i,j) = (_ONE_-fac)*flx(i,j) + fac*fhx(i,j)
         if ((AH .gt. _ZERO_) .and. (az(i,j) .gt. 0) .and. (az(i+1,j) .gt. 0)) then
            fhx(i,j) = fhx(i,j) - AH*( f(i+1,j)                                 &
                                      -f(i  ,j))/DXU*_HALF_*(Dn(i+1,j)+Dn(i,j))
         end if
      end do
   end do
!$OMP END DO

#ifndef SLICE_MODEL
!  Limiters for the v-fluxes (fac)
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-1,jmax
      do i=imin,imax
         if (fhy(i,j)-fly(i,j) .ge. _ZERO_) then
            fac=min(rm(i,j),rp(i,j+1))
         else
            fac=min(rm(i,j+1),rp(i,j))
         end if
         fhy(i,j) = (_ONE_-fac)*fly(i,j) + fac*fhy(i,j)
         if ((AH .gt. _ZERO_) .and. (az(i,j) .gt. 0) .and. (az(i,j+1) .gt. 0)) then
            fhy(i,j) = fhy(i,j) - AH*( f(i,j+1)                                 &
                                      -f(i,j  ))/DYV*_HALF_*(Dn(i,j+1)+Dn(i,j))
         end if
      end do
   end do
!$OMP END DO
#endif

! Doing the full advection in one step
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1)  then
!           CAUTION: Di(i,j) already calculated above
            advn = (                                     &
                      fhx(i,j)*DYU - fhx(i-1,j  )*DYUIM1 &
#ifndef SLICE_MODEL
                    + fhy(i,j)*DXV - fhy(i  ,j-1)*DXVJM1 &
#endif
                   )*ARCD1
            adv(i,j) = adv(i,j) + advn
            if (present(nosplit_finalise)) then
!              Note (KK): do not update f in case of nosplit_finalise=.false. !!!
               if (nosplit_finalise) then
                  f(i,j) = ( Do(i,j)*f(i,j) - dt*adv(i,j) ) / Di(i,j)
               end if
            else
               f(i,j) = ( Dio(i,j)*f(i,j) - dt*advn ) / Di(i,j)
            end if
            if (update_f) then
!              Force monotonicity, this is needed here for correcting truncations errors:
               if (f(i,j).gt.cmax(i,j)) f(i,j)=cmax(i,j)
               if (f(i,j).lt.cmin(i,j)) f(i,j)=cmin(i,j)
            end if
         end if
      end do
   end do
!$OMP END DO

!$OMP END PARALLEL

#ifdef USE_ALLOCATED_ARRAYS
#ifdef FORTRAN90
   deallocate(flx,stat=rc)    ! work array
   if (rc /= 0) stop 'adv_fct_2dh: Error de-allocating memory (flx)'
   deallocate(fhx,stat=rc)    ! work array
   if (rc /= 0) stop 'adv_fct_2dh: Error de-allocating memory (fhx)'
#ifndef SLICE_MODEL
   deallocate(fly,stat=rc)    ! work array
   if (rc /= 0) stop 'adv_fct_2dh: Error de-allocating memory (fly)'
   deallocate(fhy,stat=rc)    ! work array
   if (rc /= 0) stop 'adv_fct_2dh: Error de-allocating memory (fhy)'
#endif
   deallocate(fi,stat=rc)    ! work array
   if (rc /= 0) stop 'adv_fct_2dh: Error de-allocating memory (fi)'
   deallocate(rp,stat=rc)    ! work array
   if (rc /= 0) stop 'adv_fct_2dh: Error de-allocating memory (rp)'
   deallocate(rm,stat=rc)    ! work array
   if (rc /= 0) stop 'adv_fct_2dh: Error de-allocating memory (rm)'
#endif
#endif

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
