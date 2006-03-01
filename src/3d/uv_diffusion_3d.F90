!$Id: uv_diffusion_3d.F90,v 1.8 2006-03-01 15:54:08 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: uv_diffusion_3d - hor.\ momentum diffusion
! \label{sec-uv-diffusion-3d}
!
! !INTERFACE:
   subroutine uv_diffusion_3d(Am)
!
! !DESCRIPTION:
!
! Here, the horizontal diffusion terms are discretised in 
! momentum-conserving form 
! as well. For simplicity,
! this shown here for these terms as they appear in equations
! (\ref{uEqviCurvi}) and (\ref{vEqviCurvi}), 
! i.e.\ without multiplying them by $mn$.  
! 
! First horizontal diffusion term in (\ref{uEqviCurvi}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! \partial_{\cal X}\left(\frac{2A_Mh_k}{n} m\partial_{\cal X}u_k\right)_{i,j,k}
! \approx
! \\ \\ \displaystyle \quad
! 2A^M_{i+1,j,k} \Delta y_{i+1,j}^c h_{i+1,j,k}^c
! \frac{u_{i+1,j,k}-u_{i,j,k}}{\Delta x_{i+1,j}^c}
! -
! 2A^M_{i,j,k} \Delta y_{i,j}^c h_{i,j,k}^c
! \frac{u_{i,j,k}-u_{i-1,j,k}}{\Delta x_{i,j}^c}
! \end{array}
! \end {equation}
! 
! Second horizontal diffusion term in (\ref{uEqviCurvi}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! \partial_{\cal Y}\left(\frac{A_Mh_k}{m}
! (n\partial_{\cal Y}u_k+m\partial_{\cal X}v_k)\right)_{i,j,k}
! \approx 
! \\ \\ \displaystyle \quad 
! \,\,\,\,\,\frac14\left(
! A^M_{i+1,j+1,k}+A^M_{i,j+1,k}+A^M_{i,j,k}+A^M_{i+1,j,k}\right)
! h^+_{i,j,k}\Delta x^+_{i,j}
! \\ \\ \displaystyle \quad 
! \,\,\,\,\,\times\bigg(
! \frac{u_{i,j+1,k}-u_{i,j,k}}{\Delta y^+_{i,j}}
! -\frac{v_{i+1,j,k}-v_{i,j,k}}{\Delta x^+_{i,j}}
! \bigg)
! \\ \\ \displaystyle \quad 
! -\frac14\left(
! A^M_{i+1,j,k}+A^M_{i,j,k}+A^M_{i,j-1,k}+A^M_{i+1,j-1,k}\right)
! h^+_{i,j-1,k}\Delta x^+_{i,j-1}
! \\ \\ \displaystyle \quad 
! \,\,\,\,\,\times\bigg(
! \frac{u_{i,j,k}-u_{i,j-1,k}}{\Delta y^+_{i,j-1}}
! -\frac{v_{i+1,j-1,k}-v_{i,j-1,k}}{\Delta x^+_{i,j-1}}
! \bigg)
! \end{array}
! \end {equation}
! 
! First horizontal diffusion term in (\ref{vEqviCurvi}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! \partial_{\cal Y}\left(\frac{2A_Mh_k}{m} n\partial_{\cal Y}v_k\right)_{i,j,k}
! \approx
! \\ \\ \displaystyle \quad
! 2A^M_{i,j+1,k} \Delta x_{i,j+1}^c h_{i,j+1,k}^c
! \frac{v_{i,j+1,k}-v_{i,j,k}}{\Delta y_{i,j+1}^c}
! -
! 2A^M_{i,j,k} \Delta x_{i,j}^c h_{i,j,k}^c
! \frac{v_{i,j,k}-v_{i,j-1,k}}{\Delta y_{i,j}^c}
! \end{array}
! \end {equation}
! 
! Second horizontal diffusion term in (\ref{vEqviCurvi}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! \partial_{\cal X}\left(\frac{A_Mh_k}{n}
! (n\partial_{\cal Y}u_k+m\partial_{\cal X}v_k)\right)_{i,j,k}
! \approx 
! \\ \\ \displaystyle \quad 
! \,\,\,\,\,\frac14\left(
! A^M_{i+1,j+1,k}+A^M_{i,j+1,k}+A^M_{i,j,k}+A^M_{i+1,j,k}\right)
! h^+_{i,j,k}\Delta x^+_{i,j}
! \\ \\ \displaystyle \quad 
! \,\,\,\,\,\times\bigg(
! \frac{u_{i,j+1,k}-u_{i,j,k}}{\Delta y^+_{i,j}}
! -\frac{v_{i+1,j,k}-v_{i,j,k}}{\Delta x^+_{i,j}}
! \bigg)
! \\ \\ \displaystyle \quad 
! -\frac14\left(
! A^M_{i,j+1,k}+A^M_{i-1,j+1,k}+A^M_{i-1,j,k}+A^M_{i,j,k}\right)
! h^+_{i-1,j,k}\Delta x^+_{i-1,j}
! \\ \\ \displaystyle \quad 
! \,\,\,\,\,\times\bigg(
! \frac{u_{i-1,j+1,k}-u_{i-1,j,k}}{\Delta y^+_{i-1,j}}
! -\frac{v_{i,j,k}-v_{i-1,j,k}}{\Delta x^+_{i-1,j}}
! \bigg)
! \end{array}
! \end {equation}
! 
! It is assumed here that the horizontal momentum diffusivities 
! $A^M_{i,j,k}$ are located
! on the T-points.  
!
! For the case of a slice model simulation (compiler option {\tt SLICE\_MODEL}
! activated) the diffusive fluxes in $y$-direction are set to zero.
!
! !USES:
   use domain, only: iimin,iimax,jjmin,jjmax,kmax,az,au,av,ax
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dyc,arud1,dxx,dyx,arvd1,dxc
#else
   use domain, only: dx,dy,ard1
#endif
   use variables_3d, only: kumin,kvmin,uu,vv,ww,hn,hun,hvn,uuEx,vvEx
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
  REALTYPE, intent(in) :: Am

! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,ii,jj,kk
   REALTYPE                  :: PP(iimin-1:iimax+1,jjmin-1:jjmax+1,1:kmax)
   REALTYPE                  :: www(0:kmax)
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'D3uvDiff() # ',Ncall
#endif

! Central for dx(2*Am*dx(uu^2/hun))
   do k=1,kmax
      do j=jjmin,jjmax
         do i=iimin,iimax+1          ! PP defined on T-points
            PP(i,j,k)=_ZERO_
            if (az(i,j) .ge. 1) then
               if (k .ge. kumin(i,j)) then
                  PP(i,j,k)=2.*Am*DYC*hn(i,j,k)               &
                      *(uu(i,j,k)/hun(i,j,k)-uu(i-1,j,k)/hun(i-1,j,k))/DXC
               end if
            end if
         end do
      end do
   end do
   do k=1,kmax
      do j=jjmin,jjmax         ! uuEx defined on U-points
         do i=iimin,iimax
            if (au(i,j) .ge. 1) then
               if (k .ge. kumin(i,j)) then
                  uuEx(i,j,k)=uuEx(i,j,k)-(PP(i+1,j,k)-PP(i,j,k))*ARUD1
               end if
            end if
         end do
      end do
   end do

#ifndef SLICE_MODEL
! Central for dy(Am*(dy(uu^2/hun)+dx(vv^2/hvn)))
   do k=1,kmax
      do j=jjmin-1,jjmax          ! PP defined on X-points
         do i=iimin,iimax
            PP(i,j,k)=_ZERO_
            if (ax(i,j) .ge. 1) then
               if (k .ge. kumin(i,j)) then
                  PP(i,j,k)=Am*0.5*(hun(i,j,k)+hun(i,j+1,k))*DXX  &
                      *((uu(i,j+1,k)/hun(i,j+1,k)-uu(i,j,k)/hun(i,j,k))/DYX &
                       +(vv(i+1,j,k)/hvn(i+1,j,k)-vv(i,j,k)/hvn(i,j,k))/DXX )
               end if
            end if
         end do
      end do
   end do
   do k=1,kmax
      do j=jjmin,jjmax
         do i=iimin,iimax
            if (au(i,j) .ge. 1) then
               if (k .ge. kumin(i,j)) then
                  uuEx(i,j,k)=uuEx(i,j,k)-(PP(i,j,k)-PP(i,j-1,k))*ARUD1
               end if
            end if
         end do
      end do
   end do
#endif

! Central for dx(Am*(dy(uu^2/hun)+dx(vv^2/hvn)))
   do k=1,kmax
      do j=jjmin,jjmax          ! PP defined on X-points
         do i=iimin-1,iimax
            PP(i,j,k)=_ZERO_
            if (ax(i,j) .ge. 1) then
               if (k .ge. kumin(i,j)) then
                  PP(i,j,k)=Am*0.5*(hvn(i+1,j,k)+hvn(i,j,k))*DXX  &
                      *((uu(i,j+1,k)/hun(i,j+1,k)-uu(i,j,k)/hun(i,j,k))/DYX &
                       +(vv(i+1,j,k)/hvn(i+1,j,k)-vv(i,j,k)/hvn(i,j,k))/DXX )
               end if
            end if
         end do
      end do
   end do
   do k=1,kmax
      do j=jjmin,jjmax          ! vvEx defined on V-points
         do i=iimin,iimax
            if (av(i,j) .ge. 1) then
               if (k .ge. kvmin(i,j)) then
                  vvEx(i,j,k)=vvEx(i,j,k)-(PP(i,j,k)-PP(i-1,j,k))*ARVD1
               end if
            end if
         end do
      end do
   end do

#ifndef SLICE_MODEL
! Central for dy(2*Am*dy(vv^2/hvn))
   do k=1,kmax
      do j=jjmin,jjmax+1
         do i=iimin,iimax          ! PP defined on T-points
            if (az(i,j) .ge. 1) then
               if (k .ge. kvmin(i,j)) then
                  PP(i,j,k)=2.*Am*DXC*hn(i,j,k)               &
                      *(vv(i,j,k)/hvn(i,j,k)-vv(i,j-1,k)/hvn(i,j-1,k))/DYC
               end if
            end if
         end do
      end do
   end do

   do k=1,kmax
      do j=jjmin,jjmax          ! vvEx defined on V-points
         do i=iimin,iimax
            if (av(i,j) .ge. 1) then
               if (k .ge. kvmin(i,j)) then
                  vvEx(i,j,k)=(vvEx(i,j,k)-(PP(i,j+1,k)-PP(i,j,k))*ARVD1)
               end if
            end if
         end do
      end do
   end do
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving uv_diffusion_3d()'
   write(debug,*)
#endif
   return
   end subroutine uv_diffusion_3d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
