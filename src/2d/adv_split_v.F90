#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  adv_split_v - meridional advection of 2D quantities \label{sec-v-split-adv}
!
! !INTERFACE:
   subroutine adv_split_v(dt,f,fi,Di,adv,V,DV,   &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                          dxv,dyv,arcd1,         &
#endif
                          splitfac,scheme,AH,    &
                          mask_flux,mask_update, &
                          nvd)
!  Note (KK): Keep in sync with interface in advection.F90
!
! !DESCRIPTION:
!
! Executes an advection step in meridional direction for a 2D quantity
! in analogy to routine {\tt adv\_u\_split} (see section
! \ref{sec-u-split-adv} on page \pageref{sec-u-split-adv}).
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
#if !( defined(SPHERICAL) || defined(CURVILINEAR) )
   use domain, only: dx,dy,ard1
#endif
   use advection, only: adv_interfacial_reconstruction
   use advection, only: UPSTREAM
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                                          :: dt,splitfac,AH
   REALTYPE,dimension(E2DFIELD),intent(in)                      :: f,V,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
   REALTYPE,dimension(_IRANGE_HALO_,_JRANGE_HALO_-1),intent(in) :: dxv,dyv
   REALTYPE,dimension(E2DFIELD),intent(in)                      :: arcd1
#endif
   integer,intent(in)                                           :: scheme
   logical,dimension(_IRANGE_HALO_,_JRANGE_HALO_-1),intent(in)  :: mask_flux
   logical,dimension(E2DFIELD),intent(in)                       :: mask_update
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(inout)                   :: fi,Di,adv
   REALTYPE,dimension(:,:),pointer,intent(inout)                :: nvd
!
! !LOCAL VARIABLES:
   REALTYPE,dimension(E2DFIELD) :: vflux,vflux2
   logical            :: use_AH,calc_nvd
   integer            :: i,j,jsub
   REALTYPE           :: dti,fio,Dio,advn,adv2n,cfl,fuu,fu,fd
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'adv_split_v() # ',Ncall
#endif

   if (scheme .eq. UPSTREAM) then
      jsub = 0
   else
      jsub = 1
   end if

   use_AH = (AH .gt. _ZERO_)
   calc_nvd = associated(nvd)
   dti = splitfac*dt

!$OMP PARALLEL DEFAULT(SHARED)                                  &
!$OMP          PRIVATE(i,j,fio,Dio,advn,adv2n,cfl,fuu,fd)

! Calculating v-interface fluxes !
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO+jsub,jmax+HALO-1-jsub
      do i=imin-HALO,imax+HALO
         if (mask_flux(i,j)) then
!           Note (KK): exclude y-advection of v across N/S open bdys
            if (V(i,j) .gt. _ZERO_) then
               cfl = V(i,j)/DV(i,j)*dti/DYV
               fu = f(i,j  )               ! central
               if (mask_flux(i,j-1)) then
                  fuu = f(i,j-1)            ! upstream
               else
                  fuu = fu
               end if
               fd = f(i,j+1)            ! downstream
            else
               cfl = -V(i,j)/DV(i,j)*dti/DYV
               fu = f(i,j+1)               ! central
               if (mask_flux(i,j+1)) then
                  fuu = f(i,j+2)            ! upstream
               else
                  fuu = fu
               end if
               fd = f(i,j  )            ! downstream
            end if
            fu = adv_interfacial_reconstruction(scheme,cfl,fuu,fu,fd)
            vflux (i,j) = V(i,j)*fu
            vflux2(i,j) = vflux(i,j)*fu
            if (use_AH) then
!              Horizontal diffusion
               vflux(i,j) = vflux(i,j) - AH*DV(i,j)*(f(i,j+1)-f(i,j  ))/DYV
            end if
         else
            vflux (i,j) = _ZERO_
            vflux2(i,j) = _ZERO_
         end if
      end do
   end do
!$OMP END DO

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO+1+jsub,jmax+HALO-1-jsub
      do i=imin-HALO,imax+HALO
         if (mask_update(i,j)) then
!           Note (KK): exclude y-advection of tracer and v across N/S open bdys
            fio = fi(i,j)
            Dio = Di(i,j)
            Di(i,j) =  Dio - dti*( V(i,j  )*DXV           &
                                  -V(i,j-1)*DXVJM1)*ARCD1
            advn = splitfac*( vflux(i,j  )*DXV           &
                             -vflux(i,j-1)*DXVJM1)*ARCD1
            fi(i,j) = ( Dio*fio - dt*advn ) / Di(i,j)
            adv(i,j) = adv(i,j) + advn
            if (calc_nvd) then
               adv2n = splitfac*( vflux2(i,j  )*DXV           &
                                 -vflux2(i,j-1)*DXVJM1)*ARCD1
!               nvd(i,j) = nvd(i,j)                                       &
!                         -((Di(i,j)*fi(i,j)**2 - Dio*fio**2)/dt + adv2n)
               nvd(i,j) = ( Dio*nvd(i,j)                                   &
                           -((Di(i,j)*fi(i,j)**2 - Dio*fio**2)/dt + adv2n) &
                          )/Di(i,j)
            end if
         end if
      end do
   end do
!$OMP END DO

!$OMP END PARALLEL

#ifdef DEBUG
   write(debug,*) 'Leaving adv_split_v()'
   write(debug,*)
#endif
   return
   end subroutine adv_split_v
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
