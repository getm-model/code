#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:  adv_arakawa_j7_2dh - 2D Arakawa J7 advection 2d
! \label{sec-arakawa-j7-2dh-adv}
!
! !INTERFACE:
   subroutine adv_arakawa_j7_2dh(dt,f,Di,adv,vfU,vfV,Do,Dn,DU,DV,  &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                 dxv,dyu,dxu,dyv,arcd1,            &
#endif
                                 az,AH,                            &
                                 mask_uflux,mask_vflux,mask_xflux, &
                                 nosplit_finalise)
!  Note (KK): Keep in sync with interface in advection.F90
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
#if !( defined(SPHERICAL) || defined(CURVILINEAR) )
   use domain, only: dx,dy,ard1
#endif
   use advection, only: uflux,vflux
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)                        :: dt,AH
   REALTYPE,dimension(E2DFIELD),intent(in)    :: vfU,vfV,Do,Dn,DU,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
   REALTYPE,dimension(:,:),pointer,intent(in) :: dxu,dyu
   REALTYPE,dimension(_IRANGE_HALO_,_JRANGE_HALO_-1),intent(in) :: dxv,dyv
   REALTYPE,dimension(E2DFIELD),intent(in)    :: arcd1
#endif
   integer,dimension(E2DFIELD),intent(in)     :: az
   logical,dimension(:,:),pointer,intent(in)  :: mask_uflux,mask_xflux
   logical,dimension(_IRANGE_HALO_,_JRANGE_HALO_-1),intent(in) :: mask_vflux
   logical,intent(in),optional                :: nosplit_finalise
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(inout) :: f,Di,adv
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   logical                      :: use_AH
   integer                      :: i,j
   REALTYPE                     :: Dio,advn
   REALTYPE,dimension(E2DFIELD) :: flux_e,flux_n,flux_ne,flux_nw
   REALTYPE,dimension(E2DFIELD) :: f_e,f_n,f_ne,f_nw
   REALTYPE,parameter           :: one3rd = _ONE_/_THREE_
   REALTYPE,parameter           :: one6th = one3rd/_TWO_
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'adv_arakawa_j7_2dh() # ',Ncall
#endif

#ifdef SLICE_MODEL
 stop 'adv_arakawa_j7_2dh(): Do not use adv_arakawa_j7_2dh in SLICE_MODEL mode'
#endif

   use_AH = (AH .gt. _ZERO_)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,Dio,advn)

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO-1
      do i=imin-HALO,imax+HALO-1

         if (mask_uflux(i,j) .and. j.ne.jmin-HALO) then
            flux_e(i,j) = one3rd*( vfU(i,j-1) + vfU(i,j) )
            if (use_AH) then
!              Horizontal diffusion
               uflux(i,j) = - AH*DU(i,j)*(f(i+1,j)-f(i  ,j))/DXU
            end if
         else
            flux_e(i,j) = _ZERO_
         end if
         f_e(i,j) = _HALF_*( f(i,j) + f(i+1,j) )

         if (mask_vflux(i,j) .and. i.ne.imin-HALO) then
            flux_n(i,j) = one3rd*( vfV(i-1,j) + vfV(i,j) )
            if (use_AH) then
!              Horizontal diffusion
               vflux(i,j) = - AH*DV(i,j)*(f(i,j+1)-f(i  ,j))/DYV
            end if
         else
            flux_n(i,j) = _ZERO_
         end if
         f_n(i,j) = _HALF_*( f(i,j) + f(i,j+1) )

         if (mask_xflux(i,j)) then
            flux_ne(i,j) = one6th*( vfV(i,j) + vfU(i,j) )
            flux_nw(i,j) = one6th*( vfV(i,j) - vfU(i,j) )
         else
            flux_ne(i,j) = _ZERO_
            flux_nw(i,j) = _ZERO_
         end if
         f_ne(i,j) = _HALF_*( f(i,j) + f(i+1,j+1) )
         f_nw(i,j) = _HALF_*( f(i+1,j) + f(i,j+1) )

      end do
   end do
!$OMP END DO

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO+1,jmax+HALO-1
      do i=imin-HALO+1,imax+HALO-1
         if (az(i,j) .eq. 1)  then
            Dio = Di(i,j)
            Di(i,j) =  Dio - dt*(                                    &
                                   flux_e (i  ,j) - flux_e (i-1,j  ) &
                                 + flux_n (i  ,j) - flux_n (i  ,j-1) &
                                 + flux_ne(i  ,j) - flux_ne(i-1,j-1) &
                                 + flux_nw(i-1,j) - flux_nw(i  ,j-1) &
                                )*ARCD1
            advn = (                                                              &
                      flux_e (i  ,j)*f_e (i  ,j) - flux_e (i-1,j  )*f_e (i-1,j  ) &
                    + flux_n (i  ,j)*f_n (i  ,j) - flux_n (i  ,j-1)*f_n (i  ,j-1) &
                    + flux_ne(i  ,j)*f_ne(i  ,j) - flux_ne(i-1,j-1)*f_ne(i-1,j-1) &
                    + flux_nw(i-1,j)*f_nw(i-1,j) - flux_nw(i  ,j-1)*f_nw(i  ,j-1) &
                   )*ARCD1
            if (use_AH) then
!              Horizontal diffusion
               advn = advn + (  uflux(i,j)*DYU - uflux(i-1,j)*DYUIM1 &
                              + vflux(i,j)*DXV - vflux(i,j-1)*DXVJM1 )*ARCD1
            end if
            adv(i,j) = adv(i,j) + advn
            if (present(nosplit_finalise)) then
!              Note (KK): do not update f in case of nosplit_finalise=.false. !!!
               if (nosplit_finalise) then
                  f(i,j) = ( Do(i,j)*f(i,j) - dt*adv(i,j) ) / Di(i,j)
               end if
            else
               f(i,j) = ( Dio*f(i,j) - dt*advn ) / Di(i,j)
            end if
         end if
      end do
   end do
!$OMP END DO

!$OMP END PARALLEL

#ifdef DEBUG
   write(debug,*) 'Leaving adv_arakawa_j7_2dh()'
   write(debug,*)
#endif
   return
   end subroutine adv_arakawa_j7_2dh
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
