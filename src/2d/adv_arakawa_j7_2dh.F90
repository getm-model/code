#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:  adv_arakawa_j7_2dh - 2DH Arakawa J7 advection of 2D quantities
! \label{sec-arakawa-j7-2dh-adv}
!
! !INTERFACE:
   subroutine adv_arakawa_j7_2dh(dt,f,fi,Di,adv,vfU,vfV,Dn,DU,DV,  &
#if defined(SPHERICAL) || defined(CURVILINEAR)
                                 dxv,dyu,dxu,dyv,arcd1,            &
#endif
                                 AH,az,                            &
                                 mask_uflux,mask_vflux,mask_xflux)
!  Note (KK): Keep in sync with interface in advection.F90
!
! !DESCRIPTION:
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
   REALTYPE,intent(in)                               :: dt,AH
   REALTYPE,dimension(E2DFIELD),target,intent(in)    :: f
   REALTYPE,dimension(E2DFIELD),intent(in)           :: vfU,vfV,Dn,DU,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
   REALTYPE,dimension(:,:),pointer,intent(in)        :: dxu,dyu
   REALTYPE,dimension(_IRANGE_HALO_,_JRANGE_HALO_-1),intent(in) :: dxv,dyv
   REALTYPE,dimension(E2DFIELD),intent(in)           :: arcd1
#endif
   integer,dimension(E2DFIELD),intent(in)            :: az
   logical,dimension(:,:),pointer,intent(in)         :: mask_uflux,mask_xflux
   logical,dimension(_IRANGE_HALO_,_JRANGE_HALO_-1),intent(in) :: mask_vflux
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),target,intent(inout) :: fi,Di,adv
!
! !LOCAL VARIABLES:
   logical                      :: use_AH
   integer                      :: i,j,matsuno_it
   REALTYPE                     :: Dio,advn
   REALTYPE,dimension(:,:),pointer :: faux,p_fiaux,p_Diaux,p_advaux
   REALTYPE,dimension(E2DFIELD) :: flux_e,flux_n,flux_ne,flux_nw
   REALTYPE,dimension(E2DFIELD) :: f_e,f_n,f_ne,f_nw
   REALTYPE,dimension(E2DFIELD),target :: fiaux,Diaux,advaux
   REALTYPE,dimension(E2DFIELD) :: uflux,vflux
   REALTYPE,parameter           :: one3rd = _ONE_/_THREE_
   REALTYPE,parameter           :: one6th = one3rd/_TWO_
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
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
         else
            flux_e(i,j) = _ZERO_
         end if

         if (mask_vflux(i,j) .and. i.ne.imin-HALO) then
            flux_n(i,j) = one3rd*( vfV(i-1,j) + vfV(i,j) )
         else
            flux_n(i,j) = _ZERO_
         end if

         if (mask_xflux(i,j)) then
            flux_ne(i,j) = one6th*( vfV(i,j) + vfU(i,j) )
            flux_nw(i,j) = one6th*( vfV(i,j) - vfU(i,j) )
         else
            flux_ne(i,j) = _ZERO_
            flux_nw(i,j) = _ZERO_
         end if

      end do
   end do
!$OMP END DO

   faux     => f
   p_fiaux  => fiaux
   p_Diaux  => Diaux
   p_advaux => advaux

   do matsuno_it = 0,1

      if (matsuno_it .eq. 1) then
         p_fiaux  => fi
         p_Diaux  => Di
         p_advaux => adv
      end if

!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO,jmax+HALO-1
         do i=imin-HALO,imax+HALO-1

            f_e (i,j) = _HALF_*( faux(i  ,j) + faux(i+1,j  ) )
            f_n (i,j) = _HALF_*( faux(i  ,j) + faux(i  ,j+1) )
            f_ne(i,j) = _HALF_*( faux(i  ,j) + faux(i+1,j+1) )
            f_nw(i,j) = _HALF_*( faux(i+1,j) + faux(i  ,j+1) )

!           Horizontal diffusion
            if (use_AH) then
               if (mask_uflux(i,j) .and. j.ne.jmin-HALO) then
                  uflux(i,j) = - AH*DU(i,j)*(faux(i+1,j)-faux(i  ,j))/DXU
               else
                  uflux(i,j) = _ZERO_
               end if
               if (mask_vflux(i,j) .and. i.ne.imin-HALO) then
                  vflux(i,j) = - AH*DV(i,j)*(faux(i,j+1)-faux(i  ,j))/DYV
               else
                  vflux(i,j) = _ZERO_
               end if
            end if

         end do
      end do
!$OMP END DO

!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO+1,jmax+HALO-1
         do i=imin-HALO+1,imax+HALO-1
            if (az(i,j) .eq. 1)  then
               Dio = Di(i,j)
               p_Diaux(i,j) = Dio - dt*(                                    &
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
!                 Horizontal diffusion
                  advn = advn + (  uflux(i,j)*DYU - uflux(i-1,j)*DYUIM1 &
                                 + vflux(i,j)*DXV - vflux(i,j-1)*DXVJM1 )*ARCD1
               end if
               p_fiaux(i,j) = ( Dio*fi(i,j) - dt*advn ) / p_Diaux(i,j)
               p_advaux(i,j) = adv(i,j) + advn
            end if
         end do
      end do
!$OMP END DO

      faux => p_fiaux

   end do

!$OMP END PARALLEL

#ifdef DEBUG
   write(debug,*) 'Leaving adv_arakawa_j7_2dh()'
   write(debug,*)
#endif
   return
   end subroutine adv_arakawa_j7_2dh
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2012 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
