#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: slow_terms - calculation of slow terms \label{sec-slow-terms}
!
! !INTERFACE:
   subroutine slow_terms
!
! !DESCRIPTION:
!
! Here, the calculation of the so-called slow terms (which are the
! interaction terms between the barotropic and the baroclinic mode) is
! completed. The mathematical form of these slow terms is given by
! equations (\ref{Slowfirst}) - (\ref{Slowlast}), see section
! \ref{SectionVerticalIntegrated}.
! These calculations have been prepared in the routines
! {\tt slow\_bottom\_friction}, {\tt slow\_advection} and
! {\tt slow\_diffusion}.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax,au,av
   use variables_2d, only: Uint,Vint,UEx,VEx,Slru,Slrv,SlUx,SlVx,ru,rv
   use variables_3d, only: kumin,kvmin,uu,vv,hun,hvn,Dn,Dun,Dvn
   use variables_3d, only: uuEx,vvEx,rru,rrv
   use variables_3d, only: idpdx_m3d=>idpdx,idpdy_m3d=>idpdy,idpdx_hs,idpdy_hs,idpdx_nh,idpdy_nh
#ifdef STRUCTURE_FRICTION
   use variables_3d, only: sf
#endif
   use m3d, only: calc_ip,ip_fac
   use m2d, only: uv_advect,uv_diffusion,bottom_friction
   use internal_pressure, only: calc_ipfull
   use nonhydrostatic, only: nonhyd_method,calc_hs2d
   use halo_zones, only: update_2d_halo,wait_halo,U_TAG
   use getm_timers, only: tic, toc, TIM_SLOWTERMS
!$ use omp_lib
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   logical,save              :: first=.true.,calc_slowip
   integer                   :: i,j,k
   REALTYPE                  :: vertsum
   REALTYPE,dimension(:,:,:),pointer,save :: idpdx,idpdy
   REALTYPE,dimension(E2DFIELD) :: work2d
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'slow_terms() # ',Ncall
#endif
   call tic(TIM_SLOWTERMS)

   if (first) then
      calc_slowip = calc_ip
      idpdx => idpdx_m3d
      idpdy => idpdy_m3d
      if (calc_ipfull) then
         if (calc_hs2d) then
            idpdx => idpdx_hs
            idpdy => idpdy_hs
         end if
      else if (nonhyd_method .eq. 1) then
         if (calc_hs2d) then
            calc_slowip = .false.
         end if
      end if
      first = .false.
   end if

   work2d = _ZERO_

   if (allocated(idpdx_nh)) then
      do k=1,kmax
         work2d = work2d + idpdx_nh(:,:,k)
      end do

      call update_2d_halo(work2d,work2d,au,imin,jmin,imax,jmax,U_TAG)
      call wait_halo(U_TAG)

      do j=jmin,jmax
         call filter_1d(imin-HALO  ,imax+HALO  ,au(:,j),work2d(:,j), &
                        imin-HALO+1,imax+HALO-1,work2d(:,j))
      end do
   end if
   if (allocated(idpdx_hs)) then
      do k=1,kmax
         work2d = work2d + idpdx_hs(:,:,k)
      end do
   end if

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,vertsum)

   if (kmax .gt. 1) then

      call bottom_friction(Uint,Vint,Dun,Dvn,ru,rv)
      call uv_advect(Uint,Vint,Dun,Dvn)
      call uv_diffusion(0,Uint,Vint,Dn,Dun,Dvn) ! Has to be called after uv_advect.

!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin,imax

            if (au(i,j) .ge. 1) then

               vertsum = -UEx(i,j)
               do k=kumin(i,j),kmax
                  vertsum = vertsum + uuEx(i,j,k)
               end do
!               if (calc_slowip) then
!                  do k=kumin(i,j),kmax
!                     vertsum = vertsum - ip_fac*idpdx(i,j,k)
!                  end do
!               end if
               vertsum = vertsum - ip_fac*work2d(i,j)
               SlUx(i,j) = vertsum

#ifdef NO_SLR
               STDERR 'NO_SLR U'
               Slru(i,j)= _ZERO_
#else
               k=kumin(i,j)
               Slru(i,j) =   rru(i,j)*uu(i,j,k)/hun(i,j,k) &
                           - ru(i,j)*Uint(i,j)/Dun(i,j)
#endif

#ifdef STRUCTURE_FRICTION
               do k=kumin(i,j),kmax
                  Slru(i,j)=Slru(i,j)+uu(i,j,k)*_HALF_*(sf(i,j,k)+sf(i+1,j,k))
               end do
#endif

            end if

            if (av(i,j) .ge. 1) then

               vertsum = -VEx(i,j)
               do k=kvmin(i,j),kmax
                  vertsum = vertsum + vvEx(i,j,k)
               end do
#ifndef SLICE_MODEL
               if (calc_slowip) then
                  do k=kvmin(i,j),kmax
                     vertsum = vertsum - ip_fac*idpdy(i,j,k)
                  end do
               end if
#endif
               SlVx(i,j) = vertsum

#ifdef NO_SLR
               STDERR 'NO_SLR V'
               Slrv(i,j)= _ZERO_
#else
               k=kvmin(i,j)
               Slrv(i,j) =   rrv(i,j)*vv(i,j,k)/hvn(i,j,k) &
                           - rv(i,j)*Vint(i,j)/Dvn(i,j)
#endif

#ifdef STRUCTURE_FRICTION
               do k=kvmin(i,j),kmax
                  Slrv(i,j)=Slrv(i,j)+vv(i,j,k)*_HALF_*(sf(i,j,k)+sf(i,j+1,k))
               end do
#endif

            end if

         end do
      end do
!$OMP END DO

   else if (calc_slowip) then

!     Note (KK): here kmax=1, thus [uu|vv]=[U|V]int and only internal
!                pressure contributes to slow terms

!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin,imax
            if (au(i,j) .ge. 1) then
               SlUx(i,j) = - ip_fac*idpdx(i,j,1)
            end if
#ifndef SLICE_MODEL
            if (av(i,j) .ge. 1) then
               SlVx(i,j) = - ip_fac*idpdy(i,j,1)
            end if
#endif
         end do
      end do
!$OMP END DO

   end if

!$OMP END PARALLEL

   call toc(TIM_SLOWTERMS)
#ifdef DEBUG
   write(debug,*) 'Leaving slow_terms()'
   write(debug,*)
#endif
   return
   end subroutine slow_terms
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
