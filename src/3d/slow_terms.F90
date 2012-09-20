#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: slow_terms - calculation of slow terms \label{sec-slow-terms}
!
! !INTERFACE:
   subroutine slow_terms(runtype)
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
   use variables_3d, only: idpdx_hs,idpdy_hs,idpdx_nh,idpdy_nh
#ifdef STRUCTURE_FRICTION
   use variables_3d, only: sf
#endif
   use m3d, only: calc_ip,ip_fac
   use m2d, only: uv_advect,uv_diffusion,bottom_friction
   use internal_pressure, only: calc_ipfull
   use nonhydrostatic, only: calc_hs2d,sbnh_filter
   use halo_zones, only: update_2d_halo,wait_halo,U_TAG,V_TAG
   use getm_timers, only: tic, toc, TIM_SLOWTERMS
!$ use omp_lib
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !INPUT PARAMETERS:
   integer,intent(in)        :: runtype
!
! !LOCAL VARIABLES:
   integer                      :: i,j,k
   REALTYPE                     :: vertsum
   REALTYPE,dimension(E2DFIELD) :: work2d
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'slow_terms() # ',Ncall
#endif
#ifdef SLICE_MODEL
   j = jmax/2 ! this MUST NOT be changed!!!
#endif
   call tic(TIM_SLOWTERMS)

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j)                                         &
!$OMP          PRIVATE(i,k,vertsum)
!$OMP SINGLE

!  Note (KK): for kmax=1 [uu|vv]=[U|V]int and only internal
!             pressure contributes to slow terms
   if (kmax .gt. 1) then

      call bottom_friction(Uint,Vint,Dun,Dvn,ru,rv)
      call uv_advect(Uint,Vint,Dun,Dvn)
      call uv_diffusion(0,Uint,Vint,Dn,Dun,Dvn) ! Has to be called after uv_advect.

!$OMP END SINGLE
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin,jmax
#endif
         do i=imin,imax

            if (au(i,j) .ge. 1) then

               vertsum = -UEx(i,j)
               do k=kumin(i,j),kmax
                  vertsum = vertsum + uuEx(i,j,k)
               end do
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
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO
!$OMP SINGLE

   else

      SlUx = _ZERO_
#ifndef SLICE_MODEL
      SlVx = _ZERO_
#endif

   end if

   if (runtype .ge. 3) then ! hs PGF
!$OMP END SINGLE
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin,jmax
#endif
         do i=imin,imax

            if (au(i,j) .ge. 1) then
               vertsum = _ZERO_
               do k=1,kmax
                  vertsum = vertsum - idpdx_hs(i,j,k)
               end do
               SlUx(i,j) = SlUx(i,j) + ip_fac*vertsum
            end if

#ifndef SLICE_MODEL
            if (av(i,j) .ge. 1) then
               vertsum = _ZERO_
               do k=1,kmax
                  vertsum = vertsum - idpdy_hs(i,j,k)
               end do
               SlVx(i,j) = SlVx(i,j) + ip_fac*vertsum
            end if
#endif

         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO
!$OMP SINGLE
   end if

   if (.not. calc_hs2d) then ! nh PGF

!$OMP END SINGLE
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin,jmax
#endif
         do i=imin,imax
            if (au(i,j) .ge. 1) then
               vertsum = _ZERO_
               do k=1,kmax
                  vertsum = vertsum - idpdx_nh(i,j,k)
               end do
               work2d(i,j) = ip_fac*vertsum
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO
!$OMP SINGLE
      if (sbnh_filter) then
         call update_2d_halo(work2d,work2d,au,imin,jmin,imax,jmax,U_TAG)
         call wait_halo(U_TAG)
#ifndef SLICE_MODEL
!$OMP END SINGLE
!$OMP DO SCHEDULE(RUNTIME)
         do j=jmin,jmax
#endif
            call filter_1d(imin-HALO  ,imax+HALO  ,au(:,j),work2d(:,j), &
                           imin-HALO+1,imax+HALO-1,work2d(:,j))
#ifndef SLICE_MODEL
         end do
!$OMP END DO
!$OMP SINGLE
#endif
      end if
      SlUx = SlUx + work2d

#ifndef SLICE_MODEL
!$OMP END SINGLE
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin,jmax
         do i=imin,imax
            if (av(i,j) .ge. 1) then
               vertsum = _ZERO_
               do k=1,kmax
                  vertsum = vertsum - idpdy_nh(i,j,k)
               end do
               work2d(i,j) = ip_fac*vertsum
            end if
         end do
      end do
!$OMP END DO
!$OMP SINGLE
      if (sbnh_filter) then
         call update_2d_halo(work2d,work2d,au,imin,jmin,imax,jmax,V_TAG)
         call wait_halo(V_TAG)
!$OMP END SINGLE
!$OMP DO SCHEDULE(RUNTIME)
         do i=imin,imax
            call filter_1d(jmin-HALO  ,jmax+HALO  ,av(i,:),work2d(i,:), &
                           jmin-HALO+1,jmax+HALO-1,work2d(i,:))
         end do
!$OMP END DO
!$OMP SINGLE
      end if
      SlVx = SlVx + work2d
#endif

   end if

!$OMP END PARALLEL

#ifdef SLICE_MODEL
   SlUx(:,j+1) = SlUx(:,j)
   Slru(:,j+1) = Slru(:,j)
   SlVx(:,j-1) = SlVx(:,j)
   SlVx(:,j+1) = SlVx(:,j)
   Slrv(:,j-1) = Slrv(:,j)
   Slrv(:,j+1) = Slrv(:,j)
#endif

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
