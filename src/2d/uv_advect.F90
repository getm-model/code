#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: uv_advect - 2D advection of momentum \label{sec-uv-advect}
!
! !INTERFACE:
   subroutine uv_advect(U,V,Dold,Dnew,DU,DV,numdis)

!  Note (KK): keep in sync with interface in m2d.F90
!
! !DESCRIPTION:
!
! Wrapper to prepare and do calls to {\tt do\_advection} (see section
! \ref{sec-do-advection} on page \pageref{sec-do-advection}) to
! calculate the advection terms of the depth-averaged velocities.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,az,au,av,ax
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxv,dyu,arcd1,arud1,arvd1
#else
   use domain, only: dx,dy,ard1
#endif
   use m2d, only: vel2d_adv_split,vel2d_adv_hor
   use variables_2d, only: dtm
   use variables_2d, only: UEx,VEx
   use variables_2d, only: do_numerical_analyses_2d
   use variables_2d, only: numdis_2d_old
   use advection, only: NOADV,UPSTREAM,J7,do_advection
   use halo_zones, only: update_2d_halo,wait_halo,U_TAG,V_TAG
   use getm_timers, only: tic,toc,TIM_UVADV,TIM_UVADVH
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(in)        :: U,V
   REALTYPE,dimension(E2DFIELD),target,intent(in) :: Dold,Dnew,DU,DV
!
! !OUTPUT PARAMETERS:
   REALTYPE,dimension(:,:),pointer,intent(out),optional :: numdis
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   logical                             :: calc_numdis
   integer                             :: i,j
   REALTYPE,dimension(E2DFIELD)        :: fadv,Uadv,Vadv,DUadv,DVadv,Dires
   REALTYPE,dimension(E2DFIELD),target :: Dadv,nvd
   REALTYPE,dimension(:,:),pointer     :: pDadv,p_nvd
#ifdef _NUMERICAL_ANALYSES_OLD_
   REALTYPE,dimension(E2DFIELD)        :: work2d
#endif
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'uv_advect() # ',Ncall
#endif
#ifdef SLICE_MODEL
   j = jmax/2 ! this MUST NOT be changed!!!
#endif
   call tic(TIM_UVADV)

   if (present(numdis)) then
      calc_numdis = associated(numdis)
   else
      calc_numdis = .false.
   end if
      
   if (calc_numdis) then
      p_nvd => nvd
   else
      p_nvd => null()
   end if
   

   if (vel2d_adv_hor .ne. NOADV) then

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j)                                         &
!$OMP          PRIVATE(i)


!     Here begins dimensional split advection for u-velocity

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin-HALO,jmax+HALO
#endif
         do i=imin-HALO,imax+HALO-1
!           the velocity to be transported
            fadv(i,j) = U(i,j)/DU(i,j)
            if (vel2d_adv_hor .ne. J7) then
!              Note (KK): Uadv defined on T-points (future U-points)
!                         Vadv defined on X-points (future V-points)
               Uadv(i,j) = _HALF_*( U(i,j) + U(i+1,j) )
               Vadv(i,j) = _HALF_*( V(i,j) + V(i+1,j) )
            end if
            DUadv(i,j) = _HALF_*( Dold(i+1,j) + Dnew(i+1,j) )
!           Note (KK): DV only valid until jmax+1
!                      therefore DVadv only valid until jmax+1
            DVadv(i,j) = _HALF_*( DV(i,j) + DV(i+1,j) )
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO NOWAIT

      if (vel2d_adv_hor .eq. J7) then
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
         do j=jmin-HALO,jmax+HALO
#endif
            do i=imin-HALO,imax+HALO-1
!              Note (KK): [U|V]adv defined on T-points (future U-points)
!                         Dadv defined on X-points (future V-points)
!                         note that Dadv is shifted to j+1 !!!
               Uadv(i,j) = _HALF_*( U(i,j)*DYU + U(i+1,j)*DYUIP1 )
               if (j .ne. jmin-HALO) then
                  if (az(i+1,j) .eq. 1) then
                     Vadv(i,j) = _HALF_*( V(i+1,j-1)*DXVPM + V(i+1,j)*DXVIP1 )
                  else if(az(i+1,j) .eq. 2) then
!                    Note (KK): can be included into case above, when
!                               V is properly mirrored across n/s open bdys
                     if (av(i+1,j) .eq. 2) then ! southern open bdy
                        Vadv(i,j) = V(i+1,j)*_HALF_*( DXVPM + DXVIP1 )
                     else if (av(i+1,j-1) .eq. 2) then ! northern open bdy
                        Vadv(i,j) = V(i+1,j-1)*_HALF_*( DXVPM + DXVIP1 )
                     end if
                  end if
                  if (ax(i,j) .eq. 0) then
                     if (au(i,j-1) .ne. 0) then
                        Dadv(i,j) = DU(i,j-1)
                     else if (au(i,j) .ne. 0) then
                        Dadv(i,j) = DU(i,j)
                     end if
                  else
                     Dadv(i,j) = _HALF_*( DU(i,j-1) + DU(i,j) )
                  end if
               end if
            end do
#ifndef SLICE_MODEL
         end do
#endif
!$OMP END DO

#ifdef SLICE_MODEL
!$OMP SINGLE
         Dadv(imin-HALO:imax+HALO-1,j+1) = Dadv(imin-HALO:imax+HALO-1,j)
         Uadv(imin-HALO:imax+HALO-1,j+1) = Uadv(imin-HALO:imax+HALO-1,j)
         Vadv(imin-HALO:imax+HALO-1,j+1) = Vadv(imin-HALO:imax+HALO-1,j)
!$OMP END SINGLE
#endif

#ifndef SLICE_MODEL
!        OMP-NOTE (KK): j loop must not be changed and cannot be threaded!
         do j=jmin-HALO,jmax+HALO-1
#endif
!$OMP DO SCHEDULE(RUNTIME)
            do i=imin-HALO,imax+HALO-1
!              Note (KK): [U|V]adv defined on V-points (future X-points)
!                         no change of Uadv for northern closed bdy
!                         Dadv defined on U-points (future T-points)
!                         note that Dadv is not shifted anymore !
!              KK-TODO: Vadv for western/eastern open bdys ?
               if (az(i+1,j) .eq. 0) then ! southern closed bdy
                  Uadv(i,j) = Uadv(i,j+1)
               else if (az(i+1,j+1) .ne. 0) then
                  Uadv(i,j) = _HALF_*( Uadv(i,j) + Uadv(i,j+1) )
                  Vadv(i,j) = _HALF_*( Vadv(i,j) + Vadv(i,j+1) )
               end if
               Dadv(i,j) = _HALF_*( Dadv(i,j) + Dadv(i,j+1) )
            end do
!$OMP END DO
#ifndef SLICE_MODEL
         end do
#endif

!$OMP SINGLE
         pDadv => Dadv
!$OMP END SINGLE

      else

!$OMP SINGLE
         pDadv => DU
!$OMP END SINGLE

      end if

!$OMP SINGLE
      if (vel2d_adv_hor.ne.UPSTREAM .and. vel2d_adv_hor.ne.J7) then
!        we need to update fadv(imax+HALO,jmin-HALO:jmax+HALO)
         call tic(TIM_UVADVH)
         call update_2d_halo(fadv,fadv,au,imin,jmin,imax,jmax,U_TAG)
         call wait_halo(U_TAG)
         call toc(TIM_UVADVH)
      end if
!$OMP END SINGLE

#ifdef _NUMERICAL_ANALYSES_OLD_
      if (do_numerical_analyses_2d) then
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
         do j=jmin-HALO,jmax+HALO  ! calculate square of u-velocity before advection step
#endif
            do i=imin-HALO,imax+HALO
               work2d(i,j) = fadv(i,j)**2
            end do
#ifndef SLICE_MODEL
         end do
#endif
!$OMP END DO
      end if
#endif

!$OMP SINGLE
      call do_advection(dtm,fadv,Uadv,Vadv,DUadv,DVadv,pDadv,pDadv, &
                        vel2d_adv_split,vel2d_adv_hor,_ZERO_,U_TAG, &
                        Dires=Dires,advres=UEx,nvd=p_nvd)
!$OMP END SINGLE

      if (calc_numdis) then

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
         do j=jmin,jmax
#endif
            do i=imin,imax
               nvd(i,j) = _HALF_*nvd(i,j)*Dires(i,j)/ARUD1
            end do
#ifndef SLICE_MODEL
         end do
#endif
!$OMP END DO

!$OMP SINGLE
         call update_2d_halo(nvd,nvd,au,imin,jmin,imax,jmax,U_TAG)
         call wait_halo(U_TAG)
!$OMP END SINGLE

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
         do j=jmin,jmax
#endif
            do i=imin,imax
               if (az(i,j) .eq. 1) then
                  numdis(i,j) = _HALF_*( nvd(i-1,j) + nvd(i,j) )/Dnew(i,j)*ARCD1
               end if
            end do
#ifndef SLICE_MODEL
         end do
#endif
!$OMP END DO
      end if

#ifdef _NUMERICAL_ANALYSES_OLD_
      if (do_numerical_analyses_2d) then

!$OMP SINGLE
         call do_advection(dtm,work2d,Uadv,Vadv,DUadv,DVadv,pDadv,pDadv, &
                           vel2d_adv_split,vel2d_adv_hor,_ZERO_,U_TAG)
!$OMP END SINGLE
 
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
         do j=jmin,jmax ! calculate kinetic energy dissipaion rate for u-velocity
#endif
            do i=imin,imax
               work2d(i,j) = _HALF_*( work2d(i,j) - fadv(i,j)**2 )/dtm*Dires(i,j)/ARUD1
            end do
#ifndef SLICE_MODEL
         end do
#endif
!$OMP END DO

!$OMP SINGLE
         call update_2d_halo(work2d,work2d,au,imin,jmin,imax,jmax,U_TAG)
         call wait_halo(U_TAG)
!$OMP END SINGLE

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
         do j=jmin,jmax
#endif
            do i=imin,imax
               if (az(i,j) .eq. 1) then
                  numdis_2d_old(i,j) = _HALF_*( work2d(i-1,j) + work2d(i,j) )/Dnew(i,j)*ARCD1
               end if
            end do
#ifndef SLICE_MODEL
         end do
#endif
!$OMP END DO

      end if
#endif

!     Here begins dimensional split advection for v-velocity

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
     do j=jmin-HALO,jmax+HALO-1
#endif
         do i=imin-HALO,imax+HALO
!           the velocity to be transported
            fadv(i,j) = V(i,j)/DV(i,j)
            if (vel2d_adv_hor .ne. J7) then
!              Note (KK): Uadv defined on X-points (future U-points)
!                         Vadv defined on T-points (future V-points)
               Uadv(i,j) = _HALF_*( U(i,j) + U(i,j+1) )
               Vadv(i,j) = _HALF_*( V(i,j) + V(i,j+1) )
            end if
!           Note (KK): DU only valid until imax+1
!                      therefore DUadv only valid until imax+1
            DUadv(i,j) = _HALF_*( DU(i,j) + DU(i,j+1) )
            DVadv(i,j) = _HALF_*( Dold(i,j+1) + Dnew(i,j+1) )
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO NOWAIT

      if (vel2d_adv_hor .eq. J7) then
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
         do j=jmin-HALO,jmax+HALO-1
#endif
            do i=imin-HALO,imax+HALO
!              Note (KK): [U|V]adv defined on T-points (future V-points)
!                         Dadv defined on X-points (future U-points)
!                         note that Dadv is shifted to i+1 !!!
               if (i .ne. imin-HALO) then
                  if (az(i,j+1) .eq. 1) then
                     Uadv(i,j) = _HALF_*( U(i-1,j+1)*DYUMP + U(i,j+1)*DYUJP1 )
                  else if(az(i,j+1) .eq. 2) then
!                    Note (KK): can be included into case above, when
!                               U is properly mirrored across w/e open bdys
                     if (au(i,j+1) .eq. 2) then ! western open bdy
                        Uadv(i,j) = U(i,j+1)*_HALF_*( DYUMP + DYUJP1 )
                     else if (au(i-1,j+1) .eq. 2) then ! eastern open bdy
                        Uadv(i,j) = U(i-1,j+1)*_HALF_*( DYUMP + DYUJP1 )
                     end if
                  end if
                  if (ax(i,j) .eq. 0) then
                     if (av(i-1,j) .ne. 0) then
                        Dadv(i,j) = DV(i-1,j)
                     else if (av(i,j) .ne. 0) then
                        Dadv(i,j) = DV(i,j)
                     end if
                  else
                     Dadv(i,j) = _HALF_*( DV(i-1,j) + DV(i,j) )
                  end if
               end if
               Vadv(i,j) = _HALF_*( V(i,j)*DXV + V(i,j+1)*DXVJP1 )
            end do
#ifdef SLICE_MODEL
!$OMP END DO
#endif

!           OMP-NOTE (KK): i loop must not be changed and cannot be threaded!
            do i=imin-HALO,imax+HALO-1
!              Note (KK): [U|V]adv defined on U-points (future X-points)
!                         no change of Vadv for eastern closed bdy
!                         Dadv defined on V-points (future T-points)
!                         note that Dadv is not shifted anymore !
!              KK-TODO: Uadv for northern/southern open bdys ?
               if (az(i,j+1) .eq. 0) then ! western closed bdy
                  Vadv(i,j) = Vadv(i+1,j)
               else if (az(i+1,j+1) .ne. 0) then
                  Uadv(i,j) = _HALF_*( Uadv(i,j) + Uadv(i+1,j) )
                  Vadv(i,j) = _HALF_*( Vadv(i,j) + Vadv(i+1,j) )
               end if
               Dadv(i,j) = _HALF_*( Dadv(i,j) + Dadv(i+1,j) )
            end do
#ifndef SLICE_MODEL
         end do
!$OMP END DO
#endif

!$OMP SINGLE
         pDadv => Dadv
!$OMP END SINGLE

      else

!$OMP SINGLE
         pDadv => DV
!$OMP END SINGLE

      end if

!$OMP SINGLE
      if (vel2d_adv_hor.ne.UPSTREAM .and. vel2d_adv_hor.ne.J7) then
!        we need to update fadv(imin-HALO:imax+HALO,jmax+HALO)
         call tic(TIM_UVADVH)
         call update_2d_halo(fadv,fadv,av,imin,jmin,imax,jmax,V_TAG)
         call wait_halo(V_TAG)
         call toc(TIM_UVADVH)
      end if
!$OMP END SINGLE

#ifdef _NUMERICAL_ANALYSES_OLD_
      if (do_numerical_analyses_2d) then
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
         do j=jmin-HALO,jmax+HALO  ! calculate square of v-velocity before advection step
#endif
            do i=imin-HALO,imax+HALO
               work2d(i,j) = fadv(i,j)**2
            end do
#ifndef SLICE_MODEL
         end do
#endif
!$OMP END DO
      end if
#endif

!$OMP SINGLE
      call do_advection(dtm,fadv,Uadv,Vadv,DUadv,DVadv,pDadv,pDadv, &
                        vel2d_adv_split,vel2d_adv_hor,_ZERO_,V_TAG, &
                        Dires=Dires,advres=VEx,nvd=p_nvd)
!$OMP END SINGLE

      if (calc_numdis) then

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
         do j=jmin,jmax
#endif
            do i=imin,imax
               nvd(i,j) = _HALF_*nvd(i,j)*Dires(i,j)/ARVD1
            end do
#ifndef SLICE_MODEL
         end do
#endif
!$OMP END DO

!$OMP SINGLE
         call update_2d_halo(nvd,nvd,av,imin,jmin,imax,jmax,V_TAG)
         call wait_halo(V_TAG)
!$OMP END SINGLE

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
         do j=jmin,jmax
#endif
            do i=imin,imax
               if (az(i,j) .eq. 1) then
                  numdis(i,j) = numdis(i,j)                                      &
                               +_HALF_*( nvd(i,j-1) + nvd(i,j) )/Dnew(i,j)*ARCD1
               end if
            end do
#ifndef SLICE_MODEL
         end do
#endif
!$OMP END DO
      end if

#ifdef _NUMERICAL_ANALYSES_OLD_
      if (do_numerical_analyses_2d) then

!$OMP SINGLE
         call do_advection(dtm,work2d,Uadv,Vadv,DUadv,DVadv,pDadv,pDadv, &
                           vel2d_adv_split,vel2d_adv_hor,_ZERO_,V_TAG)
!$OMP END SINGLE

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
         do j=jmin,jmax ! calculate kinetic energy dissipaion rate for v-velocity
#endif
            do i=imin,imax
               work2d(i,j) = _HALF_*( work2d(i,j) - fadv(i,j)**2 )/dtm*Dires(i,j)/ARVD1
            end do
#ifndef SLICE_MODEL
         end do
#endif
!$OMP END DO

!$OMP SINGLE
#ifdef SLICE_MODEL
         work2d(imin:imax,j-1) = work2d(imin:imax,j)
#else
         call update_2d_halo(work2d,work2d,av,imin,jmin,imax,jmax,V_TAG)
         call wait_halo(V_TAG)
#endif
!$OMP END SINGLE

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
         do j=jmin,jmax
#endif
            do i=imin,imax
               if (az(i,j) .eq. 1) then
                  numdis_2d_old(i,j) = numdis_2d_old(i,j)                    &
                     +_HALF_*( work2d(i,j-1) + work2d(i,j) )/Dnew(i,j)*ARCD1
               end if
            end do
#ifndef SLICE_MODEL
         end do
#endif
!$OMP END DO

      end if
#endif

!$OMP END PARALLEL

   else ! if (vel2d_adv_hor .eq. NOADV)

      UEx = _ZERO_ ; VEx = _ZERO_

   end if

   call toc(TIM_UVADV)
#ifdef DEBUG
   write(debug,*) 'Leaving uv_advect()'
   write(debug,*)
#endif
   return
   end subroutine uv_advect
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
