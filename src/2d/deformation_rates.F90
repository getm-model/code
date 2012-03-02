#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: deformation_rates - \label{deformation_rates}
!
! !INTERFACE:
   subroutine deformation_rates(U,V,DU,DV,dudxC,dudxV,dudxU,         &
                                          dvdyC,dvdyU,dvdyV,         &
                                          dudyX,dvdxX,shearX,        &
                                          dvdxU,shearU,dudyV,shearV)
                                          
!  Note (KK): keep in sync with interface in m2d_general.F90
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,az,au,av,ax
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxc,dyc,dxu,dyu,dxv,dyv,dxx,dyx,arcd1,arud1,arvd1
#else
   use domain, only: dx,dy,ard1
#endif
   use domain, only: NWB,NNB,NEB,NSB
   use domain, only: wi,wfj,wlj,nj,nfi,nli,ei,efj,elj,sj,sfi,sli
   use getm_timers, only: tic,toc,TIM_DEFORM

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(in)           :: U,V,DU,DV
!
! !OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(out),optional :: dudxC,dudxV,dudxU
   REALTYPE,dimension(E2DFIELD),intent(out),optional :: dvdyC,dvdyU,dvdyV
   REALTYPE,dimension(E2DFIELD),intent(out),optional :: dudyX,dvdxX,shearX
   REALTYPE,dimension(E2DFIELD),intent(out),optional :: dvdxU,shearU
   REALTYPE,dimension(E2DFIELD),intent(out),optional :: dudyV,shearV
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
!  allocated outside a module, therefore saved
   REALTYPE,dimension(E2DFIELD)             :: work2d
   REALTYPE,dimension(:,:),allocatable,save :: u_vel,v_vel
   integer                                  :: rc
   logical,save                             :: first=.true.
   REALTYPE                                 :: dxdy,dydx,tmp,velgrad
   integer                                  :: i,j,n
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'deformation_rates() # ',Ncall
#endif
#ifdef SLICE_MODEL
   j=jmax/2
#endif
   call tic(TIM_DEFORM)

   if (first) then
      allocate(u_vel(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (u_vel)'
      u_vel=_ZERO_

      allocate(v_vel(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (v_vel)'
      v_vel=_ZERO_

      first = .false.
   end if

!  Note (KK): the quantities calculated here are elements of the velocity gradient
!             obtained by covariant differentiation of the velocity vector
!             for exact 3D deformation tensor in horizontal orthogonal coordinates
!             see e.g. Kantha and Clayson 2000, page 35

#ifndef SLICE_MODEL
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,tmp,velgrad)
#endif

!  zonal velocity
#ifndef SLICE_MODEL
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO
#endif
      do i=imin-HALO,imax+HALO-1
         if (au(i,j) .ge. 1) then
            u_vel(i,j) = U(i,j)/DU(i,j)
         end if
      end do
#ifndef SLICE_MODEL
   end do
!$OMP END DO NOWAIT
#else
   u_vel(imin-HALO:imax+HALO-1,j+1) = u_vel(imin-HALO:imax+HALO-1,j)
#endif

!  meridional velocity
#ifndef SLICE_MODEL
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO-1
#endif
      do i=imin-HALO,imax+HALO
         if (av(i,j) .ge. 1) then
            v_vel(i,j) = V(i,j)/DV(i,j)
         end if
      end do
#ifndef SLICE_MODEL
   end do
!$OMP END DO
#else
   v_vel(:,j-1) = v_vel(:,j)
   v_vel(:,j+1) = v_vel(:,j)
#endif

#ifdef _CORRECT_METRICS_
#if defined(SPHERICAL) || defined(CURVILINEAR)
#ifndef SLICE_MODEL
!$OMP MASTER
#endif
!  mirror velocities based on ouflow condition for transverse velocity
   do n = 1,NNB ! northern open boundaries (dudyX=0)
      j = nj(n)-1
      do i = nfi(n)-HALO,nli(n)+HALO-1
         if (au(i,j+1) .eq. 3) then
            dydx = ( DYVIP1 - DYV ) / DXX
            if (av(i,j) .eq. 3) then ! concave W\N
               dxdy = ( DXUJP1 - DXU ) / DYX
               u_vel(i,j+1) = (                                       &
                                 (_ONE_-_QUART_*dydx*dxdy)*u_vel(i,j) &
                               + dydx*v_vel(i+1,j)                    &
                              )                                       &
                              / (_ONE_+_QUART_*dydx*dxdy)
            else if (av(i+1,j) .eq. 3) then ! concave N/E
               dxdy = ( DXUJP1 - DXU ) / DYX
               u_vel(i,j+1) = (                                       &
                                 (_ONE_+_QUART_*dydx*dxdy)*u_vel(i,j) &
                               + dydx*v_vel(i,j)                      &
                              )                                       &
                              / (_ONE_-_QUART_*dydx*dxdy)
            else
               u_vel(i,j+1) = u_vel(i,j) + _HALF_*(v_vel(i,j)+v_vel(i+1,j))*dydx
            end if
         end if
      end do
   end do
   do n = 1,NSB ! southern open boundaries (dudyX=0)
      j = sj(n)
      do i = sfi(n)-HALO,sli(n)+HALO-1
         if (au(i,j) .eq. 3) then
            dydx = ( DYVIP1 - DYV ) / DXX
            if (av(i,j) .eq. 3) then ! concave W/S
               dxdy = ( DXUJP1 - DXU ) / DYX
               u_vel(i,j) = (                                         &
                               (_ONE_+_QUART_*dydx*dxdy)*u_vel(i,j+1) &
                             - dydx*v_vel(i+1,j)                      &
                            )                                         &
                            / (_ONE_-_QUART_*dydx*dxdy)
            else if (av(i+1,j) .eq. 3) then ! concave S\E
               dxdy = ( DXUJP1 - DXU ) / DYX
               u_vel(i,j) = (                                         &
                               (_ONE_-_QUART_*dydx*dxdy)*u_vel(i,j+1) &
                             - dydx*v_vel(i,j)                        &
                            )                                         &
                            / (_ONE_+_QUART_*dydx*dxdy)
            else
               u_vel(i,j) = u_vel(i,j+1) - _HALF_*(v_vel(i,j)+v_vel(i+1,j))*dydx
            end if
         end if
      end do
   end do
   do n = 1,NWB ! western open boundaries (dvdxX=0)
      i = wi(n)
      do j = wfj(n)-HALO,wlj(n)+HALO-1
         if (av(i,j) .eq. 3) then
            dxdy = ( DXUJP1 - DXU ) / DYX
            v_vel(i,j) = v_vel(i+1,j) - _HALF_*(u_vel(i,j)+u_vel(i,j+1))*dxdy
         end if
      end do
   end do
   do n = 1,NEB ! eastern open boundaries (dvdxX=0)
      i = ei(n)-1
      do j = efj(n)-HALO,elj(n)+HALO-1
         if (av(i+1,j) .eq. 3) then
            dxdy = ( DXUJP1 - DXU ) / DYX
            v_vel(i+1,j) = v_vel(i,j) + _HALF_*(u_vel(i,j)+u_vel(i,j+1))*dxdy
         end if
      end do
   end do
#ifndef SLICE_MODEL
!$OMP END MASTER
!$OMP BARRIER
#endif
#endif
#endif


!!!!!!!!!!!!!!!!!!!!!!!
!!!zonal strain rate!!!
!!!!!!!!!!!!!!!!!!!!!!!

   if (present(dudxC)) then
!     zonal strain rate at T-points
#ifndef SLICE_MODEL
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO+1,jmax+HALO-1
#endif
         do i=imin-HALO+1,imax+HALO-1
            if (az(i,j) .eq. 1) then
!              Note (KK): since U(au=0)=0, dudxC can also be calculated near closed boundaries
               dudxC(i,j) = (u_vel(i,j) - u_vel(i-1,j)) / DXC
#ifdef _CORRECT_METRICS_
#if defined(SPHERICAL) || defined(CURVILINEAR)
               dudxC(i,j) = dudxC(i,j) + _HALF_*(v_vel(i,j)+v_vel(i,j-1))*(DXV-DXVJM1)*ARCD1
#endif
#endif
            else if (az(i,j) .eq. 2) then
!              Note (KK): outflow condition dudxC=0 in W/E open boundaries
!                         v_vel outside N/S open boundary from outflow condition dvdyC=0
               if (av(i,j-1) .eq. 2) then ! northern open boundary
                  dudxC(i,j) = (u_vel(i,j) - u_vel(i-1,j)) / DXC
#ifdef _CORRECT_METRICS_
#if defined(SPHERICAL) || defined(CURVILINEAR)
                  tmp = v_vel(i,j-1) - _HALF_*(u_vel(i,j)+u_vel(i-1,j))*(DYU-DYUIM1)/DXC
                  dudxC(i,j) = dudxC(i,j) + _HALF_*(tmp+v_vel(i,j-1))*(DXV-DXVJM1)*ARCD1
#endif
#endif
               else if (av(i,j) .eq. 2) then ! southern open boundary
                  dudxC(i,j) = (u_vel(i,j) - u_vel(i-1,j)) / DXC
#ifdef _CORRECT_METRICS_
#if defined(SPHERICAL) || defined(CURVILINEAR)
                  tmp = v_vel(i,j) + _HALF_*(u_vel(i,j)+u_vel(i-1,j))*(DYU-DYUIM1)/DXC
                  dudxC(i,j) = dudxC(i,j) + _HALF_*(v_vel(i,j)+tmp)*(DXV-DXVJM1)*ARCD1
#endif
#endif
               end if
            end if
         end do
#ifndef SLICE_MODEL
      end do
!$OMP END DO NOWAIT
#else
      dudxC(imin-HALO+1:imax+HALO-1,j+1) = dudxC(imin-HALO+1:imax+HALO-1,j)
#endif

      if (present(dudxU)) then
   !     interpolation of zonal strain rate to U-points
#ifndef SLICE_MODEL
!$OMP BARRIER
!$OMP DO SCHEDULE(RUNTIME)
         do j=jmin-HALO+1,jmax+HALO-1
#endif
            do i=imin-HALO+1,imax+HALO-2
!              KK-TODO: outflow condition dudxU(au=2)=0 ?
               if (au(i,j).eq.1 .or. au(i,j).eq.3) then
                  dudxU(i,j) = _HALF_*(dudxC(i,j) + dudxC(i+1,j))
               else if (au(i,j) .eq. 0) then
!                 KK-TODO: is this correct (U(au=0)=0 + slip condition) ?
!                 Note (KK): in case of no-slip dudxU=0
                  if (az(i  ,j) .eq. 1) dudxU(i,j) = dudxC(i  ,j)
                  if (az(i+1,j) .eq. 1) dudxU(i,j) = dudxC(i+1,j)
               end if
            end do
#ifndef SLICE_MODEL
         end do
!$OMP END DO NOWAIT
#else
         dudxU(imin-HALO+1:imax+HALO-2,j+1) = dudxU(imin-HALO+1:imax+HALO-2,j)
#endif
      end if

   end if


   if (present(dudxV)) then
!     zonal strain rate at V-points
#ifndef SLICE_MODEL
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO,jmax+HALO-1
#endif
!        calculate u_velX
         do i=imin-HALO,imax+HALO-1
            if (ax(i,j) .ne. 0) then
               work2d(i,j) = _HALF_ * ( u_vel(i,j) + u_vel(i,j+1) )
            else
               work2d(i,j) = _ZERO_
            end if
         end do
!        set dudxV
         do i=imin-HALO+1,imax+HALO-1
!           Note (KK): outflow condition dudxV(av=3)=0
            if (av(i,j).eq.1 .or. av(i,j).eq.2) then
               dudxV(i,j) = ( work2d(i,j) - work2d(i-1,j) )/DXV
#ifdef _CORRECT_METRICS_
#if defined(SPHERICAL) || defined(CURVILINEAR)
               dudxV(i,j) = dudxV(i,j) + v_vel(i,j)*(DXCJP1-DXC)*ARVD1
#endif
#endif
            else if (av(i,j) .eq. 0) then
!              Note (KK): V(av=0)=0 and slip condition dudyV(av=0)=0 at N/S closed bdys
               if (az(i,j  ).eq.1 .or. az(i,j  ).eq.2) dudxV(i,j) = (u_vel(i,j  ) - u_vel(i-1,j  ))/DXV
               if (az(i,j+1).eq.1 .or. az(i,j+1).eq.2) dudxV(i,j) = (u_vel(i,j+1) - u_vel(i-1,j+1))/DXV
            end if
         end do
#ifndef SLICE_MODEL
      end do
!$OMP END DO NOWAIT
#else
      dudxV(imin-HALO+1:imax+HALO-1,j-1) = dudxV(imin-HALO+1:imax+HALO-1,j)
      dudxV(imin-HALO+1:imax+HALO-1,j+1) = dudxV(imin-HALO+1:imax+HALO-1,j)
#endif
   end if

#ifndef SLICE_MODEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!meridional strain rate!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (present(dvdyC)) then
!     meridional strain rate at T-points
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO+1,jmax+HALO-1
         do i=imin-HALO+1,imax+HALO-1
            if (az(i,j) .eq. 1) then
!              Note (KK): since V(av=0)=0, dvdyC can also be calculated near closed boundaries
               dvdyC(i,j) = (v_vel(i,j) - v_vel(i,j-1)) / DYC
#ifdef _CORRECT_METRICS_
#if defined(SPHERICAL) || defined(CURVILINEAR)
               dvdyC(i,j) = dvdyC(i,j) + _HALF_*(u_vel(i,j)+u_vel(i-1,j))*(DYU-DYUIM1)*ARCD1
#endif
#endif
            else if(az(i,j) .eq. 2) then
!              Note (KK): outflow condition dvdyC=0 in N/S open boundaries
!                         u_vel outside W/E open boundary from outflow condition dudxC=0
!                         dvdyC in W/E open boundary needed for dvdyV(av=3)
               if (au(i-1,j) .eq. 2) then ! eastern open boundary
                  dvdyC(i,j) = (v_vel(i,j) - v_vel(i,j-1)) / DYC
#ifdef _CORRECT_METRICS_
#if defined(SPHERICAL) || defined(CURVILINEAR)
                  tmp = u_vel(i-1,j) - _HALF_*(v_vel(i,j)+v_vel(i,j-1))*(DXV-DXVJM1)/DYC
                  dvdyC(i,j) = dvdyC(i,j) + _HALF_*(tmp+u_vel(i-1,j))*(DYU-DYUIM1)*ARCD1
#endif
#endif
               else if (au(i,j) .eq. 2) then ! western open boundary
                  dvdyC(i,j) = (v_vel(i,j) - v_vel(i,j-1)) / DYC
#ifdef _CORRECT_METRICS_
#if defined(SPHERICAL) || defined(CURVILINEAR)
                  tmp = u_vel(i,j) + _HALF_*(v_vel(i,j)+v_vel(i,j-1))*(DXV-DXVJM1)/DYC
                  dvdyC(i,j) = dvdyC(i,j) + _HALF_*(u_vel(i,j)+tmp)*(DYU-DYUIM1)*ARCD1
#endif
#endif
               end if
            end if
         end do
      end do
!$OMP END DO NOWAIT

      if (present(dvdyV)) then
!        interpolation of meridional strain rate to V-points
!$OMP BARRIER
!$OMP DO SCHEDULE(RUNTIME)
         do j=jmin-HALO+1,jmax+HALO-2
            do i=imin-HALO+1,imax+HALO-1
!              KK-TODO: outflow condition dvdyV(av=2)=0 ?
               if (av(i,j).eq.1 .or. av(i,j).eq.3) then
!                 KK-TODO: for dvdyV(av=3) average of nonlinear metric
!                          correction might cause violation of outflow
!                          condition dvdx=0?
                  dvdyV(i,j) = _HALF_*(dvdyC(i,j) + dvdyC(i,j+1))
               else if (av(i,j) .eq. 0) then
!                 KK-TODO: is this correct (V(av=0)=0 + slip condition) ?
!                 Note (KK): in case of no-slip dvdyV=0
                  if (az(i,j  ) .eq. 1) dvdyV(i,j) = dvdyC(i,j  )
                  if (az(i,j+1) .eq. 1) dvdyV(i,j) = dvdyC(i,j+1)
               end if
            end do
         end do
!$OMP END DO NOWAIT
      end if

   end if

   if (present(dvdyU)) then
!     meriodional strain rate at U-points
!     calculate v_velX
!$OMP BARRIER
!     Note (KK): we might need omp barrier because work2d was used before with nowait
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO,jmax+HALO-1
         do i=imin-HALO,imax+HALO-1
            if (ax(i,j) .ne. 0) then
               work2d(i,j) = _HALF_ * ( v_vel(i,j) + v_vel(i+1,j) )
            else
               work2d(i,j) = _ZERO_
            end if
         end do
      end do
!$OMP END DO
!     set dvdyU
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO+1,jmax+HALO-1
         do i=imin-HALO,imax+HALO-1
!           Note (KK): outflow condition dvdyU(au=3)=0
            if (au(i,j).eq.1 .or. au(i,j).eq.2) then
               dvdyU(i,j) = ( work2d(i,j) - work2d(i,j-1) ) / DYU
#ifdef _CORRECT_METRICS_
#if defined(SPHERICAL) || defined(CURVILINEAR)
               dvdyU(i,j) = dvdyU(i,j) + u_vel(i,j)*(DYCIP1-DYC)*ARUD1
#endif
#endif
            else if (au(i,j) .eq. 0) then
!              Note (KK): U(au=0)=0 and slip condition dvdxU(au=0)=0
               if (az(i  ,j) .ge. 1) dvdyU(i,j) = (v_vel(i  ,j) - v_vel(i  ,j-1))/DYU
               if (az(i+1,j) .ge. 1) dvdyU(i,j) = (v_vel(i+1,j) - v_vel(i+1,j-1))/DYU
            end if
         end do
      end do
!$OMP END DO NOWAIT
   end if

#endif


!!!!!!!!!!!!!!!!
!!!shear rate!!!
!!!!!!!!!!!!!!!!

!  shear rate at X-points

   if (present(dvdxX) .or. present(dudyX) .or. present(shearX)) then

#ifndef SLICE_MODEL
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO,jmax+HALO-1
#endif
         do i=imin-HALO,imax+HALO-1
!           Note (KK): at closed boundaries slip condition for tangential
!                      velocity in terms of covariant derivative
!                      no normal flow through closed boundaries does only imply
!                      vanishing partial derivative!!!
!                      at concave and at convex corners shearX(ax=0) is set to zero
!                      (despite discontinuous velocities for the latter!?)
            if (ax(i,j) .eq. 1) then
               tmp = _ZERO_
               if (au(i,j).eq.1 .or. au(i,j+1).eq.1) then
!                 Note (KK): excludes concave and W/E open boundaries (dvdxX=0)
!                            includes convex open boundaries (no mirroring)
                  velgrad = (v_vel(i+1,j) - v_vel(i,j)) / DXX
                  if (present(dvdxX)) then
                     dvdxX(i,j) = velgrad
#ifdef _CORRECT_METRICS_
#if defined(SPHERICAL) || defined(CURVILINEAR)
                     dvdxX(i,j) = dvdxX(i,j) - _HALF_*(u_vel(i,j)+u_vel(i,j+1))/DXX*(DXUJP1-DXU)/DYX
#endif
#endif
                  end if
                  tmp = velgrad
               end if

#ifndef SLICE_MODEL
               if (av(i,j).eq.1 .or. av(i+1,j).eq.1) then
!                 Note (KK): excludes concave and N/S open boundaries (dudyX=0)
!                            includes convex open boundaries (no mirroring)
                  velgrad = (u_vel(i,j+1) - u_vel(i,j)) / DYX
                  if (present(dudyX)) then
                     dudyX(i,j) = velgrad
#ifdef _CORRECT_METRICS_
#if defined(SPHERICAL) || defined(CURVILINEAR)
                     dudyX(i,j) = dudyX(i,j) - _HALF_*(v_vel(i,j)+v_vel(i+1,j))/DYX*(DYVIP1-DYV)/DXX
#endif
#endif
                  end if
                  tmp = tmp + velgrad
               end if
#endif

               if (present(shearX)) then
#if !(defined(_CORRECT_METRICS_) && (defined(SPHERICAL) || defined(CURVILINEAR)))
                  shearX(i,j) = tmp
#else
!                 Note (KK): although sum of dvdxX and dudyX would give shearX
!                            a combined formula that does not involve averages
!                            is used
                  shearX(i,j) = DYX * (v_vel(i+1,j)/DYVIP1 - v_vel(i,j)/DYV) / DXX
#ifndef SLICE_MODEL
                  shearX(i,j) = shearX(i,j) + DXX * (u_vel(i,j+1)/DXUJP1 - u_vel(i,j)/DXU) / DYX
#endif
#endif
               end if

#ifdef _CORRECT_METRICS_
#if defined(SPHERICAL) || defined(CURVILINEAR)
            else

               if (av(i,j).eq.0 .and. av(i+1,j).eq.0) then
!                 Note (KK): exclude convex corners
!                            at W/E closed boundaries slip condition dvdx=0
                  if (au(i,j) .eq. 1) then ! northern closed boundary
                     velgrad = - u_vel(i,j)/DXX*(DXUJP1-DXU)/DYX
                     if (present(dvdxX)) dvdxX(i,j) = velgrad
                     if (present(shearX)) shearX(i,j) = velgrad
                  end if
                  if (au(i,j+1) .eq. 1) then ! southern closed boundary
                     velgrad = - u_vel(i,j+1)/DXX*(DXUJP1-DXU)/DYX
                     if (present(dvdxX)) dvdxX(i,j) = velgrad
                     if (present(shearX)) shearX(i,j) = velgrad
                  end if
               end if

#ifndef SLICE_MODEL
               if (au(i,j).eq.0 .and. au(i,j+1).eq.0) then
!                 Note (KK): exclude convex corners
!                            at N/S closed boundaries slip condition dudy=0
                  if (av(i,j) .eq. 1) then ! eastern closed boundary
                     velgrad = - v_vel(i,j)/DYX*(DYVIP1-DYV)/DXX
                     if (present(dudyX)) dudyX(i,j) = velgrad
                     if (present(shearX)) shearX(i,j) = shearX(i,j) + velgrad
                  end if
                  if (av(i+1,j) .eq. 1) then ! western closed boundary
                     velgrad = - v_vel(i+1,j)/DYX*(DYVIP1-DYV)/DXX
                     if (present(dudyX)) dudyX(i,j) = velgrad
                     if (present(shearX)) shearX(i,j) = shearX(i,j) + velgrad
                  end if
               end if
#endif

#endif
#endif

            end if
         end do
#ifndef SLICE_MODEL
      end do
!$OMP END DO
#else
      if (present(dvdxX)) then
         dvdxX(imin-HALO:imax+HALO-1,j-1) = dvdxX(imin-HALO:imax+HALO-1,j)
         dvdxX(imin-HALO:imax+HALO-1,j+1) = dvdxX(imin-HALO:imax+HALO-1,j)
      end if
      if (present(shearX)) then
         shearX(imin-HALO:imax+HALO-1,j-1) = shearX(imin-HALO:imax+HALO-1,j)
         shearX(imin-HALO:imax+HALO-1,j+1) = shearX(imin-HALO:imax+HALO-1,j)
      end if
#endif

   end if

!  shear rate at U-points
   if (present(dvdxU)) then
#ifndef SLICE_MODEL
!$OMP BARRIER
!     Note (KK): we need omp barrier here, because work2d was used with nowait
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO+1,jmax+HALO-1
#endif
!        calculate v_velC
         do i=imin-HALO,imax+HALO
!           Note (KK): we only need v_velC(az=1)
            if (az(i,j) .eq. 1) then
               work2d(i,j) = _HALF_ * ( v_vel(i,j-1) + v_vel(i,j) )
            end if
         end do
!        calculate dvdxU
         do i=imin-HALO,imax+HALO-1
!           Note (KK): slip condition dvdxU(au=0)=0
!                      outflow condition at W/E open bdys dvdxU(au=2)=0
!           KK-TODO: metric correction
            if (au(i,j) .eq. 1) then
               dvdxU(i,j) = ( work2d(i+1,j) - work2d(i,j) ) /DXU
            else if (au(i,j) .eq. 3) then
               if (au(i,j-1) .eq. 1) then ! northern open bdy
                  dvdxU(i,j) = (v_vel(i+1,j-1) - v_vel(i  ,j-1))/DXU
               else ! southern open bdy
                  dvdxU(i,j) = (v_vel(i+1,j  ) - v_vel(i  ,j  ))/DXU
               end if
            end if
         end do
#ifndef SLICE_MODEL
      end do
!$OMP END DO NOWAIT
#else
      dvdxU(imin-HALO:imax+HALO-1,j+1) = dvdxU(imin-HALO:imax+HALO-1,j)
#endif
   end if
   if (present(shearX) .and. present(shearU)) then
!     interpolation of shear rate to U-points
#ifndef SLICE_MODEL
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO+1,jmax+HALO-1
#endif
         do i=imin-HALO,imax+HALO-1
!           Note (KK): slip condition dvdxU(au=0)=0
!                      prolonged outflow condition dudyU(au=3)=0
!                      shearU(au=3) would require shearX outside open boundary
!                      (however shearU(au=3) not needed, therefore not calculated)
            if (az(i,j).eq.1 .or. az(i+1,j).eq.1) then
               shearU(i,j) = _HALF_*(shearX(i,j-1) + shearX(i,j))
            end if
         end do
#ifndef SLICE_MODEL
      end do
!$OMP END DO
#else
      shearU(imin-HALO:imax+HALO-1,j+1) = shearU(imin-HALO:imax+HALO-1,j)
#endif
   end if

!  shear rate at V-points
#ifndef SLICE_MODEL
   if (present(dudyV)) then
!     calculate u_velC
!$OMP BARRIER
!     Note (KK): we need omp barrier here, because work2d was used with nowait
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO+1,imax+HALO-1
!           Note (KK): we only need u_velC(az=1)
            if (az(i,j) .eq. 1) then
               work2d(i,j) = _HALF_ * ( u_vel(i-1,j) + u_vel(i,j) )
            end if
         end do
      end do
!$OMP END DO
!     calculate dudyV
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO,jmax+HALO-1
         do i=imin-HALO+1,imax+HALO-1
!           Note (KK): slip condition dudyV(av=0)=0
!                      outflow condition at N/S open bdys dudyV(av=2)=0
!           KK-TODO: metric correction
            if (av(i,j) .eq. 1) then
               dudyV(i,j) = ( work2d(i,j+1) - work2d(i,j) ) / DYV
            else if (av(i,j) .eq. 3) then
               if (av(i-1,j) .eq. 1) then ! eastern open bdy
                  dudyV(i,j) = (u_vel(i-1,j+1) - u_vel(i-1,j  ))/DYV
               else ! western open bdy
                  dudyV(i,j) = (u_vel(i  ,j+1) - u_vel(i  ,j  ))/DYV
               end if
            end if
         end do
      end do
!$OMP END DO NOWAIT
   end if
#endif
   if (present(shearX) .and. present(shearV)) then
!     interpolation of shear rate to V-points
#ifndef SLICE_MODEL
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO,jmax+HALO-1
#endif
         do i=imin-HALO+1,imax+HALO-1
!           Note (KK): slip condition dudyV(av=0)=0
!                      prolonged outflow condition dvdxV(av=3)=0
!                      shearV(av=3) would require shearX outside open boundary
!                      (however shearV(av=3) not needed, therefore not calculated)
            if (az(i,j).eq.1 .or. az(i,j+1).eq.1) then
               shearV(i,j) = _HALF_*(shearX(i-1,j) + shearX(i,j))
            end if
         end do
#ifndef SLICE_MODEL
      end do
!$OMP END DO NOWAIT
#else
      shearV(imin-HALO+1:imax+HALO-1,j-1) = shearV(imin-HALO+1:imax+HALO-1,j)
      shearV(imin-HALO+1:imax+HALO-1,j+1) = shearV(imin-HALO+1:imax+HALO-1,j)
#endif
   end if

#ifndef SLICE_MODEL
!$OMP END PARALLEL
#endif

   call toc(TIM_DEFORM)
#ifdef DEBUG
   write(debug,*) 'Leaving deformation_rates()'
   write(debug,*)
#endif
   return
   end subroutine deformation_rates

!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
