#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: tracer_diffusion - lateral tracer diffusion
! \label{sec-tracer-diffusion}
!
! !INTERFACE:
   subroutine tracer_diffusion(f,hn,AH_method,AH_const,AH_Prt,AH_stirr_const, &
                               phymix)
!  Note (KK): keep in sync with interface in getm_fabm.F90
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax,az,au,av,ax
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxu,dyu,dxv,dyv,arcd1,arud1,arvd1
#else
   use domain, only: dx,dy,ard1
#endif
   use variables_3d, only: dt
   use variables_3d, only: diffxx,diffxy,diffyx,diffyy
   use variables_les, only: AmU_3d,AmV_3d
   use variables_3d, only: do_numerical_analyses_3d
   use getm_timers, only: tic,toc,TIM_TRACEDIFF
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,intent(in)           :: hn(I3DFIELD)
   integer,intent(in)            :: AH_method
   REALTYPE,intent(in)           :: AH_const,AH_Prt,AH_stirr_const
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,intent(inout)        :: f(I3DFIELD)
!
! !OUTPUT PARAMETERS:
   REALTYPE,intent(out),optional :: phymix(I3DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
   REALTYPE,dimension(I2DFIELD)             :: work2d
!  must be saved because of allocation outside a module
   REALTYPE,dimension(:,:),allocatable,save :: delfxU,delfyV,delfxV,delfyU
   REALTYPE,dimension(:,:),allocatable,save :: phymixU,phymixV
   REALTYPE                                 :: dfdxU,dfdyV
   logical,save                             :: first=.true.
   logical                                  :: calc_phymix
   integer                                  :: rc
   integer                                  :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'tracer_diffusion() # ',Ncall
#endif
#ifdef SLICE_MODEL
   j = jmax/2 ! this MUST NOT be changed!!!
#endif
   call tic(TIM_TRACEDIFF)

   if (first) then
      allocate(delfxU(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'tracer_diffusion: Error allocating memory (delfxU)'
      delfxU=_ZERO_

#ifndef SLICE_MODEL
      allocate(delfyV(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'tracer_diffusion: Error allocating memory (delfyV)'
      delfyV=_ZERO_

      allocate(delfxV(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'tracer_diffusion: Error allocating memory (delfxV)'
      delfxV=_ZERO_

      allocate(delfyU(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'tracer_diffusion: Error allocating memory (delfyU)'
      delfyU=_ZERO_
#endif

      if (do_numerical_analyses_3d) then

         allocate(phymixU(I2DFIELD),stat=rc)
         if (rc /= 0) stop 'tracer_diffusion: Error allocating memory (phymixU)'
         phymixU=_ZERO_

#ifndef SLICE_MODEL
         allocate(phymixV(I2DFIELD),stat=rc)
         if (rc /= 0) stop 'tracer_diffusion: Error allocating memory (phymixV)'
         phymixV=_ZERO_
#endif

      end if

      first=.false.
   end if

   calc_phymix = (do_numerical_analyses_3d .and. present(phymix))

!  Note (KK): diffusion only within new layer heigts

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j)                                         &
!$OMP          PRIVATE(i,k)

   do k=1,kmax

!     x-change at U-points
!$OMP DO SCHEDULE (RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin-HALO,jmax+HALO
#endif
         do i=imin-HALO,imax+HALO-1
            if (au(i,j) .ge. 1) then ! zero gradient across closed boundary
               delfxU(i,j) = f(i+1,j,k) - f(i,j,k)
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!OMP END DO NOWAIT

#ifndef SLICE_MODEL
!     y-change at V-points
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO,jmax+HALO-1
         do i=imin-HALO,imax+HALO
            if (av(i,j) .ge. 1) then ! zero normal gradient across closed boundary
               delfyV(i,j) = f(i,j+1,k) - f(i,j,k)
            end if
         end do
      end do
!$OMP END DO NOWAIT

      if (AH_method .eq. 3) then
!$OMP BARRIER

!$OMP DO SCHEDULE(RUNTIME)
         do j=jmin-HALO,jmax+HALO-1
!           interpolation of delfxU to X-points
            do i=imin-HALO,imax+HALO-1
               if(ax(i,j) .ge. 1) then
                  work2d(i,j) = _HALF_ * ( delfxU(i,j) + delfxU(i,j+1) )
               else
!                 Note (KK): delfxX at corners set to 0
!                            delfxX at closed N/S open bdys not 0, but not needed
                  work2d(i,j) = _ZERO_
               end if
            end do

!           interpolation of delfxX to V-points
            do i=imin-HALO+1,imax+HALO-1
!              Note (KK): we only need delfxV(av=[1|2]), therefore settings for
!                         av=0 (due to zero gradient use of dfdxC) and for
!                         av=3 (zero gradient???) can be skipped
               if (av(i,j).eq.1 .or. av(i,j).eq.2) then
                  delfxV(i,j) = _HALF_ * ( work2d(i-1,j) + work2d(i,j) )
               end if
            end do
         end do
!$OMP END DO

!        interpolation of delfyV to X-points
!$OMP DO SCHEDULE(RUNTIME)
         do j=jmin-HALO,jmax+HALO-1
            do i=imin-HALO,imax+HALO-1
               if(ax(i,j) .ge. 1) then
                  work2d(i,j) = _HALF_ * ( delfyV(i,j) + delfyV(i+1,j) )
               else
!                 Note (KK): delfyX at corners set to 0
!                            delfyX at closed W/E open bdys not 0, but not needed
                  work2d(i,j) = _ZERO_
               end if
            end do
         end do
!$OMP END DO
!        interpolation of delfyX to U-points
!$OMP DO SCHEDULE(RUNTIME)
         do j=jmin-HALO+1,jmax+HALO-1
            do i=imin-HALO,imax+HALO-1
!              Note (KK): we only need delfyU(au=[1|2]), therefore settings for
!                         au=0 (due to zero gradient use of dfdyC) and for
!                         au=3 (zero gradient???) can be skipped
               if (au(i,j).eq.1 .or. au(i,j).eq.2) then
                  delfyU(i,j) = _HALF_ * ( work2d(i,j-1) + work2d(i,j) )
               end if
            end do
         end do
!$OMP END DO
      end if
#endif

!$OMP BARRIER

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin,jmax
#endif
         do i=imin-1,imax
!           flux at U-points
!           Note (KK): flux at au=3 is not needed (AmU would be trash anyway),
!                      at closed boundaries delfxU and diffxy are already 0
!                      therefore only flux calculation at au=[1|2]
            if(au(i,j).eq.1 .or. au(i,j).eq.2) then
               dfdxU = delfxU(i,j) / DXU
               select case (AH_method)
                  case(1)
                     delfxU(i,j) = AH_const * dfdxU
                  case(2)
                     delfxU(i,j) = AmU_3d(i,j,k) / AH_Prt * dfdxU
                  case(3)
                     delfxU(i,j) = AH_stirr_const                    &
                                  * (                                &
                                       diffxx(i,j,k)*dfdxU           &
#ifndef SLICE_MODEL
                                     + diffxy(i,j,k)*delfyU(i,j)/DYU &
#endif
                                    )
               end select
               if (calc_phymix) then
                  phymixU(i,j) = _TWO_ * delfxU(i,j) * dfdxU
               end if
               delfxU(i,j) = DYU*_HALF_*(hn(i,j,k)+hn(i+1,j,k))*delfxU(i,j)
!              now delfxU is flux!
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO NOWAIT

#ifndef SLICE_MODEL
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-1,jmax
         do i=imin,imax
!           flux at V-points
!           Note (KK): flux at av=3 is not needed (AmV would be trash anyway),
!                      at closed boundaries delfyV and diffyx are already 0,
!                      therefore only flux calculation at av=[1|2]
            if(av(i,j).eq.1 .or. av(i,j).eq.2) then
               dfdyV = delfyV(i,j) / DYV
               select case (AH_method)
                  case(1)
                     delfyV(i,j) = AH_const * dfdyV
                  case(2)
                     delfyV(i,j) = AmV_3d(i,j,k) / AH_Prt * dfdyV
                  case(3)
                     delfyV(i,j) = AH_stirr_const                    &
                                  * (                                &
                                       diffyy(i,j,k)*dfdyV           &
                                     + diffyx(i,j,k)*delfxV(i,j)/DXV &
                                    )
               end select
               if (calc_phymix) then
                  phymixV(i,j) = _TWO_ * delfyV(i,j) * dfdyV
               end if
               delfyV(i,j) = DXV*_HALF_*(hn(i,j,k)+hn(i,j+1,k))*delfyV(i,j)
!              now delfyV is flux!
            end if
         end do
      end do
!$OMP END DO NOWAIT
#endif

!$OMP BARRIER
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin,jmax
#endif
         do i=imin,imax
!           Note (KK): tracers at az=2 will be set by external boundary data
!                      therefore dfdx(au=3) and dfdy(av=3) not needed
            if (az(i,j) .eq. 1) then
               if (calc_phymix) then
                  phymix(i,j,k) = _HALF_*(                                  &
                                            phymixU(i-1,j  ) + phymixU(i,j) &
#ifndef SLICE_MODEL
                                          + phymixV(i  ,j-1) + phymixV(i,j) &
#endif
                                         )
               end if
               f(i,j,k) =  f(i,j,k)                              &
                         + dt * (                                &
                                   delfxU(i,j) - delfxU(i-1,j  ) &
#ifndef SLICE_MODEL
                                 + delfyV(i,j) - delfyV(i  ,j-1) &
#endif
                                )                                &
                                * ARCD1 / hn(i,j,k)
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO

#ifdef SLICE_MODEL
!$OMP SINGLE
      if (calc_phymix) then
         phymix(imin:imax,j+1,k) = phymix(imin:imax,j,k)
      end if
      f(imin:imax,j+1,k) = f(imin:imax,j,k)
!$OMP END SINGLE
#endif

   end do

!$OMP END PARALLEL

   call toc(TIM_TRACEDIFF)
#ifdef DEBUG
   write(debug,*) 'Leaving tracer_diffusion()'
   write(debug,*)
#endif
   return
   end subroutine tracer_diffusion
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2011 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
