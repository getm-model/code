#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: tracer_stirring - hor.\ SGS-advection
! \label{sec-tracer-stirring}
!
! !INTERFACE:
   subroutine tracer_stirring()
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax,az,au,av
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxu,dyu,dxv,dyv
#else
   use domain, only: dx,dy
#endif
   use variables_3d, only: dudxC_3d,dvdyC_3d,shearX_3d
   use variables_3d, only: dudxV_3d,dvdyU_3d,shearU_3d
   use variables_3d, only: diffxx,diffxy,diffyx,diffyy
   use getm_timers, only: tic,toc,TIM_STIRR
!$ use omp_lib
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
   REALTYPE :: velgrad,dvdyU=_ZERO_,shear,tension,deformation,factor
   integer  :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'tracer_stirring() # ',Ncall
#endif
#ifdef SLICE_MODEL
   j = jmax/2 ! this MUST NOT be changed!!!
#endif
   call tic(TIM_STIRR)

!  Note (KK): in case of SLICE_MODEL only diffxx
!             is needed in tracer_diffusion

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j,dvdyU)                                   &
!$OMP          PRIVATE(i,k,velgrad,shear,tension,deformation,factor)

   do k=1,kmax

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin-1,jmax+1
#endif
         do i=imin-1,imax
!           Note (KK): we only need diffx[x|y](au=[1|2]) in tracer_diffusion
            if(au(i,j).eq.1 .or. au(i,j).eq.2) then
!              interpolate dudxC to U-points (see deformation_rates)
               if (au(i,j) .eq. 1) then
                  velgrad = _HALF_*(dudxC_3d(i,j,k) + dudxC_3d(i+1,j,k))
               else
                  velgrad = _ZERO_
               end if
#ifndef SLICE_MODEL
               dvdyU = dvdyU_3d(i,j,k)
#endif
               tension = velgrad - dvdyU
               deformation = sqrt(tension**2 + shearU_3d(i,j,k)**2)
               factor = _HALF_*DXU*DYU*(_ONE_+(velgrad+dvdyU)/max(SMALL,deformation))
               diffxx(i,j,k) = factor*(tension+deformation)
#ifndef SLICE_MODEL
               diffxy(i,j,k) = factor*shearU_3d(i,j,k)
#endif
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO NOWAIT

#ifdef SLICE_MODEL
!$OMP BARRIER
!$OMP SINGLE
      diffxx(imin-1:imax,j+1,k) = diffxx(imin-1:imax,j,k)
!$OMP END SINGLE
#endif

#ifndef SLICE_MODEL
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-1,jmax
         do i=imin-1,imax+1
!           Note (KK): we only need diffy[x|y](av=[1|2]) in tracer_diffusion
            if(av(i,j).eq.1 .or. av(i,j).eq.2) then
!              interpolate dvdyC and shearX to V-points (see deformation_rates)
               if (av(i,j) .eq. 1) then
                  velgrad = _HALF_*(dvdyC_3d(i,j,k) + dvdyC_3d(i,j+1,k))
               else
                  velgrad = _ZERO_
               end if
               shear = _HALF_*(shearX_3d(i-1,j,k) + shearX_3d(i,j,k))
               tension = dudxV_3d(i,j,k) - velgrad
               deformation = sqrt(tension**2 + shear**2)
               factor = _HALF_*DXV*DYV*(_ONE_+(dudxV_3d(i,j,k)+velgrad)/max(SMALL,deformation))
               diffyx(i,j,k) = factor*shear
               diffyy(i,j,k) = factor*(-tension+deformation)
            end if
         end do
      end do
!$OMP END DO NOWAIT
#endif

   end do

!$OMP END PARALLEL

   call toc(TIM_STIRR)
#ifdef DEBUG
   write(debug,*) 'Leaving tracer_stirring()'
   write(debug,*)
#endif
   return
   end subroutine tracer_stirring
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2011 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
