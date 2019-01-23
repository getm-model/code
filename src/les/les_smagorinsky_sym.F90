#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: les_smagorinsky_sym - \label{les_smagorinsky_sym}
!
! !INTERFACE:
   subroutine les_smagorinsky_sym(dudxC,dudxV,   &
#ifndef SLICE_MODEL
                                  dvdyC,dvdyU,   &
#endif
                                  shearX,shearU, &
                                  AmC,AmX,AmU,AmV)

!  Note (KK): keep in sync with interface in les.F90
!
! !DESCRIPTION:
!
!
! !USES:
   use variables_les, only: SmagC2_2d,SmagX2_2d,SmagU2_2d,SmagV2_2d
   use domain, only: imin,imax,jmin,jmax,az,ax,au,av
   use domain, only: dxc,dyc,dxx,dyx,dxu,dyu,dxv,dyv
   use getm_timers, only: tic,toc,TIM_SMAG2D
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(in) :: dudxC,dudxV
#ifndef SLICE_MODEL
   REALTYPE,dimension(E2DFIELD),intent(in) :: dvdyC,dvdyU
#endif
   REALTYPE,dimension(E2DFIELD),intent(in) :: shearX,shearU
!
! !OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(out),optional :: AmC,AmX,AmU,AmV
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
   REALTYPE :: dudx,dvdy
   integer  :: i,j

!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'les_smagorinsky_sym() # ',Ncall
#endif
#ifdef SLICE_MODEL
   j = jmax/2 ! this MUST NOT be changed!!!
#endif
   call tic(TIM_SMAG2D)

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j)                                         &
!$OMP          PRIVATE(i,dudx,dvdy)

   if (present(AmC)) then
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin-1,jmax+1
#endif
         do i=imin-1,imax+1
!           Note (KK): AmC(az=1) needed in uv_diffusion
            if (az(i,j) .eq. 1) then
!              interpolate shearC
!              Note (KK): in W/E open boundary cells shearC(az=2) would
!                         require shearU outside open boundary
!                         in N/S open boundary cells shearC(az=2) would
!                         require shearU(au=3)
!                         however shearC(az=2) not needed
               AmC(i,j) =  (                                           &
                              dudxC(i,j)                               &
#ifndef SLICE_MODEL
                            - dvdyC(i,j)                               &
#endif
                           )**2                                        &
                         + (_HALF_*(shearU(i-1,j) + shearU(i,j)))**2
               AmC(i,j) = SmagC2_2d(i,j)*DXC*DYC*sqrt(AmC(i,j))
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO NOWAIT

#ifdef SLICE_MODEL
!$OMP BARRIER
!$OMP SINGLE
      AmC(imin-1:imax+1,j+1) =  AmC(imin-1:imax+1,j)
!$OMP END SINGLE
#endif

   end if

   if (present(AmX)) then
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin-1,jmax
#endif
         do i=imin-1,imax
            if (ax(i,j) .eq. 1) then
!              interpolate dudxX and dvdyX
               if (av(i,j).eq.3 .or. av(i+1,j).eq.3) then
!                 Note (KK): western/eastern open bdy
                  dudx = _ZERO_
               else
                  dudx = _HALF_*(dudxV(i,j) + dudxV(i+1,j  ))
               end if
#ifndef SLICE_MODEL
               if (au(i,j).eq.3 .or. au(i,j+1).eq.3) then
!                 Note (KK): northern/southern open bdy
                  dvdy = _ZERO_
               else
                  dvdy = _HALF_*(dvdyU(i,j) + dvdyU(i  ,j+1))
               end if
#endif
               AmX(i,j) =  (                                           &
                              dudx                                     &
#ifndef SLICE_MODEL
                            - dvdy                                     &
#endif
                           )**2                                        &
                         + shearX(i,j)**2
               AmX(i,j) = SmagX2_2d(i,j)*DXX*DYX*sqrt(AmX(i,j))
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO NOWAIT

#ifdef SLICE_MODEL
!$OMP BARRIER
!$OMP SINGLE
      AmX(imin-1:imax,j-1) = AmX(imin-1:imax,j)
      AmX(imin-1:imax,j+1) = AmX(imin-1:imax,j)
!$OMP END SINGLE
#endif

   end if

   if (present(AmU)) then
!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
      do j=jmin-1,jmax+1
#endif
         do i=imin-1,imax
!           Note (KK): shearU(au=3) not available
!                      (however we only need AmU(au=[1|2]) in tracer_diffusion)
            if(au(i,j).eq.1 .or. au(i,j).eq.2) then
!              interpolate dudxU (see deformation_rates)
               if (au(i,j) .eq. 1) then
                  dudx = _HALF_*(dudxC(i,j) + dudxC(i+1,j))
               else
                  dudx = _ZERO_
               end if
               AmU(i,j) =  (                                           &
                              dudx                                     &
#ifndef SLICE_MODEL
                            - dvdyU(i,j)                               &
#endif
                           )**2                                        &
                         + shearU(i,j)**2
               AmU(i,j) = SmagU2_2d(i,j)*DXU*DYU*sqrt(AmU(i,j))
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
!$OMP END DO NOWAIT

#ifdef SLICE_MODEL
!$OMP BARRIER
!$OMP SINGLE
      AmU(imin-1:imax,j+1) = AmU(imin-1:imax,j)
!$OMP END SINGLE
#endif

   end if

#ifndef SLICE_MODEL
   if (present(AmV)) then
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-1,jmax
         do i=imin-1,imax+1
!           Note (KK): we only need AmV(av=[1|2]) in tracer_diffusion
!                      (shearV(av=3) cannot be calculated anyway)
            if(av(i,j).eq.1 .or. av(i,j).eq.2) then
!              interpolate dvdyV and shearV (see deformation_rates)
               if (av(i,j) .eq. 1) then
                  dvdy = _HALF_*(dvdyC(i,j) + dvdyC(i,j+1))
               else
                  dvdy = _ZERO_
               end if
               AmV(i,j) =  (                                           &
                              dudxV(i,j)                               &
                            - dvdy                                     &
                           )**2                                        &
                         + (_HALF_*(shearX(i-1,j) + shearX(i,j)))**2
               AmV(i,j) = SmagV2_2d(i,j)*DXV*DYV*sqrt(AmV(i,j))
            end if
         end do
      end do
!$OMP END DO NOWAIT
   end if
#endif

!$OMP END PARALLEL

   call toc(TIM_SMAG2D)
#ifdef DEBUG
   write(debug,*) 'Leaving les_smagorinsky_sym()'
   write(debug,*)
#endif
   return
   end subroutine les_smagorinsky_sym
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2011 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
