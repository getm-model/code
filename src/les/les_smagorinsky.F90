#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: les_smagorinsky - \label{les_smagorinsky}
!
! !INTERFACE:
   subroutine les_smagorinsky(dudxC,dudxV,   &
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
   use les, only: smag_const
   use domain, only: imin,imax,jmin,jmax,az,ax,au,av
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxc,dyc,dxx,dyx,dxu,dyu,dxv,dyv
#else
   use domain, only: dx,dy
#endif
   use getm_timers, only: tic,toc,TIM_SMAG2D

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
   REALTYPE :: dudx
#ifndef SLICE_MODEL
   REALTYPE :: dvdy
#endif
   integer  :: i,j

!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'les_smagorinsky() # ',Ncall
#endif
   call tic(TIM_SMAG2D)

   if (present(AmC)) then
#ifdef SLICE_MODEL
      j=jmax/2
#else
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
               AmC(i,j) =  dudxC(i,j)**2                                    &
#ifndef SLICE_MODEL
                         + dvdyC(i,j)**2                                    &
#endif
                         + _HALF_*(_HALF_*(shearU(i-1,j) + shearU(i,j)))**2
               AmC(i,j) = (smag_const**2)*DXC*DYC*sqrt(_TWO_*AmC(i,j))
            end if
         end do
#ifndef SLICE_MODEL
      end do
#else
      AmC(imin-1:imax+1,jmax/2+1) =  AmC(imin-1:imax+1,jmax/2)
#endif
   end if

   if (present(AmX)) then
#ifdef SLICE_MODEL
      j=jmax/2
#else
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
               AmX(i,j) =  dudx**2 &
#ifndef SLICE_MODEL
                         + dvdy**2 &
#endif
                         + _HALF_*shearX(i,j)**2
               AmX(i,j) = (smag_const**2)*DXX*DYX*sqrt(_TWO_*AmX(i,j))
#ifdef _CORRECT_METRICS_
#if defined(SPHERICAL) || defined(CURVILINEAR)
            else
!              Note (KK): shearX at corners already set to 0, there no AmX needed
               if (av(i,j).eq.0 .and. av(i+1,j).eq.0) then
                  if (au(i,j) .eq. 1) then ! northern closed boundary
                     AmX(i,j) =  (_HALF_*(dudxV(i,j) + dudxV(i+1,j)))**2 &
#ifndef SLICE_MODEL
                               + dvdyU(i,j  )**2                         &
#endif
                               + _HALF_*shearX(i,j)**2
                     AmX(i,j) = (smag_const**2)*DXX*DYX*sqrt(_TWO_*AmX(i,j))
                  end if
                  if (au(i,j+1) .eq. 1) then ! southern closed boundary
                     AmX(i,j) =  (_HALF_*(dudxV(i,j) + dudxV(i+1,j)))**2 &
#ifndef SLICE_MODEL
                               + dvdyU(i,j+1)**2                         &
#endif
                               + _HALF_*shearX(i,j)**2
                     AmX(i,j) = (smag_const**2)*DXX*DYX*sqrt(_TWO_*AmX(i,j))
                  end if
               end if
               if (au(i,j).eq.0 .and. au(i,j+1).eq.0) then
                  if (av(i,j) .eq. 1) then ! eastern closed boundary
                     AmX(i,j) =  dudxV(i  ,j)**2                         &
#ifndef SLICE_MODEL
                               + (_HALF_*(dvdyU(i,j) + dvdyU(i,j+1)))**2 &
#endif
                               + _HALF_*shearX(i,j)**2
                     AmX(i,j) = (smag_const**2)*DXX*DYX*sqrt(_TWO_*AmX(i,j))
                  end if
                  if (av(i+1,j) .eq. 1) then ! western closed boundary
                     AmX(i,j) =  dudxV(i+1,j)**2                         &
#ifndef SLICE_MODEL
                               + (_HALF_*(dvdyU(i,j) + dvdyU(i,j+1)))**2 &
#endif
                               + _HALF_*shearX(i,j)**2
                     AmX(i,j) = (smag_const**2)*DXX*DYX*sqrt(_TWO_*AmX(i,j))
                  end if
               end if
#endif
#endif
            end if
         end do
#ifndef SLICE_MODEL
      end do
#else
      AmX(imin-1:imax,jmax/2-1) = AmX(imin-1:imax,jmax/2)
      AmX(imin-1:imax,jmax/2+1) = AmX(imin-1:imax,jmax/2)
#endif
   end if

   if (present(AmU)) then
#ifdef SLICE_MODEL
      j=jmax/2
#else
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
               AmU(i,j) =  dudx**2                                     &
#ifndef SLICE_MODEL
                         + dvdyU(i,j)**2                               &
#endif
                         + _HALF_*shearU(i,j)**2
               AmU(i,j) = (smag_const**2)*DXU*DYU*sqrt(_TWO_*AmU(i,j))
            end if
         end do
#ifndef SLICE_MODEL
      end do
#else
      AmU(imin-1:imax,jmax/2+1) = AmU(imin-1:imax,jmax/2)
#endif
   end if

#ifndef SLICE_MODEL
   if (present(AmV)) then
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
               AmV(i,j) =  dudxV(i,j)**2                                    &
                         + dvdy**2                                          &
                         + _HALF_*(_HALF_*(shearX(i-1,j) + shearX(i,j)))**2
               AmV(i,j) = (smag_const**2)*DXV*DYV*sqrt(_TWO_*AmV(i,j))
            end if
         end do
      end do
   end if
#endif

   call toc(TIM_SMAG2D)
#ifdef DEBUG
   write(debug,*) 'Leaving les_smagorinsky()'
   write(debug,*)
#endif
   return
   end subroutine les_smagorinsky

!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
