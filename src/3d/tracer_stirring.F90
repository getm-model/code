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
   use domain, only: imin,imax,jmin,jmax,kmax,az,au,av,ax
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxc,dyc,arud1,dxx,dyx,arvd1
#else
   use domain, only: dx,dy,ard1
#endif
   use variables_3d, only: dudxC_3d,dvdyC_3d,shearX_3d
   use variables_3d, only: dudxU_3d,dvdyV_3d,shearU_3d
   use variables_3d, only: diffxx,diffxy,diffyx,diffyy

   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
   REALTYPE :: dudxV,dvdyU=_ZERO_,shearV,tension,deformation,factor
   integer  :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'tracer_stirring() # ',Ncall
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!data ranges                                    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!diffxx(imin-1   :imax     ,jmin-1   :jmax+1   )!!!
!!!diffxy(imin-1   :imax     ,jmin-1   :jmax+1   )!!!
!!!diffyx(imin-1   :imax+1   ,jmin-1   :jmax     )!!!
!!!diffyy(imin-1   :imax+1   ,jmin-1   :jmax     )!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!requires:                                      !!!
!!!f     (imin-1   :imax+1   ,jmin-1   :jmax+1   )!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do k=1,kmax

#ifdef SLICE_MODEL
      j=jmax/2
#else
      do j=jmin-1,jmax+1
#endif
         do i=imin-1,imax
#ifndef SLICE_MODEL
            if (au(i,j) .eq. 0) then
               dvdyU=_ZERO_
            else
               dvdyU = _HALF_*(dvdyC_3d(i,j,k) + dvdyC_3d(i+1,j,k))
            end if
#endif
            tension = dudxU_3d(i,j,k)-dvdyU
            deformation = sqrt(tension**2 + shearU_3d(i,j,k)**2)
            factor = _HALF_*DXU*DYU*(_ONE_+(dudxU_3d(i,j,k)+dvdyU)/deformation)
            diffxx(i,j,k) = factor*(tension+deformation)
            diffxy(i,j,k) = factor*shearU_3d(i,j,k)
         end do
#ifndef SLICE_MODEL
      end do
#endif

#ifndef SLICE_MODEL
      do j=jmin-1,jmax
         do i=imin-1,imax+1
            if (av(i,j) .eq. 0) then
               dudxV = _ZERO_
               shearV = _ZERO_
            else
               dudxV = _HALF_*(dudxC_3d(i,j,k) + dudxC_3d(i,j+1,k))
               shearV = _HALF_*(shearX_3d(i-1,j,k) + shearX_3d(i,j,k)) ! not valid for av=3
!              however diffusivities are not used at av=3 in tracer_diffusion
            end if
            tension=dudxV-dvdyV_3d(i,j,k)
            deformation = sqrt(tension**2 + shearV**2)
            factor = _HALF_*DXV*DYV*(_ONE_+(dudxV+dvdyV_3d(i,j,k))/deformation)
            diffyx(i,j,k) = factor*shearV
            diffyy(i,j,k) = factor*(-tension+deformation)
         end do
      end do
#endif
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving tracer_stirring()'
   write(debug,*)
#endif
   return
   end subroutine tracer_stirring
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
