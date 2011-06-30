!$Id: les_smagorinsky.F90,v 1.11 2009-09-30 11:28:45 bjb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: les_smagorinsky - \label{les_smagorinsky}
!
! !INTERFACE:
   subroutine les_smagorinsky(dudxC,dudxU,dvdyC,dvdyV,shearX,shearU,&
                              Am,AmX,AmU,AmV)
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

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(in)           :: dudxC,dudxU
   REALTYPE,dimension(E2DFIELD),intent(in)           :: dvdyC,dvdyV
   REALTYPE,dimension(E2DFIELD),intent(in)           :: shearX,shearU
   REALTYPE,dimension(E2DFIELD),intent(out),optional :: Am,AmX,AmU,AmV
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
   integer                                  :: i,j

!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'les_smagorinsky() # ',Ncall
#endif

   if (present(Am)) then
#ifdef SLICE_MODEL
      j=jmax/2
#else
      do j=jmin-1,jmax+1
#endif
         do i=imin-1,imax+1
            if (az(i,j) .ge. 1) then
               Am(i,j) =  dudxC(i,j)**2                           &
#ifndef SLICE_MODEL
                        + dvdyC(i,j)**2                           &
#endif
                        + _HALF_*(_HALF_*(shearU(i-1,j) + shearU(i,j)))**2
               Am(i,j) = (smag_const**2)*DXC*DYC*sqrt(_TWO_*Am(i,j))
            end if
         end do
#ifndef SLICE_MODEL
      end do
#else
      Am(imin-1:imax+1,jmax/2+1) =  Am(imin-1:imax+1,jmax/2)
#endif
   end if

   if (present(AmX)) then
#ifdef SLICE_MODEL
      j=jmax/2
#else
      do j=jmin-1,jmax
#endif
         do i=imin-1,imax
            if (ax(i,j) .ge. 1) then
               AmX(i,j) =  (_HALF_*(dudxU(i,j) + dudxU(i  ,j+1)))**2 &
#ifndef SLICE_MODEL
                         + (_HALF_*(dvdyV(i,j) + dvdyV(i+1,j  )))**2 &
#endif
                         + _HALF_*shearX(i,j)**2
               AmX(i,j) = (smag_const**2)*DXX*DYX*sqrt(_TWO_*AmX(i,j))
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
            if (au(i,j) .ge. 1) then
               AmU(i,j) =  dudxU(i,j)**2                           &
#ifndef SLICE_MODEL
                         + (_HALF_*(dvdyC(i,j) + dvdyC(i+1,j)))**2 &
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

   if (present(AmV)) then
#ifdef SLICE_MODEL
      j=jmax/2
#else
      do j=jmin-1,jmax
#endif
         do i=imin-1,imax+1
!           average of shearX to shearV is not straight-forward at av=3
!           however AmV(av=3) is not used in tracer_diffusion
            if (av(i,j).eq.1 .or. av(i,j).eq.2) then
!                             average of dudxC to dudxV is straight-forward
                  AmV(i,j) =  (_HALF_*(dudxC(i,j) + dudxC(i,j+1)))**2          &
#ifndef SLICE_MODEL
                            + dvdyV(i,j)**2                                    &
#endif
!                              average of shearX to shearV is straight-forward
                            + _HALF_*(_HALF_*(shearX(i-1,j) + shearX(i,j)))**2
                  AmV(i,j) = (smag_const**2)*DXV*DYV*sqrt(_TWO_*AmV(i,j))
            end if
         end do
#ifndef SLICE_MODEL
      end do
#else
      AmV(imin-1:imax+1,jmax/2-1) = AmV(imin-1:imax+1,jmax/2)
      AmV(imin-1:imax+1,jmax/2+1) = AmV(imin-1:imax+1,jmax/2)
#endif
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving les_smagorinsky()'
   write(debug,*)
#endif
   return
   end subroutine les_smagorinsky

!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2011 - Knut Klingbeil                                  !
!-----------------------------------------------------------------------
