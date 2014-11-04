#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: rs_force - calculates layer-integrated Radiation Stress
!
! !INTERFACE:
   subroutine rs_force(SE,gradterms,UEx,VEx)
!
! !USES:
   use domain         , only: imin,imax,jmin,jmax,au,av,ax
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain         , only: dxc,dyc,dxx,dyx,dxu,dyv,arud1,arvd1
#else
   use domain         , only: dx,dy,ard1
#endif
   use variables_waves, only: coswavedir,sinwavedir
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(in)    :: SE,gradterms
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(inout) :: UEx,VEx
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  Original author(s): Ulf Graewe
!                      Saeed Moghimi
!                      Knut Klingbeil
!
! !LOCAL VARIABLES
   REALTYPE,dimension(E2DFIELD) :: Sxx,Sxy,SxyX
#ifndef SLICE_MODEL
   REALTYPE,dimension(E2DFIELD) :: Syy
#endif
   integer                      :: i,j
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'rs_force() # ',Ncall
#endif

!  parts of the layer-integrated Radiation Stress
   Sxx = coswavedir**2                 * SE
   Sxy = coswavedir    * sinwavedir    * SE
#ifndef SLICE_MODEL
   Syy =                 sinwavedir**2 * SE
#endif

!  average Sxy to V-points
!  Note (KK): avg. to U-points disables threading of j loop in final
!             calculation of SxyX (see OMP-NOTE below)
   do j=jmin-HALO,jmax+HALO-1
      do i=imin-HALO,imax+HALO
         if ( av(i,j) .ne. 0 ) then
            SxyX(i,j) = _HALF_ * ( Sxy(i,j) + Sxy(i,j+1) )
         end if
      end do
   end do

!  average Sxy to X-points
!  OMP-NOTE (KK): i loop must not be changed and cannot be threaded!
   do j=jmin-HALO,jmax+HALO-1
      do i=imin-HALO,imax+HALO-1
         if ( ax(i,j) .eq. 1 ) then
            SxyX(i,j) = _HALF_ * ( SxyX(i,j) + SxyX(i+1,j) )
         else
            SxyX(i,j) = _ZERO_
         end if
      end do
   end do

!  layer-integrated force in x-direction
   do j=jmin-HALO+1,jmax+HALO-1
      do i=imin-HALO,imax+HALO-1
         if ( au(i,j).eq.1 .or. au(i,j).eq.2 ) then
            UEx(i,j) =  UEx(i,j)                                   &
!           first part originates from divergence
                      + (                                          &
                           DYCIP1*Sxx (i+1,j) - DYC   *Sxx (i,j  ) &
#ifndef SLICE_MODEL
                         + DXX   *SxyX(i  ,j) - DXXJM1*SxyX(i,j-1) &
#endif
                        ) * ARUD1                                  &
!           second (and third) part originate from gradient
                      + (                                          &
                           gradterms(i+1,j) - gradterms(i,j)       &
                        ) / DXU
         end if
      end do
   end do

!  layer-integrated force in y-direction
   do j=jmin-HALO,jmax+HALO-1
      do i=imin-HALO+1,imax+HALO-1
         if ( av(i,j).eq.1 .or. av(i,j).eq.2 ) then
            VEx(i,j) =  VEx(i,j)                                   &
!           first part originates from divergence
                      + (                                          &
                           DYX   *SxyX(i,j  ) - DYXIM1*SxyX(i-1,j) &
#ifndef SLICE_MODEL
                         + DXCJP1*Syy (i,j+1) - DXC   *Syy (i  ,j) &
#endif
                        ) * ARVD1                                  &
#ifndef SLICE_MODEL
!           second (and third) part originate from gradient
                      + (                                          &
                           gradterms(i,j+1) - gradterms(i,j)       &
                        ) / DYV
#else
                      & !comment only needed for line continuation :-)
#endif
         end if
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving rs_force()'
   write(debug,*)
#endif
   return
   end subroutine rs_force
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2013 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
