!$Id: les_smagorinsky.F90,v 1.11 2009-09-30 11:28:45 bjb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: les_smagorinsky - \label{les_smagorinsky}
!
! !INTERFACE:
   subroutine les_smagorinsky(U,V,DU,DV,Am,AmX,AmU,AmV)
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

   REALTYPE,intent(in)           :: U(E2DFIELD),V(E2DFIELD)
   REALTYPE,intent(in)           :: DU(E2DFIELD),DV(E2DFIELD)
   REALTYPE,intent(out)          :: Am(E2DFIELD),AmX(E2DFIELD)
   REALTYPE,intent(out),optional :: AmU(E2DFIELD),AmV(E2DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
!  must be saved because of dynamic allocation outside a module
   REALTYPE,dimension(:,:),allocatable,save :: dxdu,dxduU,shearX,shearU
#ifndef SLICE_MODEL
   REALTYPE,dimension(:,:),allocatable,save :: dydv,dydvU
#endif
   logical,save                             :: first=.true.
   integer                                  :: rc
   integer                                  :: i,j


!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'les_smagorinsky() # ',Ncall
#endif

   if (first) then
      allocate(dxdu(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'les_smagorinsky: Error allocating memory (dxdu)'
      dxdu   = _ZERO_

      allocate(dxduU(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'les_smagorinsky: Error allocating memory (dxduU)'
      dxduU  = _ZERO_

#ifndef SLICE_MODEL
      allocate(dydv(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'les_smagorinsky: Error allocating memory (dydv)'
      dydv   = _ZERO_

      allocate(dydvU(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'les_smagorinsky: Error allocating memory (dydvU)'
      dydvU  = _ZERO_
#endif

      allocate(shearX(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'les_smagorinsky: Error allocating memory (shearX)'
      shearX = _ZERO_

      allocate(shearU(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'les_smagorinsky: Error allocating memory (shearU)'
      shearU = _ZERO_

      first  = .false.
   end if

!  strain terms at center points
#ifdef SLICE_MODEL
   j=jmax/2
#else
   do j=jmin-1,jmax+1
#endif
      do i=imin-1,imax+1
!        this condition also excludes tangential strain inside open boundary cells! ok?
         if (az(i,j) .eq. 1) then
            dxdu(i,j) = (U(i,j)/DU(i,j) - U(i-1,j)/DU(i-1,j))/DXC
#ifndef SLICE_MODEL
            dydv(i,j) = (V(i,j)/DV(i,j) - V(i,j-1)/DV(i,j-1))/DYC
#endif
         end if
      end do
#ifndef SLICE_MODEL
   end do
#else
   dxdu(:,jmax/2+1) = dxdu(:,jmax/2)
#endif

!  interpolation to U-points
#ifdef SLICE_MODEL
   j=jmax/2
#else
   do j=jmin-1,jmax+1
#endif
      do i=imin-1,imax
         if (au(i,j) .ge. 1) then
            dxduU(i,j) = _HALF_*(dxdu(i,j) + dxdu(i+1,j))
#ifndef SLICE_MODEL
            dydvU(i,j) = _HALF_*(dydv(i,j) + dydv(i+1,j))
#endif
         end if
      end do
#ifndef SLICE_MODEL
   end do
#else
   dxduU(:,jmax/2+1) = dxduU(:,jmax/2)
#endif

!  shear term at X-points
#ifdef SLICE_MODEL
   j=jmax/2
#else
   do j=jmin-1,jmax+1
#endif
      do i=imin-1,imax+1
         if (ax(i,j) .ge. 1) then
            shearX(i,j) =                                &
#ifndef SLICE_MODEL
               (U(i,j+1)/DU(i,j+1) - U(i,j)/DU(i,j))/DYX &
             +                                           &
#endif
               (V(i+1,j)/DV(i+1,j) - V(i,j)/DV(i,j))/DXX
         end if
      end do
#ifndef SLICE_MODEL
   end do
#else
   shearX(:,jmax/2-1) = shearX(:,jmax/2)
#endif

!  interpolation to U-points
#ifdef SLICE_MODEL
   j=jmax/2
#else
   do j=jmin,jmax+1
#endif
      do i=imin-1,imax+1
!        this condition excludes shear inside open boundary cells! ok?
         if (au(i,j) .eq. 1 .or. au(i,j) .eq. 2) then
            shearU(i,j) = _HALF_*(shearX(i,j-1) + shearX(i,j))
         end if
      end do
#ifndef SLICE_MODEL
   end do
#endif

!  Finalise
#ifdef SLICE_MODEL
   j=jmax/2
#else
   do j=jmin,jmax+1
#endif
      do i=imin,imax+1
!        this condition excludes shear inside open boundary cells! ok?
         if (az(i,j) .eq. 1) then
            Am(i,j) =  dxdu(i,j)**2                           &
#ifndef SLICE_MODEL
                     + dydv(i,j)**2                           &
#endif
                     + _HALF_*(shearU(i-1,j) + shearU(i,j))**2
            Am(i,j) = (smag_const**2)*DXC*DYC*sqrt(_TWO_*Am(i,j))
         end if
      end do
#ifndef SLICE_MODEL
   end do
#else
   Am(:,jmax/2+1) =  Am(:,jmax/2)
#endif

#ifdef SLICE_MODEL
   j=jmax/2
#else
   do j=jmin-1,jmax
#endif
      do i=imin-1,imax
         if (ax(i,j) .ge. 1) then
            AmX(i,j) =  (_HALF_*(dxduU(i,j) + dxduU(i,j+1)))**2 &
#ifndef SLICE_MODEL
                      + (_HALF_*(dydvU(i,j) + dydvU(i,j+1)))**2 &
#endif
                      + _HALF_*shearX(i,j)**2
            AmX(i,j) = (smag_const**2)*DXX*DYX*sqrt(_TWO_*AmX(i,j))
         end if
      end do
#ifndef SLICE_MODEL
   end do
#else
   AmX(:,jmax/2-1) = AmX(:,jmax/2)
   AmX(:,jmax/2+1) = AmX(:,jmax/2)
#endif

   if (present(AmU)) then
#ifdef SLICE_MODEL
      j=jmax/2
#else
      do j=jmin,jmax
#endif
         do i=imin-1,imax
            if (au(i,j) .ge. 1) then
               AmU(i,j) =  dxduU(i,j)**2        &
#ifndef SLICE_MODEL
                         + dydvU(i,j)**2        &
#endif
                         + _HALF_*shearU(i,j)**2
               AmU(i,j) = (smag_const**2)*DXU*DYU*sqrt(_TWO_*AmU(i,j))
            end if
         end do
#ifndef SLICE_MODEL
      end do
#else
      AmU(:,jmax/2+1) = AmU(:,jmax/2)
#endif
   end if

   if (present(AmV)) then
#ifdef SLICE_MODEL
      j=jmax/2
#else
      do j=jmin-1,jmax
#endif
         do i=imin,imax
!           this condition excludes shear inside open boundary cells! ok?
            if (av(i,j) .eq. 1 .or. av(i,j) .eq. 2) then
               AmV(i,j) =  (_HALF_*(dxdu(i,j) + dxdu(i,j+1)))**2  &
#ifndef SLICE_MODEL
                         + (_HALF_*(dydv(i,j) + dydv(i,j+1)))**2  &
#endif
                         + _HALF_*(shearX(i-1,j) + shearX(i,j))**2
               AmV(i,j) = (smag_const**2)*DXV*DYV*sqrt(_TWO_*AmV(i,j))
            end if
         end do
#ifndef SLICE_MODEL
      end do
#else
      AmV(:,jmax/2-1) = AmV(:,jmax/2)
      AmV(:,jmax/2+1) = AmV(:,jmax/2)
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
