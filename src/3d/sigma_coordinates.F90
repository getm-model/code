!$Id: sigma_coordinates.F90,v 1.2 2007-06-07 10:25:19 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:  equidistant and zoomed sigma-coordinates
! \label{sec-sigam-coordinates}
!
! !INTERFACE:
   subroutine sigma_coordinates(first)
!
! !DESCRIPTION:
!
! Here, the vertical layer distribution in T-, U- and V-points is updated
! during every macro time step. This is done for the old and the new
! layer thicknesses at every point. Calculation of the layer distribution
! in the U- and V-points is done indepently from the calculation in the
! T-points, since different methods for the calculation of the 
! bathymetry values in the U- and V-points are possible, see routine
! {\tt uv\_depths} described on page \pageref{sec-uv-depth}.
!
! Here, three different methods for the vertical layer distribution
! are coded:
!
! \begin{enumerate}
! \item Classical $\sigma$ coordinates where layer interfaces for each 
! layer index have a fixed relative position $\sigma_k$ in the water column,
! which may be even equidistant or non-equidistant, see equations 
! (\ref{sigma}) and (\ref{formula_Antoine}). 
! The surface and bottom zooming factors 
! $d_u$ and $d_l$ are read in via the {\tt domain} namelist in {\tt getm.inp}
! as {\tt ddu} and {\tt ddl}.
! In the first call to coordinates, the relative interface positions
! {\tt dga} are calculated as a one-dimensional vector (in case of
! non-equidistant $\sigma$ coordinates), and those are then multiplied with
! the water depths in all T-, U- and V-points to get the layer thicknesses. 
! \item Also $z$- (i.e.\ geopotential) coordinates are enabled in GETM
! in principle. However, they may not yet work and need further
! development. First of all, fixed $z$-levels are defined by means of
! zooming factors and the maximum water depth $H_{\max}$:
!
! \begin{equation}\label{formula_Antoine_zlevels}
! z_k = H_{\max}\left(\frac{\mbox{tanh}\left( (d_l+d_u)(1+\sigma_k)-d_l\right)
! +\mbox{tanh}(d_l)}{\mbox{tanh}(d_l)+\mbox{tanh}(d_u)}-1\right),
! \qquad k=0,\dots,N\qquad
! \end{equation}
!
! Then, layers are from the surface down filled into the T-point 
! water column locally.
! When the last layer is shallower than {\tt hnmin} (hard coded as local
! variable), the two last layers are combined. The index of the lowest 
! layer is then stored in the integer field {\tt kmin\_pmz}.
! layer thicknesses in U- and V-points are then taken as the minimum 
! values of adjacent thicknesses in T-points, and bottom indices
! {\tt kumin\_pmz} and  {\tt kvmin\_pmz} are taken as the maximum
! of adjacent  {\tt kmin\_pmz} indices.
! \item The third and so far most powerful method are the genral
! vertical coordinates, discussed in section \ref{SectionGeneralCoordinates},
! see equations (\ref{sigma}) - (\ref{MLDTransform}), which is basically
! an interpolation between equidistant and non-equaidistant $\sigma$
! coordinates. During the first call, a three-dimensional field
! {\tt gga} containing the relative interface positions is calculated,
! which further down used together with the actual water depth in the 
! T-, U- and V-points for calculating the updated old and new layer
! thicknesses.
!\end{enumerate}
! 
! A fourth option will soon be the adaptive grids which have been
! conceptionally developed by \cite{BURCHARDea04}.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax,H,HU,HV
   use domain, only: ga,ddu,ddl
   use variables_3d, only: kmin,kumin,kvmin,ho,hn,huo,hun,hvo,hvn
   use variables_3d, only: sseo,ssen,ssuo,ssun,ssvo,ssvn
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical, intent(in)                  :: first
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer          :: i,j,k,rc
   REALTYPE         :: kmaxm1
   logical, save    :: equiv_sigma=.false.
   REALTYPE, save, dimension(:), allocatable  :: dga
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'coordinates() # ',Ncall
#endif

   if (first) then
      if (.not. allocated(ga)) allocate(ga(0:kmax),stat=rc)
      if (rc /= 0) stop 'coordinates: Error allocating (ga)'
      if (ddu .le. _ZERO_ .and. ddl .le. _ZERO_) then
         equiv_sigma=.true.
         ga(0) = -_ONE_
         do k=1,kmax
            ga(k) = ga(k-1) + _ONE_/kmax
         end do
         ga(kmax) = _ZERO_
      else
         ! Non-equidistant sigma coordinates
         ! This zooming routine is from Antoine Garapon, ICCH, DK
         if (ddu .lt. _ZERO_) ddu=_ZERO_
         if (ddl .lt. _ZERO_) ddl=_ZERO_
         allocate(dga(0:kmax),stat=rc)
         if (rc /= 0) STOP 'coordinates: Error allocating (dga)'
         ga(0)= -_ONE_
         dga(0)= _ZERO_
         do k=1,kmax
            ga(k)=tanh((ddl+ddu)*k/float(kmax)-ddl)+tanh(ddl)
            ga(k)=ga(k)/(tanh(ddl)+tanh(ddu)) - _ONE_
            dga(k)=ga(k)-ga(k-1)
         end do
      end if
      kmin=1
      kumin=1
      kvmin=1
   end if ! first

   if (equiv_sigma) then
      kmaxm1= _ONE_/float(kmax)
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO
            ho(i,j,:)=(sseo(i,j)+H(i,j))*kmaxm1
            hn(i,j,:)=(ssen(i,j)+H(i,j))*kmaxm1
         end do
      end do

      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO-1
            huo(i,j,:)=(ssuo(i,j)+HU(i,j))*kmaxm1
            hun(i,j,:)=(ssun(i,j)+HU(i,j))*kmaxm1
         end do
      end do

      do j=jmin-HALO,jmax+HALO-1
         do i=imin-HALO,imax+HALO
            hvo(i,j,:)=(ssvo(i,j)+HV(i,j))*kmaxm1
            hvn(i,j,:)=(ssvn(i,j)+HV(i,j))*kmaxm1
         end do
      end do

   else ! non-equivdistant

      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO
            ho(i,j,1:kmax)=(sseo(i,j)+H(i,j))*dga(1:kmax)
            hn(i,j,1:kmax)=(ssen(i,j)+H(i,j))*dga(1:kmax)
         end do
      end do

      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO-1
            huo(i,j,1:kmax)=(ssuo(i,j)+HU(i,j))*dga(1:kmax)
            hun(i,j,1:kmax)=(ssun(i,j)+HU(i,j))*dga(1:kmax)
         end do
      end do

      do j=jmin-HALO,jmax+HALO-1
         do i=imin-HALO,imax+HALO
            hvo(i,j,1:kmax)=(ssvo(i,j)+HV(i,j))*dga(1:kmax)
            hvn(i,j,1:kmax)=(ssvn(i,j)+HV(i,j))*dga(1:kmax)
         end do
      end do
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving sigma_coordinates()'
   write(debug,*)
#endif
   return
   end subroutine sigma_coordinates
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
