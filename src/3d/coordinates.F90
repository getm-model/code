#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:  coordinates - defines the vertical coordinate
! \label{sec-coordinates}
!
! !INTERFACE:
   subroutine coordinates(hotstart)
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
#ifdef SLICE_MODEL
   use domain, only: imin,imax,jmin,jmax,kmax
   use variables_3d, only: kvmin,hvo,hvn
#endif
   use getm_timers, only: tic, toc,TIM_COORDS
   use m3d
   use domain, only: vert_cord
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!   integer, intent(in)                 :: cord_type
!   REALTYPE, intent(in)                :: cord_relax
!   REALTYPE, intent(in)                :: maxdepth
   logical, intent(in)                 :: hotstart
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   logical, save   :: first=.true.
   integer         :: ii
!   integer         :: preadapt=0
#ifdef SLICE_MODEL
   integer          :: i,j,k
#endif
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'coordinates() # ',Ncall
#endif
   call tic(TIM_COORDS)

   if (first) then
      select case (vert_cord)
         case (_SIGMA_COORDS_) ! sigma coordinates
            LEVEL2 'using sigma vertical coordinates'
            call sigma_coordinates(.true.)
         case (_Z_COORDS_) ! z-level
         case (_GENERAL_COORDS_) ! general vertical coordinates
            LEVEL2 'using general vertical coordinates'
            call general_coordinates(.true.,cord_relax,maxdepth)
         case (_HYBRID_COORDS_) ! hybrid vertical coordinates
            LEVEL2 'using hybrid vertical coordinates'
            call hybrid_coordinates(.true.)
STDERR 'coordinates(): hybrid_coordinates not coded yet'
stop
         case (_ADAPTIVE_COORDS_) ! adaptive vertical coordinates
            LEVEL2 'using adaptive vertical coordinates'
#ifndef NO_BAROCLINIC
            call adaptive_coordinates(.true.,hotstart)
#endif
         case default
      end select
      first = .false.
   else
      select case (vert_cord)
         case (_SIGMA_COORDS_) ! sigma coordinates
            call sigma_coordinates(.false.)
         case (_Z_COORDS_) ! z-level
         case (_GENERAL_COORDS_) ! general vertical coordinates
            call general_coordinates(.false.,cord_relax,maxdepth)
         case (_HYBRID_COORDS_) ! hybrid vertical coordinates
            call hybrid_coordinates(.false.)
         case (_ADAPTIVE_COORDS_) ! adaptive vertical coordinates
#ifndef NO_BAROCLINIC
            call adaptive_coordinates(.false.,hotstart)
#endif
         case default
      end select
   end if ! first

#ifdef SLICE_MODEL
   do i=imin,imax
      do k=kvmin(i,2),kmax
         hvo(i,1,k)=hvo(i,2,k)
         hvo(i,3,k)=hvo(i,2,k)
         hvn(i,1,k)=hvn(i,2,k)
         hvn(i,3,k)=hvn(i,2,k)
      end do
   end do
#endif

   call toc(TIM_COORDS)
#ifdef DEBUG
   write(debug,*) 'Leaving coordinates()'
   write(debug,*)
#endif
   return
   end subroutine coordinates
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
