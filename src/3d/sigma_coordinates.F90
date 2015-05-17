#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  equidistant and zoomed sigma-coordinates
! \label{sec-sigma-coordinates}
!
! !INTERFACE:
   subroutine sigma_coordinates(first,hotstart)
!
! !DESCRIPTION:
!
! Here, the sigma coordinates layer distribution in T-, U- and V-points is calculated.
! The layer interfaces for each
! layer index have a fixed relative position $\sigma_k$ in the water column,
! which may be even equidistant or non-equidistant, see equations
! (\ref{sigma}) and (\ref{formula_Antoine}).
! The surface and bottom zooming factors
! $d_u$ and $d_l$ are read in via the {\tt domain} namelist in {\tt getm.inp}
! as {\tt ddu} and {\tt ddl}.
! In the first call to the {\tt sigma\_coordinates}, the relative interface positions
! {\tt dga} are calculated as a one-dimensional vector (in case of
! non-equidistant $\sigma$ coordinates), and those are then multiplied with
! the water depths in all T-, U- and V-points to get the layer thicknesses.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax,H
   use domain, only: ga,ddu,ddl
   use variables_3d, only: kmin,kumin,kvmin,ho,hn,hun,hvn
   use variables_3d, only: Dn,Dun,Dvn,sseo
   use vertical_coordinates,only: restart_with_ho,restart_with_hn
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical, intent(in)                  :: first
   logical, intent(in)                  :: hotstart
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer          :: i,j,k,rc
   REALTYPE, save, dimension(:), allocatable  :: dga
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'coordinates() # ',Ncall
#endif

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
!$OMP SINGLE

   if (first) then

      if (.not. allocated(ga)) then
         allocate(ga(0:kmax),stat=rc)
         if (rc /= 0) stop 'coordinates: Error allocating (ga)'
      end if
      allocate(dga(0:kmax),stat=rc)
      if (rc /= 0) STOP 'coordinates: Error allocating (dga)'
      ga(0) = -_ONE_
      ga(kmax) = _ZERO_
      dga(0)= _ZERO_
      if (ddu .le. _ZERO_ .and. ddl .le. _ZERO_) then
         dga(1:kmax) = _ONE_/kmax
         do k=1,kmax-1
            ga(k) = ga(k-1) + dga(k)
         end do
      else
         ! Non-equidistant sigma coordinates
         ! This zooming routine is from Antoine Garapon, ICCH, DK
         if (ddu .lt. _ZERO_) ddu=_ZERO_
         if (ddl .lt. _ZERO_) ddl=_ZERO_
         do k=1,kmax
            ga(k)=tanh((ddl+ddu)*k/float(kmax)-ddl)+tanh(ddl)
            ga(k)=ga(k)/(tanh(ddl)+tanh(ddu)) - _ONE_
            dga(k)=ga(k)-ga(k-1)
         end do
      end if

      if (.not. restart_with_hn) then
         if (hotstart) then
            LEVEL2 'WARNING: assume sigma coordinates for hn'
         end if
         do k=1,kmax
            do j=jmin-HALO,jmax+HALO
               do i=imin-HALO,imax+HALO
                  hn(i,j,k) = Dn(i,j) * dga(k)
               end do
            end do
         end do
      end if
!     only for backward compatibility
      if (.not. restart_with_ho) then
         if (hotstart) then
            LEVEL2 'WARNING: assume sigma coordinates for ho'
         end if
         do k=1,kmax
            do j=jmin-HALO,jmax+HALO
               do i=imin-HALO,imax+HALO
                  ho(i,j,k) = (sseo(i,j)+H(i,j)) * dga(k)
               end do
            end do
         end do
      end if

      kmin=1
      kumin=1
      kvmin=1

   else

!$OMP END SINGLE
      do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
         do j=jmin-HALO,jmax+HALO
            do i=imin-HALO,imax+HALO
               hn(i,j,k)= Dn(i,j) * dga(k)
            end do
         end do
!$OMP END DO NOWAIT
      end do
!$OMP SINGLE

   end if

!$OMP END SINGLE
      do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
         do j=jmin-HALO,jmax+HALO
            do i=imin-HALO,imax+HALO-1
               hun(i,j,k) = Dun(i,j) * dga(k)
            end do
         end do
!$OMP END DO NOWAIT
      end do

      do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
         do j=jmin-HALO,jmax+HALO-1
            do i=imin-HALO,imax+HALO
               hvn(i,j,k) = Dvn(i,j) * dga(k)
            end do
         end do
!$OMP END DO NOWAIT
      end do
!$OMP SINGLE

!$OMP END PARALLEL

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
