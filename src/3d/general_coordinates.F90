#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:  general vertical coordinates
! \label{sec-general-coordinates}
!
! !INTERFACE:
   subroutine general_coordinates(first,cord_relax,maxdepth)
!
! !DESCRIPTION:
!
! Here, the general vertical coordinates layer distribution
! in T-, U- and V-points is calculated. The general vertical coordinates as
! discussed in section \ref{SectionGeneralCoordinates},
! see equations (\ref{sigma}) - (\ref{MLDTransform}), are basically
! an interpolation between equidistant and non-equaidistant $\sigma$
! coordinates. During the first call, a three-dimensional field
! {\tt gga} containing the relative interface positions is calculated,
! which further down used together with the actual water depth in the
! T-, U- and V-points for calculating the updated old and new layer
! thicknesses.
!
!
! !USES:
   use domain, only: ga,ddu,ddl,d_gamma,gamma_surf
   use domain, only: imin,imax,jmin,jmax,kmax,H,HU,HV,az,au,av,min_depth
   use variables_3d, only: dt,kmin,kumin,kvmin,ho,hn,huo,hun,hvo,hvn
   use variables_3d, only: sseo,ssen,ssuo,ssun,ssvo,ssvn
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical, intent(in)                 :: first
   REALTYPE, intent(in)                :: cord_relax
   REALTYPE, intent(in)                :: maxdepth
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer         :: i,j,k,rc,kk
   REALTYPE        :: alpha
   REALTYPE        :: HH,zz,r
   REALTYPE, save, dimension(:),     allocatable  :: dga,be,sig
   REALTYPE, save, dimension(:,:,:), allocatable  :: gga
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
      do k=0,kmax
         ga(k) = k
      end do
      allocate(sig(0:kmax),stat=rc)    ! dimensionless sigma-coordinate
      if (rc /= 0) STOP 'coordinates: Error allocating (sig)'
      allocate(be(0:kmax),stat=rc)     ! dimensionless beta-coordinate
      if (rc /= 0) STOP 'coordinates: Error allocating (be)'
      allocate(gga(I3DFIELD),stat=rc)  ! dimensionless gamma-coordinate
      if (rc /= 0) stop 'coordinates: Error allocating memory (gga)'
      be(0)=  -_ONE_
      sig(0)= -_ONE_
      if (ddu .lt. _ZERO_) ddu=_ZERO_
      if (ddl .lt. _ZERO_) ddl=_ZERO_
      do k=1,kmax
         be(k)=tanh((ddl+ddu)*k/float(kmax)-ddl)+tanh(ddl)
         be(k)=be(k)/(tanh(ddl)+tanh(ddu))-_ONE_
         sig(k)=k/float(kmax)-_ONE_
      end do
      if (gamma_surf) then
         kk=kmax
      else
         kk=1
      end if
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO
            HH=max(sseo(i,j)+H(i,j),min_depth)
            alpha=min(&
                     ((be(kk)-be(kk-1))-D_gamma/HH&
                      *(sig(kk)-sig(kk-1)))&
                      /((be(kk)-be(kk-1))-(sig(kk)-sig(kk-1))),_ONE_)
            gga(i,j,0)=-_ONE_
            do k=1,kmax
               gga(i,j,k)=alpha*sig(k)+(1.-alpha)*be(k)
               if (gga(i,j,k) .lt. gga(i,j,k-1)) then
                  STDERR kk,(be(kk)-be(kk-1)),(sig(kk)-sig(kk-1))
                  STDERR D_gamma,HH
                  STDERR alpha
                  STDERR k-1,gga(i,j,k-1),be(k-1),sig(k-1)
                  STDERR k,gga(i,j,k),be(k),sig(k)
                  stop 'coordinates'
               end if
            end do
         end do
      end do

      kmin=1
      kumin=1
      kvmin=1

!     Here, the initial layer distribution is calculated.
      do k=1,kmax
         do j=jmin-HALO,jmax+HALO
            do i=imin-HALO,imax+HALO
               HH=max(sseo(i,j)+H(i,j),min_depth)
               hn(i,j,k)=HH*(gga(i,j,k)-gga(i,j,k-1))
            end do
         end do
      end do

! BJB-TODO: Change 0.5 -> _HALF_ (and 1. -> _ONE_) in this file.
      do k=1,kmax
         do j=jmin-HALO,jmax+HALO
            do i=imin-HALO,imax+HALO-1
               HH=max(ssuo(i,j)+HU(i,j),min_depth)
               huo(i,j,k)=HH*0.5*            &
                (gga(i,j,k)-gga(i,j,k-1)+gga(i+1,j,k)-gga(i+1,j,k-1))
               hun(i,j,k)=huo(i,j,k)
            end do
         end do
      end do

      do k=1,kmax
         do j=jmin-HALO,jmax+HALO-1
            do i=imin-HALO,imax+HALO
               HH=max(ssvo(i,j)+HV(i,j),min_depth)
               hvo(i,j,k)=HH*0.5*            &
                (gga(i,j,k)-gga(i,j,k-1)+gga(i,j+1,k)-gga(i,j+1,k-1))
               hvn(i,j,k)=hvo(i,j,k)
            end do
         end do
      end do
   end if ! first

! The general vertical coordinates can be relaxed towards the new layer
! thicknesses by the following relaxation time scale r. This should
! later be generalised also for sigma coordinates.

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,rc,kk,HH,zz,r)
!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO
         r=cord_relax/dt*H(i,j)/maxdepth
         HH=ssen(i,j)+H(i,j)
         if (HH .lt. D_gamma) then
            do k=1,kmax
               ho(i,j,k)=hn(i,j,k)
               hn(i,j,k)=HH/kmax
            end do
         else
            zz=-H(i,j)
            do k=1,kmax-1
               ho(i,j,k)=hn(i,j,k)
               hn(i,j,k)=(ho(i,j,k)*r+HH*(gga(i,j,k)-gga(i,j,k-1)))/(r+1.)
               zz=zz+hn(i,j,k)
            end do
            ho(i,j,kmax)=hn(i,j,kmax)
            hn(i,j,kmax)=ssen(i,j)-zz
         end if
      end do
   end do
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO-1
!KBK         if (au(i,j) .gt. 0) then
            r=cord_relax/dt*HU(i,j)/maxdepth
            zz=-HU(i,j)
            HH=ssun(i,j)+HU(i,j)
            do k=1,kmax-1
               huo(i,j,k)=hun(i,j,k)
               hun(i,j,k)=(huo(i,j,k)*r+HH*0.5*(gga(i,j,k)-gga(i,j,k-1) &
                         +gga(i+1,j,k)-gga(i+1,j,k-1)))/(r+1.)
               zz=zz+hun(i,j,k)
            end do
            huo(i,j,kmax)=hun(i,j,kmax)
            hun(i,j,kmax)=ssun(i,j)-zz
!KBK         end if
      end do
   end do
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(RUNTIME)
   do j=jmin-HALO,jmax+HALO-1
      do i=imin-HALO,imax+HALO
!KBK         if (av(i,j).gt.0) then
            r=cord_relax/dt*HV(i,j)/maxdepth
            zz=-HV(i,j)
            HH=ssvn(i,j)+HV(i,j)
            do k=1,kmax-1
               hvo(i,j,k)=hvn(i,j,k)
               hvn(i,j,k)=(hvo(i,j,k)*r+HH*0.5*(gga(i,j,k)-gga(i,j,k-1) &
                         +gga(i,j+1,k)-gga(i,j+1,k-1)))/(r+1.)
               zz=zz+hvn(i,j,k)
            end do
            hvo(i,j,kmax)=hvn(i,j,kmax)
            hvn(i,j,kmax)=ssvn(i,j)-zz
!KBK         end if
      end do
   end do
!$OMP END DO
!$OMP END PARALLEL

#ifdef DEBUG
   write(debug,*) 'Leaving general_coordinates()'
   write(debug,*)
#endif
   return
   end subroutine general_coordinates
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
