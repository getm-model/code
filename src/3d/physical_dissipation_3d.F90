#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: physical_dissipation_3d()
!
! !INTERFACE:
   subroutine physical_dissipation_3d()
!
! !DESCRIPTION:
!
! Here, the physical dissipation of mean kinetic energy is calculated as 
! \begin{equation}
! D^{phy} = D^{phy}_h + D^{phy}_v = 2 \alpha A^M_h S_h^2 +(\nu_t+\nu) S_v^2
! \end{equation}
! with the horizontal eddy viscosity, $A^M_h$, the vertical eddy viscosity,
! $\nu_t$, the molecular viscosity, $\nu$, the drying parameter, $\alpha$
! (see equations (\ref{uEq}) and (\ref{vEq})), the horizontal shear tensor,
! \begin{equation}
! S_{ij} = \frac12 \left(\partial_{x_i} u_j + \partial_{x_j}u_i \right),
! \quad i,j=1,2
! \end{equation}
! (with $x_1=x$, $x_2=y$, $u_1=u$ and $u_2=u$) and the scalar vertical shear
! \begin{equation}
! S_v = \left[\left(\partial_zu\right)^2
! +\left(\partial_zv\right)^2\right]^{1/2}.
! \end{equation}
! The horizontal shear square is then denoted as
! \begin{equation}
! S_h^2 = \sum_{i,j=1}^2 
! \left[\frac12 \left(\partial_{x_i} u_j + \partial_{x_j}u_i \right)\right]^2.
! \end{equation}
! Note that the dissipation term is obtained from the multiplication of the
! $u$-equation (\ref{uEq}) by $u$, the multiplicationof the $v$-equation
! (\ref{vEq}) by $v$ and subsequent addition of the two resulting equations. 
!
! This routine is called only for {\tt save\_numerical\_analyses = .true.}
!
!
! !USES:
   use domain,       only: imin,imax,jmin,jmax,kmax,H,au,av,ax,az
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxc,dyc,dxx,dyx
#else
   use domain, only: dx,dy
#endif
   use domain, only: dry_z
   use m2d, only: Am
   use variables_3d, only: phydis_3d,phydis_int,num,uu,vv,hn,hun,hvn,SS
   use parameters, only: avmmol

   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard
!
! !LOCAL VARIABLES:
   REALTYPE                  :: dupper,dlower
   integer                   :: i,j,k
   REALTYPE                  :: aux(I3DFIELD)
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'physical_dissipation() # ',Ncall
#endif
#ifdef SLICE_MODEL
   j = jmax/2 ! this MUST NOT be changed!!!
#endif
   phydis_int = _ZERO_

   if (Am .gt. _ZERO_) then
      do k=1,kmax
         call physical_dissipation(uu(:,:,k),vv(:,:,k),hun(:,:,k),hvn(:,:,k),Am,phydis_3d(:,:,k))
      end do
   else
      phydis_3d = _ZERO_
   end if

   ! Av * ( (du/dz)**2 + (dv/dz)**2 ) on W-POINTS
   aux(:,:,kmax)=_ZERO_
   aux(:,:,0)   =_ZERO_
   do k=1,kmax-1
#ifndef SLICE_MODEL
      do j=jmin,jmax
#endif
         do i=imin,imax
            if (az(i,j).eq.1) then
               aux(i,j,k)=(num(i,j,k)+avmmol)*SS(i,j,k)
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
   end do 
   ! Add Av * ( (du/dz)**2 + (dv/dz)**2 ) on T-POINTS
   do k=1,kmax
#ifndef SLICE_MODEL
      do j=jmin,jmax
#endif
         do i=imin,imax
            if (az(i,j).eq.1) then
               phydis_3d(i,j,k)=phydis_3d(i,j,k)+_HALF_*(aux(i,j,k-1)+aux(i,j,k))
               phydis_int(i,j)=phydis_int(i,j)+phydis_3d(i,j,k)*hn(i,j,k)
            end if
         end do
#ifndef SLICE_MODEL
      end do
#endif
   end do 

#ifdef DEBUG
   write(debug,*) 'Leaving physical_dissipation_3d()'
   write(debug,*)
#endif
   return
   end subroutine physical_dissipation_3d
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2012 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
