#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: physical_dissipation_3d()
!
! !INTERFACE:
   subroutine physical_dissipation_3d(Am,pd3d,pd2d)
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
   use domain, only: dxx,dyx,dxc,dyc
#else
   use domain, only: dx,dy,dry_z
#endif
   use domain, only: dry_z
   use variables_3d, only: num,uu,vv,hn,hun,hvn,SS
   use parameters, only: avmmol

   IMPLICIT NONE

! !INPUT PARAMETERS
   REALTYPE, intent(in) :: Am

! !INPUT/OUTPUT PARAMETERS
   REALTYPE, intent(out) :: pd3d(I3DFIELD)
   REALTYPE, intent(out) :: pd2d(I2DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Hannes Rennau
!
! !LOCAL VARIABLES:
   REALTYPE                  :: dupper,dlower
   integer                   :: i,j,k
   REALTYPE                  :: aux(I3DFIELD)
!EOP
!-----------------------------------------------------------------------
!BOC
   pd3d=_ZERO_
   pd2d=_ZERO_
   if (Am .gt. _ZERO_) then
      ! 0.5*(dv/dx + du/dy)**2 on X-points
      do k=1,kmax
         do j=jmin-1,jmax
            do i=imin-1,imax
               if (ax(i,j).gt.0) then
                  aux(i,j,k)=                                                  &
                     _HALF_*(                                                  &
                      (vv(i+1,j,k)/hvn(i+1,j,k)-vv(i,j,k)/hvn(i,j,k))/DXX      &
                     +(uu(i,j+1,k)/hun(i,j+1,k)-uu(i,j,k)/hun(i,j,k))/DYX      &
                                                                     )**2               
               else 
                  aux(i,j,k)=_ZERO_
               end if
            end do
         end do
      end do 
      ! 2*Am*((du/dx)**2+0.5*(dv/dx + du/dy)**2+(dv/dy)**2) on T-points
      do k=1,kmax
         do j=jmin,jmax
            do i=imin,imax
               if (az(i,j).eq.1) then
                  pd3d(i,j,k)=_TWO_*dry_z(i,j)*Am*(                            &
                     _QUART_*(                                                 &
                           aux(i,j,k)+aux(i-1,j,k)+aux(i,j-1,k)+aux(i-1,j-1,k))&
                    +((uu(i,j,k)/hun(i,j,k)-uu(i-1,j,k)/hun(i-1,j,k))/DXC)**2  &
                    +((vv(i,j,k)/hvn(i,j,k)-vv(i,j-1,k)/hvn(i,j-1,k))/DYC)**2)
               end if
            end do
         end do
      end do 
   end if
   ! Av * ( (du/dz)**2 + (dv/dz)**2 ) on W-POINTS
   aux(:,:,kmax)=_ZERO_
   aux(:,:,0)   =_ZERO_
   do k=1,kmax-1
      do j=jmin,jmax
         do i=imin,imax
            if (az(i,j).eq.1) then
               aux(i,j,k)=(num(i,j,k)+avmmol)*SS(i,j,k)
            end if
         end do
      end do
   end do 
   ! Add Av * ( (du/dz)**2 + (dv/dz)**2 ) on T-POINTS
   do k=1,kmax
      do j=jmin,jmax
         do i=imin,imax
            if (az(i,j).eq.1) then
               pd3d(i,j,k)=pd3d(i,j,k)+_HALF_*(aux(i,j,k-1)+aux(i,j,k))
               pd2d(i,j)=pd2d(i,j)+pd3d(i,j,k)*hn(i,j,k)
            end if
         end do
      end do
   end do 

   return
   end subroutine physical_dissipation_3d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2010 - Hannes Rennau, Richard Hofmeister               !
!-----------------------------------------------------------------------
