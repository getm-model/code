!$Id: stresses_3d.F90,v 1.9 2009-04-22 10:02:36 lars Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: stresses_3d - bottom and surface stresses\label{sec-stresses-3d}
!
! !INTERFACE:
   subroutine stresses_3d
!
! !DESCRIPTION:
!
! As preparation of the call to {\tt do\_turbulence} in the routine {\tt gotm},
! see section \ref{sec-gotm}, the normalised surface and bottom stresses,
! $\tau_s/\rho_0$ (variable {\tt taus}) and $\tau_b/\rho_0$ (variable {\tt
! taub}), respectively, are calculated and interpolated to the T-points.
! Input parameters to this routine are {\tt rru} and {tt rrv}, which
! contain $r\sqrt{u^2+v^2}$ for the U- and V-points, respectively. 
! The modules of the surface and bottom stress vectors are calculated
! then by means of taking the square root of the sum of the squares of
! the stess components. In a similar way also the $x$- and $y$-components
! of the bottom stress are computed for output. 
!
! !USES:
   use parameters, only: rho_0
   use domain, only: az,au,av,imin,imax,jmin,jmax
   use variables_3d, only: kumin,kvmin,uu,vv,hun,hvn,rru,rrv
   use variables_3d, only: taus,taubx,tauby,taub
   use meteo, only: tausx,tausy
   use halo_zones, only : update_2d_halo,wait_halo,z_TAG
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,k1,k2,k3,k4
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'stresses_3d() # ',Ncall
#endif

!  we need to know rru and rrv in the halos as well
   call update_2d_halo(rru,rru,au,imin,jmin,imax,jmax,10)
   call wait_halo(10)
   call update_2d_halo(rrv,rrv,av,imin,jmin,imax,jmax,10)
   call wait_halo(10)


!  x-component of bottom momentum flux at U-points
   do j=jmin,jmax    
      do i=imin-1,imax

         k          = kumin(i,j)                      ! bottom index
         taubx(i,j) = - uu(i,j,k)/hun(i,j,k)*rru(i,j) ! momentum flux

      enddo
   enddo

!  y-component of bottom momentum flux at V-points  
   do j=jmin-1,jmax    
      do i=imin,imax
 
         k          = kvmin(i,j)                      ! bottom index
         tauby(i,j) = - vv(i,j,k)/hvn(i,j,k)*rrv(i,j) ! momentum flux

      enddo
   enddo

!  stress magnitude
   do j=jmin,jmax    
      do i=imin,imax

!        lower indices at U- and V-points
         k1=kumin(i-1,j  )
         k2=kumin(i  ,j  )
         k3=kvmin(i  ,j-1)
         k4=kvmin(i  ,j  )


!        total bottom stress at T-points
         taub(i,j)=0.5*(                                            &
                   (uu(i-1,j  ,k1)/hun(i-1,j  ,k1)*rru(i-1,j  ))**2 &
                  +(uu(i  ,j  ,k2)/hun(i  ,j  ,k2)*rru(i  ,j  ))**2 &
                  +(vv(i  ,j-1,k3)/hvn(i  ,j-1,k3)*rrv(i  ,j-1))**2 &
                  +(vv(i  ,j  ,k4)/hvn(i  ,j  ,k4)*rrv(i  ,j  ))**2)

         taub(i,j)=sqrt(taub(i,j))

!        total surface stress at T-points
         taus(i,j)=0.5*(                                &
                        tausx(i,j)**2+tausx(i-1,j)**2   &
                      + tausy(i,j)**2+tausy(i,j-1)**2)
         taus(i,j)=sqrt(taus(i,j))/rho_0

      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving stresses_3d()'
   write(debug,*)
#endif
   return
   end subroutine stresses_3d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
