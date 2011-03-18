#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: uv_depths - calculate depths in u and v points.\label{sec-uv-depth}
!
! !INTERFACE:
   subroutine uv_depths(vel_depth_method)
!
! !DESCRIPTION:
!
! In this routine which is called once during the model initialisation,
! the bathymetry value in the U- and the V-points are calculated from the
! bathymetry values in the T-points. The interpolation depends on the value
! which is given to {\tt vel\_depth\_method}:
!
! \begin{equation}
! H^u_{i,j} = \left\{
! \begin{array}{ll}
! \displaystyle
! \frac12 \left(H_{i,j}+H_{i+1,j}\right), & 
! \displaystyle
! \mbox{ for {\tt vel\_depth\_method}} =0, \\ \\ 
! \displaystyle
! \min\left\{H_{i,j}+H_{i+1,j}\right\}, & 
! \displaystyle
! \mbox{ for {\tt vel\_depth\_method}} =1, \\ \\ 
! \displaystyle
! \min\left\{H_{i,j}+H_{i+1,j}\right\}, & 
! \displaystyle
! \mbox{ for {\tt vel\_depth\_method}} =2 \mbox{ and } \min\{H_{i,j}i,H_{i+1,j}\}<D_{crit} \\ \\ 
! \displaystyle
! \frac12 \left(H_{i,j}+H_{i+1,j}\right), & 
! \displaystyle
! \mbox{ for {\tt vel\_depth\_method}} =2 \mbox{ and } \min\{H_{i,j},H_{i+1,j}\}\geq D_{crit} \\ \\ 
! \end{array}
! \right.
! \end{equation}
!
! The calculation of $H^v_{i,j}$ is done accordingly.
!
! The options 1 and 2 for {\tt vel\_depth\_method} may help to stabilise
! calculations when drying and flooding is involved.
!
! !USES:
   use exceptions
   use domain, only: imin,imax,jmin,jmax,az,au,av,H,HU,HV
   use variables_2d, only: DU,DV
   use getm_timers, only: tic,toc,TIM_UVDEPTHS
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: vel_depth_method
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: i,j
   REALTYPE                  :: d_crit=2.0 
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(0,*) 'uv_depths() # ',Ncall
#endif
   CALL tic(TIM_UVDEPTHS)

   do j=jmin-HALO,jmax+HALO
      do i=imin-HALO,imax+HALO-1
         select case (vel_depth_method)
            case (0)
               HU(i,j)=0.5*(H(i,j)+H(i+1,j))
            case (1)
               HU(i,j)=min(H(i,j),H(i+1,j))
            case (2)
               if (H(i,j) .lt. d_crit .or. H(i+1,j) .lt. d_crit) then
                  HU(i,j)=min(H(i,j),H(i+1,j))
               else
                  HU(i,j)=0.5*(H(i,j)+H(i+1,j))
               end if
            case default
               call getm_error("uv_depths()", &
                               "vel_depth_method must be 0, 1 or 2")
         end select
      end do
   end do

   do j=jmin-HALO,jmax+HALO-1
      do i=imin-HALO,imax+HALO
         select case (vel_depth_method)
            case (0)
               HV(i,j)=0.5*(H(i,j)+H(i,j+1))
            case (1)
               HV(i,j)=min(H(i,j),H(i,j+1))
            case (2)
               if (H(i,j) .lt. d_crit .or. H(i,j+1) .lt. d_crit) then
                  HV(i,j)=min(H(i,j),H(i,j+1))
               else
                  HV(i,j)=0.5*(H(i,j)+H(i,j+1))
               end if
            case default
               call getm_error("uv_depths()", &
                               "vel_depth_method must be 0, 1 or 2")
         end select
      end do
   end do

   CALL toc(TIM_UVDEPTHS)
#ifdef DEBUG
   write(0,*) 'Leaving uv_depths()'
   write(0,*)
#endif
   return
   end subroutine uv_depths
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
