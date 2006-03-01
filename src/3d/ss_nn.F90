!$Id: ss_nn.F90,v 1.6 2006-03-01 14:45:12 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ss_nn - calculates shear and buoyancy frequency
! \label{sec-ss-nn}
!
! !INTERFACE:
   subroutine ss_nn()
!
! !DESCRIPTION:
!
! Here, the shear frequency squared, 
! $M^2=\left(\partial_z u\right)^2+\left(\partial_z v\right)^2$, 
! and the buoyancy frequency squared, $N^2=\partial_z b$, with buoyancy $b$ from
! (\ref{bdef}) are calculated. 
! For both calculations, two alternative methods are coded.
! The two straight-forward methods which are explained first, do
! both have the disadvantage of generating numerical instabilities.
! The straight-forward way for calculating $M^2$ is as follows:
!
! \begin{equation}\label{ShearSquaredOld}
! \begin{array}{l}
! \displaystyle
! (M^2)_{i,j,k}\approx
! \frac12
! \Bigg(\left(\frac{u_{i,j,k+1}-u_{i,j,k}}
! {\frac12(h^u_{i,j,k+1}+h^u_{i,j,k})}\right)^2
! +
! \left(\frac{u_{i-1,j,k+1}-u_{i-1,j,k}}
! {\frac12(h^u_{i-1,j,k+1}+h^u_{i-1,j,k})}\right)^2
! \\ \\ \qquad\qquad\quad
! \displaystyle
! +
! \left(\frac{v_{i,j,k+1}-v_{i,j,k}}
! {\frac12(h^v_{i,j,k+1}+h^v_{i,j,k})}\right)^2
! +
! \left(\frac{v_{i,j-1,k+1}-v_{i,j-1,k}}
! {\frac12(h^v_{i,j-1,k+1}+h^v_{i,j-1,k})}\right)^2
! \Bigg)
! \end{array}
! \end{equation}
!
! \cite{BURCHARD01c} developed a new scheme, which guarantees that
! the mean kinetic energy which is dissipated from the mean flow
! equals the shear production of turbulent kinetic energy. Therefore,
! this scheme should be numerically more stable than (\ref{ShearSquaredOld}):
!
! \begin{equation}\label{ShearSquaredNew}
! \begin{array}{l}
! \displaystyle
! (M^2)_{i,j,k}\approx
! \frac12
! \Bigg(\frac{\frac12(\nu_{i,j,k}+\nu_{i+1,j,k})
! (u_{i,j,k+1}-u_{i,j,k})^2}{\frac12(h^u_{i,j,k+1}+h^u_{i,j,k})}
! \\ \\ \qquad\qquad\quad
! \displaystyle
! +
! \frac{\frac12(\nu_{i-1,j,k}+\nu_{i,j,k})
! (u_{i-1,j,k+1}-u_{i-1,j,k})^2}{\frac12(h^u_{i-1,j,k+1}+h^u_{i-1,j,k})}
! \\ \\ \qquad\qquad\quad
! \displaystyle
! +
! \frac{\frac12(\nu_{i,j,k}+\nu_{i,j+1,k})
! (v_{i,j,k+1}-v_{i,j,k})^2}{\frac12(h^v_{i,j,k+1}+h^v_{i,j,k})}
! \\ \\ \qquad\qquad\quad
! \displaystyle
! +
! \frac{\frac12(\nu_{i,j-1,k}+\nu_{i,j,k})
! (v_{i,j-1,k+1}-v_{i,j-1,k})^2}{\frac12(h^v_{i,j-1,k+1}+h^v_{i,j-1,k})}
! \Bigg)
! \\ \\ \qquad\qquad\quad
! \displaystyle
! \cdot
! \left(\frac12\left(h^c_{i,j,k}+h^c_{i,j,k+1}\right)\nu_{i,j,k}\right)^{-1}
! \end{array}
! \end{equation}
! 
! The straight-forward discretisation of
! $N^2$ is given by
! 
! \begin{equation}\label{Nstraight}
! \begin{array}{l}
! \displaystyle
! (N^2)_{i,j,k}\approx
! \frac{b_{i,j,k+1}-b_{i,j,k}}{\frac12(h^t_{i,j,k+1}+h^t_{i,j,k})}.
! \end{array}
! \end{equation}
! 
! In some cases, together with the straight-forward discretisation
! of the shear squared, (\ref{ShearSquaredOld}), this
! did not produce stable numerical results. The reason for this might be that
! the velocities involved in the calculation for the shear squared do depend
! on the buoyancies in the two
! neighbouring T-points such that the straight-forward
! method (\ref{Nstraight}) leads to an inconsistency.
! However, other experiments with the energy-conserving discretisation
! of the shear stress squared, (\ref{ShearSquaredNew})
! and the straight-forward discretisation of
! $N^2$, (\ref{Nstraight}),  produced numerically stable results.
!
! Most stable results have been obtained with a weighted average 
! for the $N^2$ calculation:
!
! \begin{equation}\label{Naveraged}
! \begin{array}{l}
! \displaystyle
! (N^2)_{i,j,k}\approx
! \frac16 \Bigg(
! 2\frac{b_{i,j,k+1}-b_{i,j,k}}{\frac12(h^t_{i,j,k+1}+h^t_{i,j,k})}
! \\ \\ \qquad\qquad\qquad
! \displaystyle
! +
! \frac{b_{i+1,j,k+1}-b_{i+1,j,k}}{\frac12(h^t_{i+1,j,k+1}+h^t_{i+1,j,k})}
! +
! \frac{b_{i-1,j,k+1}-b_{i-1,j,k}}{\frac12(h^t_{i-1,j,k+1}+h^t_{i-1,j,k})}
! \\ \\ \qquad\qquad\qquad
! \displaystyle
! +
! \frac{b_{i,j+1,k+1}-b_{i,j+1,k}}{\frac12(h^t_{i,j+1,k+1}+h^t_{i,j+1,k})}
! +
! \frac{b_{i,j-1,k+1}-b_{i,j-1,k}}{\frac12(h^t_{i,j-1,k+1}+h^t_{i,j-1,k})}
! \Bigg).
! \end{array}
! \end{equation}
!
! These stability issues need to be further investigated in the future. 
!
! !USES:
   use domain, only: iimin,iimax,jjmin,jjmax,kmax,au,av,az
   use variables_3d, only: kmin,kumin,hn,uu,hun,kvmin,vv,hvn,SS,num
#ifndef NO_BAROCLINIC
   use variables_3d, only: NN,rho
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,nb
   REALTYPE                  :: dz,NNc,NNe,NNw,NNn,NNs
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'ss_nn() # ',Ncall
#endif

#undef NEW_SS
#define NEW_SS

!  Prandtl frequency
   do j=jjmin,jjmax
      do i=iimin,iimax
         if (az(i,j) .eq. 1 ) then
            do k=1,kmax-1
! This is an older version which we should keep here.
#ifndef NEW_SS
              SS(i,j,k)=0.5* (                                               &
                   ( (uu(i,j,k+1)/hun(i,j,k+1)-uu(i,j,k)/hun(i,j,k))         &
                   /(0.5*(hun(i,j,k+1)+hun(i,j,k))) )**2                     &
                +  ( (uu(i-1,j,k+1)/hun(i-1,j,k+1)-uu(i-1,j,k)/hun(i-1,j,k)) &
                   /(0.5*(hun(i-1,j,k+1)+hun(i-1,j,k))) )**2                 &
                +  ( (vv(i,j,k+1)/hvn(i,j,k+1)-vv(i,j,k)/hvn(i,j,k))         &
                   /(0.5*(hvn(i,j,k+1)+hvn(i,j,k))) )**2                     &
                +  ( (vv(i,j-1,k+1)/hvn(i,j-1,k+1)-vv(i,j-1,k)/hvn(i,j-1,k)) &
                   /(0.5*(hvn(i,j-1,k+1)+hvn(i,j-1,k))) )**2                 &
                            )
#else
! This version should better conserve energy.
              SS(i,j,k)=0.5* (                                                &
                   (uu(i,j,k+1)/hun(i,j,k+1)-uu(i,j,k)/hun(i,j,k))**2         &
                   /(0.5*(hun(i,j,k+1)+hun(i,j,k)))                           &
                    *0.5*(num(i,j,k)+num(i+1,j,k))                            &
               +  (uu(i-1,j,k+1)/hun(i-1,j,k+1)-uu(i-1,j,k)/hun(i-1,j,k))**2 &
                  /(0.5*(hun(i-1,j,k+1)+hun(i-1,j,k)))                       &
                   *0.5*(num(i-1,j,k)+num(i,j,k))                            &
                +  (vv(i,j,k+1)/hvn(i,j,k+1)-vv(i,j,k)/hvn(i,j,k))**2         &
                   /(0.5*(hvn(i,j,k+1)+hvn(i,j,k)))                           &
                    *0.5*(num(i,j,k)+num(i,j+1,k))                            &
                +  (vv(i,j-1,k+1)/hvn(i,j-1,k+1)-vv(i,j-1,k)/hvn(i,j-1,k))**2 &
                   /(0.5*(hvn(i,j-1,k+1)+hvn(i,j-1,k)))                       &
                    *0.5*(num(i,j-1,k)+num(i,j,k))                            &
                            )/(0.5*(hn(i,j,k)+hn(i,j,k+1)))/num(i,j,k)
#endif
            end do
         end if
      end do
   end do

#ifndef NO_BAROCLINIC
#define NEW_NN
#undef NEW_NN
   do j=jjmin,jjmax
      do i=iimin,iimax
         if (az(i,j) .eq. 1 ) then
            do k=kmax-1,1,-1
               dz=0.5*(hn(i,j,k+1)+hn(i,j,k))
               NNc =(rho(i,j,k+1)-rho(i,j,k))/dz
#ifndef NEW_NN
               if (az(i+1,j) .eq. 1) then
                  dz=0.5*(hn(i+1,j,k+1)+hn(i+1,j,k))
                  NNe=(rho(i+1,j,k+1)-rho(i+1,j,k))/dz
               else
                  NNe=NNc
               end if 
               if (az(i-1,j) .eq. 1) then
                  dz=0.5*(hn(i-1,j,k+1)+hn(i-1,j,k))
                  NNw=(rho(i-1,j,k+1)-rho(i-1,j,k))/dz
               else
                  NNw=NNc
               end if
               if (az(i,j+1) .eq. 1) then
                  dz=0.5*(hn(i,j+1,k+1)+hn(i,j+1,k))
                  NNn=(rho(i,j+1,k+1)-rho(i,j+1,k))/dz
               else
                  NNn=NNc
               end if
               if (az(i,j-1) .eq. 1) then
                  dz=0.5*(hn(i,j-1,k+1)+hn(i,j-1,k))
                  NNs=(rho(i,j-1,k+1)-rho(i,j-1,k))/dz
               else
                  NNs=NNc
               end if
               NN(i,j,k)=0.3333333*NNc+0.1666666*(NNe+NNw+NNn+NNs)
#else
               NN(i,j,k)=NNc
#endif
            end do
         end if
      end do
   end do
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving ss_nn()'
   write(debug,*)
#endif
   return
   end subroutine ss_nn
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
