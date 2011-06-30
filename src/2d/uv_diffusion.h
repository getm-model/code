   interface
!     explicit interface needed due to optional keyword arguments
      subroutine uv_diffusion(Am_method,An_method,UEx,VEx,         &
                              U,V,dudxC,dvdyC,shearX,D,DU,DV,Am,AmX)
         use domain, only: imin,imax,jmin,jmax,az,au,av,ax
#if defined(SPHERICAL) || defined(CURVILINEAR)
         use domain, only: dyc,arud1,dxx,dyx,arvd1,dxc
#else
         use domain, only: dx,dy,ard1
#endif
         use m2d, only: Am_const
         use variables_2d, only: PP,An,AnX
         use getm_timers,  only: tic,toc,TIM_UVDIFFUS
!$       use omp_lib
         IMPLICIT NONE
         integer,intent(in)                               :: Am_method,An_method
         REALTYPE,dimension(E2DFIELD),intent(in),optional :: U,V
         REALTYPE,dimension(E2DFIELD),intent(in),optional :: dudxC,dvdyC,shearX
         REALTYPE,dimension(E2DFIELD),intent(in),optional :: D,DU,DV
         REALTYPE,dimension(E2DFIELD),intent(in),optional :: Am,AmX
         REALTYPE,dimension(E2DFIELD),intent(inout)       :: UEx,VEx
      end subroutine uv_diffusion
   end interface