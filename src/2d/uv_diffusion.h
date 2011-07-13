   interface
!     explicit interface needed due to optional keyword arguments
      subroutine uv_diffusion(An_method,UEx,VEx,U,V,D,DU,DV, &
                              dudxC,dvdyC,shearX,AmC,AmX)
         use domain, only: imin,imax,jmin,jmax,az,au,av,ax
#if defined(SPHERICAL) || defined(CURVILINEAR)
         use domain, only: dyc,arud1,dxx,dyx,arvd1,dxc
#else
         use domain, only: dx,dy,ard1
#endif
         use m2d, only: Am_method,Am_const,AM_CONSTANT,AM_LES,An_const
         use variables_2d, only: PP,AnC,AnX
         use getm_timers,  only: tic,toc,TIM_UVDIFFUS
!$       use omp_lib
         IMPLICIT NONE
         integer,intent(in)                               :: An_method
         REALTYPE,dimension(E2DFIELD),intent(in),optional :: U,V,D,DU,DV
         REALTYPE,dimension(E2DFIELD),intent(in),optional :: dudxC,dvdyC,shearX
         REALTYPE,dimension(E2DFIELD),intent(in),optional :: AmC,AmX
         REALTYPE,dimension(E2DFIELD),intent(inout)       :: UEx,VEx
      end subroutine uv_diffusion
   end interface
