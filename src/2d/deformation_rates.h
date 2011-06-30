   interface
!     explicit interface needed due to optional keyword arguments
      subroutine deformation_rates(U,V,DU,DV,dudxC,dvdyC,&
                                             dudyX,dvdxX,shearX,&
                                             dudxU,dvdyU,shearU,&
                                             dudxV,dvdyV,shearV )
         use domain, only: imin,imax,jmin,jmax,az,au,av,ax
#if defined(SPHERICAL) || defined(CURVILINEAR)
         use domain, only: dxc,dyc,dxu,dyu,dxv,dyv,dxx,dyx
#else
         use domain, only: dx,dy
#endif
         IMPLICIT NONE
         REALTYPE,dimension(E2DFIELD),intent(in)           :: U,V,DU,DV
         REALTYPE,dimension(E2DFIELD),intent(out)          :: dudxC,dvdyC
         REALTYPE,dimension(E2DFIELD),intent(out),optional :: dudyX,dvdxX,shearX
         REALTYPE,dimension(E2DFIELD),intent(out),optional :: dudxU,dvdyU,shearU
         REALTYPE,dimension(E2DFIELD),intent(out),optional :: dudxV,dvdyV,shearV
      end subroutine deformation_rates
   end interface
