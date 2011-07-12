   interface
!     explicit interface needed due to optional keyword arguments
      subroutine deformation_rates(U,V,DU,DV,dudxC,dudxV,dudxU,  &
                                             dvdyC,dvdyU,dvdyV,  &
                                             dudyX,dvdxX,        &
                                             shearX,shearU,shearV)
         use domain, only: imin,imax,jmin,jmax,az,au,av,ax
#if defined(SPHERICAL) || defined(CURVILINEAR)
         use domain, only: dxc,dyc,dxu,dyu,dxv,dyv,dxx,dyx,arcd1,arud1,arvd1
#else
         use domain, only: dx,dy,ard1
#endif
         use domain, only: NWB,NNB,NEB,NSB
         use domain, only: wi,wfj,wlj,nj,nfi,nli,ei,efj,elj,sj,sfi,sli
         use getm_timers, only: tic,toc,TIM_DEFORM
         IMPLICIT NONE
         REALTYPE,dimension(E2DFIELD),intent(in)           :: U,V,DU,DV
         REALTYPE,dimension(E2DFIELD),intent(out),optional :: dudxC,dudxV,dudxU
         REALTYPE,dimension(E2DFIELD),intent(out),optional :: dvdyC,dvdyU,dvdyV
         REALTYPE,dimension(E2DFIELD),intent(out),optional :: dudyX,dvdxX
         REALTYPE,dimension(E2DFIELD),intent(out),optional :: shearX,shearU,shearV
      end subroutine deformation_rates
   end interface
