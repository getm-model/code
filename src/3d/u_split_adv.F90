!$Id: u_split_adv.F90,v 1.1 2004-01-06 15:04:00 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:  u_split_adv()
!
! !INTERFACE:
   subroutine u_split_adv(dt,f,uu,hun, &
                          delxu,delyu,area_inv,au,splitfac,method,az,AH)
! !DESCRIPTION:
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
   use domain, only: iimin,iimax,jjmin,jjmax,kmax
   use advection_3d, only: hi,hio,cu
   use advection_3d, only: UPSTREAM_SPLIT,P2,SUPERBEE,MUSCL,P2_PDM
   use advection_3d, only: ONE,TWO,ONE6TH
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)      :: uu(I3DFIELD),hun(I3DFIELD)
   REALTYPE, intent(in)      :: delxu(I2DFIELD),delyu(I2DFIELD)
   REALTYPE, intent(in)      :: area_inv(I2DFIELD),dt
   integer, intent(in)       :: au(E2DFIELD), az(E2DFIELD)
   REALTYPE, intent(in)      :: splitfac
   integer, intent(in)       :: method
   REALTYPE, intent(in)      :: AH
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout)   :: f(I3DFIELD)
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: u_split_adv.F90,v $
!  Revision 1.1  2004-01-06 15:04:00  kbk
!  FCT advection + split of advection_3d.F90 + extra adv. input checks
!
!
! !LOCAL VARIABLES:
   integer         :: i,ii,j,jj,k,kk
   REALTYPE        :: c,alpha,beta,x,r,Phi,limit,fu,fc,fd
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'u_split_adv() # ',Ncall
#endif

   cu = _ZERO_

! Calculating u-interface fluxes !
   select case (method)
      case (UPSTREAM_SPLIT)
         do k=1,kmax
            do j=jjmin,jjmax
               do i=iimin-1,iimax
                  cu(i,j,k) = _ZERO_
                  if (au(i,j) .gt. 0) then
                     if (uu(i,j,k) .gt. _ZERO_) then
                        cu(i,j,k)=uu(i,j,k)*f(i,j,k)
                     else
                        cu(i,j,k)=uu(i,j,k)*f(i+1,j,k)
                     end if
                  end if
               end do
            end do
         end do
      case ((P2),(Superbee),(MUSCL),(P2_PDM))
         do k=1,kmax
            do j=jjmin,jjmax
               do i=iimin-1,iimax
                  cu(i,j,k) = _ZERO_
                  if (au(i,j) .gt. 0) then
                     if (uu(i,j,k) .gt. _ZERO_) then
                        c=uu(i,j,k)/hun(i,j,k)*dt/delxu(i,j)
                        if (au(i-1,j) .gt. 0) then
                           fu=f(i-1,j,k)         ! upstream
                        else
                           fu=f(i  ,j,k)
                        end if
                        fc=f(i  ,j,k)            ! central
                        fd=f(i+1,j,k)            ! downstream
                        if (abs(fd-fc) .gt. 1e-10) then
                           r=(fc-fu)/(fd-fc)     ! slope ratio
                        else
                           r=(fc-fu)*1.e10
                        end if
                     else
                        c=-uu(i,j,k)/hun(i,j,k)*dt/delxu(i,j)
                        if (au(i+1,j) .gt. 0) then
                           fu=f(i+2,j,k)         ! upstream
                        else
                           fu=f(i+1,j,k)
                        end if
                        fc=f(i+1,j,k)            ! central
                        fd=f(i  ,j,k)            ! downstream
                        if (abs(fc-fd) .gt. 1e-10) then
                           r=(fu-fc)/(fc-fd)     ! slope ratio
                        else
                           r=(fu-fc)*1.e10
                        end if
                     end if
                     select case (method)
                        case ((P2),(P2_PDM))
                           x = one6th*(ONE-TWO*c)
                           Phi=(0.5+x)+(0.5-x)*r
                           if (method.eq.P2) then
                              limit=Phi
                           else
                              limit=max(_ZERO_,min(Phi,2./(1.-c), &
                                        2.*r/(c+1.e-10)))
                           end if
                        case (Superbee)
                           limit=max(_ZERO_,min(ONE,TWO*r),min(r,TWO))
                        case (MUSCL)
                           limit=max(_ZERO_,min(TWO,TWO*r,0.5*(ONE+r)))
                        case default
                           FATAL 'Not so good - do_advection_3d()'
                           stop 'u_split_adv'
                     end select
                     cu(i,j,k)=uu(i,j,k)*(fc+0.5*limit*(1-c)*(fd-fc))
!Horizontal diffusion
                     if ( AH.gt.0. .and. az(i,j).gt.0 .and. az(i+1,j).gt.0 ) then
                        cu(i,j,k) = cu(i,j,k)-AH*hun(i,j,k) &
                                     *(f(i+1,j,k)-f(i,j,k))/delxu(i,j)
                     end if
                  end if
               end do
            end do
         end do
   end select

   do k=1,kmax   ! Doing the u-advection step
      do j=jjmin,jjmax
         do i=iimin,iimax
            if (az(i,j) .eq. 1) then
               hio(i,j,k)=hi(i,j,k)
               hi(i,j,k)=hio(i,j,k)                           &
                        -splitfac*dt*(uu(i,j,k)*delyu(i,j)    &
                                     -uu(i-1,j,k)*delyu(i-1,j))*area_inv(i,j)
               f(i,j,k)=(f(i,j,k)*hio(i,j,k)           &
                        -splitfac*dt*(cu(i,j,k)*delyu(i,j)    &
                        -cu(i-1,j,k)*delyu(i-1,j))*area_inv(i,j) &
                        )/hi(i,j,k)
            end if
         end do
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving u_split_adv()'
   write(debug,*)
#endif
   return
   end subroutine u_split_adv

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
