!$id: advection_3d.F90,v 1.18 2001/09/19 13:53:08 bbh Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  3D advection
!
! !INTERFACE:
   module advection_3d
!
! !DESCRIPTION:
!  This module do advection of scalars.  The module follows the same
!  convention as the other modules in 'getm'. The module is initialised
!  by calling 'init\_advection\_3d()'. In the time-loop 'do\_advection\_3d' is
!  called. 'do\_advection\_3d' is a wrapper routine which - dependent on the
!  actual advection scheme chosen - makes a call to the appropriate subroutine.
!  New advection schemes are easily implemented - at least from a program
!  point of view - since only this module needs to be changed.
!  Additional work arrays can easily be added following the stencil given
!  below. To add a new advection scheme 3 things must be done: 1) define
!  a unique constant to identify the scheme (see e.g. UPSTREAM and TVD)
!  2) adopt the 'select case' in 'do\_advection\_3d()' and 3) write the actual
!  subroutine.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
   use domain, only: iimin,iimax,jjmin,jjmax,kmax

   IMPLICIT NONE
!
   private
!
! !PUBLIC DATA MEMBERS:
   public init_advection_3d, do_advection_3d
!
! !PRIVATE DATA MEMBERS:
#ifdef STATIC
   REALTYPE	:: cu(I3DFIELD),hi(I3DFIELD),hio(I3DFIELD)
#else
   REALTYPE, dimension(:,:,:), allocatable	:: cu,hi,hio
#endif
   REALTYPE, parameter		:: one6th=1./6.
   integer, parameter 		:: UPSTREAM=1,UPSTREAM_SPLIT=2,P2=3
   integer, parameter 		:: Superbee=4,MUSCL=5,P2_PDM=6
   logical, parameter 		:: STRANG=.true.
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: advection_3d.F90,v $
!  Revision 1.1  2002-05-02 14:00:58  gotm
!  Initial revision
!
!  Revision 1.18  2001/09/19 13:53:08  bbh
!  Typo
!
!  Revision 1.17  2001/09/19 13:46:52  bbh
!  Nasty memory leak in upstream_adv() - only if -DFORTRAN90
!
!  Revision 1.16  2001/09/03 20:13:08  bbh
!  Corrected order of select (strang)
!
!  Revision 1.15  2001/09/03 20:04:21  bbh
!  Allow individual advection settings for momentum, salinity and temperature
!
!  Revision 1.14  2001/09/03 12:56:58  bbh
!  Advection can now be split into different schemes for each direction
!
!  Revision 1.13  2001/08/29 14:22:59  bbh
!  Dimensions for masks are needed
!
!  Revision 1.12  2001/08/29 12:44:39  bbh
!  masks are E2DFIELD
!
!  Revision 1.11  2001/08/27 11:51:45  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.10  2001/08/01 08:31:22  bbh
!  CURVILINEAR now implemented
!
!  Revision 1.9  2001/07/26 12:49:33  bbh
!  Added a number of new advection schemes - UPSTREAM_SPLIT, P2,
!  Superbee, MUSCL and P2_PDM.
!
!  Revision 1.8  2001/05/22 08:24:21  bbh
!  Added adv in call to advection routines
!
!  Revision 1.7  2001/05/18 09:45:16  bbh
!  Removed some LEVEL2 statements
!
!  Revision 1.6  2001/05/18 09:37:07  bbh
!  IMPLICITE NONE was commented out  at one place
!
!  Revision 1.5  2001/05/18 07:59:46  bbh
!  UPSTREAM advection is ready
!
!  Revision 1.4  2001/05/11 13:47:53  bbh
!  First do UPSTREAM
!
!  Revision 1.3  2001/05/04 13:15:05  bbh
!  Public members was commented out
!
!  Revision 1.2  2001/05/04 09:53:06  bbh
!  Added stubs for upstream_adv() and tdv_adv + description
!
!  Revision 1.1  2001/05/03 20:20:33  bbh
!  Stubs for baroclinicity
!
!
! !LOCAL VARIABLES:
   integer :: advection_method
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_advection_3d
!
! !INTERFACE:
   subroutine init_advection_3d(method)
!
! !DESCRIPTION:
!  Reads the namelist and makes calls to the init functions of the
!  various model components.
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)	:: method
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer	:: rc
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_advection_3d() # ',Ncall
#endif

   LEVEL1 'init_advection_3d()'

#ifdef STATIC
#else
   allocate(cu(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection_3d: Error allocating memory (cu)'

   allocate(hi(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection_3d: Error allocating memory (hi)'

   allocate(hio(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_advection_3d: Error allocating memory (hio)'
#endif

   cu  = _ZERO_
   hi  = _ZERO_
   hio = _ZERO_

#ifdef DEBUG
   write(debug,*) 'Leaving init_advection_3d()'
   write(debug,*)
#endif
   return
   end subroutine init_advection_3d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_advection_3d()
!
! !INTERFACE:
   subroutine do_advection_3d(dt,f,uu,vv,ww,hun,hvn,ho,hn,      &
                              delxu,delxv,delyu,delyv,area_inv,  &
			      az,au,av,hor_adv,ver_adv,strang,AH)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
  REALTYPE, intent(in)		:: uu(I3DFIELD),vv(I3DFIELD),ww(I3DFIELD)
  REALTYPE, intent(in)		:: hun(I3DFIELD),hvn(I3DFIELD)
  REALTYPE, intent(in)		:: ho(I3DFIELD),hn(I3DFIELD)
  REALTYPE, intent(in)		:: delxu(I2DFIELD),delxv(I2DFIELD)
  REALTYPE, intent(in)		:: delyu(I2DFIELD),delyv(I2DFIELD)
  REALTYPE, intent(in)		:: area_inv(I2DFIELD),dt,AH
  integer, intent(in)		:: az(E2DFIELD),au(E2DFIELD),av(E2DFIELD)
  integer, intent(in)		:: hor_adv,ver_adv,strang
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
  REALTYPE, intent(out)	:: f(I3DFIELD)
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   REALTYPE, parameter	:: a1=0.5,a2=1.0
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_advection_3d() # ',Ncall
#endif

   select case (hor_adv)
      case (UPSTREAM)
         call upstream_adv(dt,f,uu,vv,ww,ho,hn,delxv,delyu,delxu,delyv,area_inv,az,AH)
      case ((UPSTREAM_SPLIT),(P2),(Superbee),(MUSCL),(P2_PDM))
         hi=ho
         select case (strang)
            case (0)
               call u_split_adv(dt,f,uu,hun,delxu,delyu,area_inv,au,a2,   &
                                hor_adv,az,AH)
               call v_split_adv(dt,f,vv,hvn,delxv,delyv,area_inv,av,a2,   &
                                hor_adv,az,AH)
               if (kmax.gt.1) then
#ifdef ITERATE_VERT_ADV
                  call w_split_it_adv(dt,f,ww,az,a2,ver_adv)
#else
                  call w_split_adv(dt,f,ww,az,a2,ver_adv)
#endif
               end if
            case (1)
               call u_split_adv(dt,f,uu,hun,delxu,delyu,area_inv,au,a1,   &
                                hor_adv,az,AH)
               call v_split_adv(dt,f,vv,hvn,delxv,delyv,area_inv,av,a1,   &
                                hor_adv,az,AH)
               if (kmax.gt.1) then
#ifdef ITERATE_VERT_ADV
                  call w_split_it_adv(dt,f,ww,az,a2,ver_adv)
#else
                  call w_split_adv(dt,f,ww,az,a2,ver_adv)
#endif
	       end if
               call v_split_adv(dt,f,vv,hvn,delxv,delyv,area_inv,av,a1,   &
	                        hor_adv,az,AH)
               call u_split_adv(dt,f,uu,hun,delxu,delyu,area_inv,au,a1,   &
	                        hor_adv,az,AH)
           case default
              FATAL 'Not valid strang parameter'
         end select
      case default
         FATAL 'This is not so good - do_advection_3d()'
	 stop
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving do_advection_3d()'
   write(debug,*)
#endif
   return
   end subroutine do_advection_3d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  upstream_adv()
!
! !INTERFACE:
   subroutine upstream_adv(dt,f,uu,vv,ww,ho,hn,delxv,delyu,delxu,delyv,area_inv,az,AH)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE             :: uu(I3DFIELD),vv(I3DFIELD),ww(I3DFIELD)
   REALTYPE             :: ho(I3DFIELD),hn(I3DFIELD)
   REALTYPE		:: delxv(I2DFIELD),delyu(I2DFIELD)
   REALTYPE		:: delxu(I2DFIELD),delyv(I2DFIELD)
   REALTYPE		:: area_inv(I2DFIELD),dt,AH
   integer, intent(in)	:: az(E2DFIELD)
   REALTYPE             :: f(I3DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer	:: rc,i,ii,j,jj,k,kk
   REALTYPE, dimension(:,:,:), allocatable:: adv
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'upstream_adv() # ',Ncall
#endif

   allocate(adv(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'upstream_adv: Error allocating memory (adv)'

   cu(iimin-1:iimin-1,jjmin:jjmax,0:kmax)= _ZERO_
   cu(iimax:iimax,jjmin:jjmax,0:kmax)= _ZERO_

   do k=1,kmax   ! Calculating u-interface fluxes !
      do j=jjmin,jjmax
         do i=iimin,iimax-1
            if (uu(i,j,k) .gt. _ZERO_) then
               cu(i,j,k)=uu(i,j,k)*f(i,j,k)
	    else
               cu(i,j,k)=uu(i,j,k)*f(i+1,j,k)
	    end if
            if ((AH.gt.0.).and.(az(i,j).gt.0).and.(az(i+1,j).gt.0))  &
               cu(i,j,k)=cu(i,j,k)-AH*(f(i+1,j,k)-f(i,j,k))/delxu(i,j)  &
	                 *0.5*(hn(i+1,j,k)+hn(i,j,k))
         end do
      end do
   end do

   do k=1,kmax   ! Updating the advection term for u-advection !
      do j=jjmin,jjmax
         do i=iimin,iimax
	    adv(i,j,k)=(cu(i,j,k)*delyu(i,j)    &
	               -cu(i-1,j,k)*delyu(i-1,j))*area_inv(i,j)
         end do
      end do
   end do

   cu(iimin:iimax,jjmin-1:jjmin-1,0:kmax)= _ZERO_
   cu(iimin:iimax,jjmax:jjmax,0:kmax)= _ZERO_

   do k=1,kmax   ! Calculating v-interface fluxes !
      do j=jjmin,jjmax-1
         do i=iimin,iimax
            if (vv(i,j,k) .gt. _ZERO_) then
               cu(i,j,k)=vv(i,j,k)*f(i,j,k)
	    else
               cu(i,j,k)=vv(i,j,k)*f(i,j+1,k)
	    end if
            if ((AH.gt.0.).and.(az(i,j).gt.0).and.(az(i,j+1).gt.0))   &
            cu(i,j,k)=cu(i,j,k)-AH*(f(i,j+1,k)-f(i,j,k))/delyv(i,j)   &
	        *0.5*(hn(i,j+1,k)+hn(i,j,k))
         end do
      end do
   end do

   do k=1,kmax   ! Updating the advection term for v-advection !
      do j=jjmin,jjmax
         do i=iimin,iimax
	    adv(i,j,k)=adv(i,j,k)+(cu(i,j,k)*delxv(i,j)   &
	                          -cu(i,j-1,k)*delxv(i,j-1))*area_inv(i,j)
         end do
      end do
   end do

   cu(iimin:iimax,jjmin:jjmax,0)= _ZERO_
   cu(iimin:iimax,jjmin:jjmax,kmax)= _ZERO_

   if (kmax.gt.1) then
      do k=1,kmax-1   ! Calculating w-interface fluxes !
         do j=jjmin,jjmax
            do i=iimin,iimax
               if (ww(i,j,k) .gt. _ZERO_) then
                  cu(i,j,k)=ww(i,j,k)*f(i,j,k)
	       else
                  cu(i,j,k)=ww(i,j,k)*f(i,j,k+1)
	       end if
            end do
         end do
      end do

      do k=1,kmax   ! Updating the advection term for w-advection !
         do j=jjmin,jjmax
            do i=iimin,iimax
	       adv(i,j,k)=adv(i,j,k)+(cu(i,j,k)-cu(i,j,k-1))
            end do
         end do
      end do
   end if

   do k=1,kmax   ! Doing the full advection in one step
      do j=jjmin,jjmax
         do i=iimin,iimax
            if (az(i,j).eq.1)                                        &
	    f(i,j,k)=(f(i,j,k)*ho(i,j,k)-dt*adv(i,j,k))/hn(i,j,k)
         end do
      end do
   end do

#ifdef FORTRAN90
   deallocate(adv,stat=rc)    ! work array
   if (rc /= 0) stop 'upstream_adv: Error de-allocating memory (adv)'
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving upstream_adv()'
   write(debug,*)
#endif
   return
   end subroutine upstream_adv
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  u_split_adv()
!
! !INTERFACE:
   subroutine u_split_adv(dt,f,uu,hun,delxu,delyu,area_inv,au,splitfac,method,az,AH)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in) :: uu(I3DFIELD),hun(I3DFIELD)
   REALTYPE, intent(in)	:: delxu(I2DFIELD),delyu(I2DFIELD)
   REALTYPE, intent(in)	:: area_inv(I2DFIELD),dt
   integer, intent(in)	:: au(E2DFIELD), az(E2DFIELD)
   REALTYPE, intent(in)	:: splitfac
   integer, intent(in)	:: method
   REALTYPE, intent(in) :: AH
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)	:: f(I3DFIELD)
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer	:: i,ii,j,jj,k,kk
   REALTYPE	:: c,alpha,beta,x,r,Phi,limit,fu,fc,fd
!
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
                        x = one6th*(1.-2.0*c)
                        Phi=(0.5+x)+(0.5-x)*r
			if (method.eq.P2) then
			   limit=Phi
			else
		           limit=max(_ZERO_,min(Phi,2./(1.-c),2.*r/(c+1.e-10)))
			end if
		     case (Superbee)
		        limit=max(_ZERO_, min(1.0, 2.0*r), min(r,2.0) )
		     case (MUSCL) 	
		        limit=max(_ZERO_,min(2.0,2.0*r,0.5*(1.0+r)))
                     case default
                        FATAL 'This is not so good - do_advection_3d()'
	                stop
		     end select
	             cu(i,j,k)=uu(i,j,k)*(fc+0.5*limit*(1-c)*(fd-fc))
                     !Horizontal diffusion
                     if ((AH.gt.0.).and.(az(i,j).gt.0).and.(az(i+1,j).gt.0))   then
		        cu(i,j,k)=cu(i,j,k)-AH*hun(i,j,k)*(f(i+1,j,k)-f(i,j,k))/delxu(i,j)
		     end if
		  end if
               end do
            end do
         end do
   end select

   do k=1,kmax   ! Doing the u-advection step
      do j=jjmin,jjmax
         do i=iimin,iimax
	    if (az(i,j).eq.1) then
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
!BOP
!
! !IROUTINE:  v_split_adv()
!
! !INTERFACE:
   subroutine v_split_adv(dt,f,vv,hvn,delxv,delyv,area_inv,av,splitfac,method,az,AH)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE , intent(in) :: vv(I3DFIELD),hvn(I3DFIELD)
   REALTYPE, intent(in)	:: delxv(I2DFIELD),delyv(I2DFIELD)
   REALTYPE, intent(in)	:: area_inv(I2DFIELD),dt
   integer, intent(in)	:: av(E2DFIELD),az(E2DFIELD)
   REALTYPE, intent(in)	:: splitfac
   integer, intent(in)	:: method
   REALTYPE, intent(in) :: AH
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)	:: f(I3DFIELD)
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer	:: i,ii,j,jj,k,kk
   REALTYPE	:: c,alpha,beta,x,r,Phi,limit,fu,fc,fd
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'v_split_adv() # ',Ncall
#endif

   cu = _ZERO_

! Calculating v-interface fluxes !

   select case (method)
      case (UPSTREAM_SPLIT)
         do k=1,kmax
            do j=jjmin-1,jjmax
               do i=iimin,iimax
                  cu(i,j,k) = _ZERO_
                  if(av(i,j) .gt. 0) then
                     if (vv(i,j,k) .gt. _ZERO_) then
                        cu(i,j,k)=vv(i,j,k)*f(i,j,k)
	             else
                        cu(i,j,k)=vv(i,j,k)*f(i,j+1,k)
	             end if
                  end if
               end do
            end do
         end do
      case ((P2),(Superbee),(MUSCL),(P2_PDM))
         do k=1,kmax
            do j=jjmin-1,jjmax
               do i=iimin,iimax
                  cu(i,j,k) = _ZERO_
	          if (av(i,j) .gt. 0) then
                     if (vv(i,j,k) .gt. _ZERO_) then
                        c=vv(i,j,k)/hvn(i,j,k)*dt/delyv(i,j)
		        if (av(i,j-1) .gt. 0) then
		           fu=f(i,j-1,k)         ! upstream
			else
		           fu=f(i,j  ,k)
			end if
		        fc=f(i,j  ,k)            ! central
		        fd=f(i,j+1,k)            ! downstream
		        if (abs(fd-fc) .gt. 1e-10) then
		           r=(fc-fu)/(fd-fc)     ! slope ratio
		        else
		           r=(fc-fu)*1.e10
		        end if
		     else
                        c=-vv(i,j,k)/hvn(i,j,k)*dt/delyv(i,j)
		        if (av(i,j+1) .gt. 0) then
		           fu=f(i,j+2,k)         ! upstream
			else
		           fu=f(i,j+1,k)
			end if
		        fc=f(i,j+1,k)            ! central
		        fd=f(i,j  ,k)            ! downstream
		        if (abs(fc-fd) .gt. 1e-10) then
		           r=(fu-fc)/(fc-fd)     ! slope ratio
		        else
		           r=(fu-fc)*1.e10
		        end if
		     end if
                     x = one6th*(1.-2.0*c)
                     Phi=(0.5+x)+(0.5-x)*r
		     select case (method)
		     case ((P2),(P2_PDM))
                        x = one6th*(1.-2.0*c)
                        Phi=(0.5+x)+(0.5-x)*r
			if (method.eq.P2) then
			   limit=Phi
			else
		           limit=max(_ZERO_,min(Phi,2./(1.-c),2.*r/(c+1.e-10)))
			end if
		     case (Superbee)
		        limit=max(_ZERO_, min(1.0, 2.0*r), min(r,2.0) )
		     case (MUSCL) 	
		        limit=max(_ZERO_,min(2.0,2.0*r,0.5*(1.0+r)))
                     case default
                        FATAL 'This is not so good - do_advection_3d()'
	                stop
		     end select
                     cu(i,j,k)=vv(i,j,k)*(fc+0.5*limit*(1-c)*(fd-fc))
		     !Horizontal diffusion
		     if ((AH.gt.0.).and.(az(i,j).gt.0).and.(az(i,j+1).gt.0)) &
		        cu(i,j,k)=cu(i,j,k)-AH*hvn(i,j,k)*(f(i,j+1,k)-f(i,j,k))/delyv(i,j)
		  end if
               end do
            end do
         end do
   end select

   do k=1,kmax   ! Doing the v-advection step
      do j=jjmin,jjmax
         do i=iimin,iimax
	    if (az(i,j).eq.1) then
               hio(i,j,k)=hi(i,j,k)
               hi(i,j,k)=hio(i,j,k)              &
	             -splitfac*dt*(vv(i,j,k)*delxv(i,j)    &
		                  -vv(i,j-1,k)*delxv(i,j-1))*area_inv(i,j)
	       f(i,j,k)=(f(i,j,k)*hio(i,j,k)        &
                   -splitfac*dt*(cu(i,j,k)*delxv(i,j)    &
       	                        -cu(i,j-1,k)*delxv(i,j-1))*area_inv(i,j))/hi(i,j,k)
	    end if
         end do
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving v_split_adv()'
   write(debug,*)
#endif
   return
   end subroutine v_split_adv

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  w_split_adv()
!
! !INTERFACE:
   subroutine w_split_adv(dt,f,ww,az,splitfac,method)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE , intent(in) :: ww(I3DFIELD),dt
   integer , intent(in)  :: az(E2DFIELD)
   REALTYPE, intent(in)	:: splitfac
   integer, intent(in)	:: method
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)	:: f(I3DFIELD)
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer	:: i,ii,j,jj,k,kk
   REALTYPE	:: c,alpha,beta,x,r,Phi,limit,fu,fc,fd
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'w_split_adv() # ',Ncall
#endif

   cu = _ZERO_

! Calculating w-interface fluxes !

   select case (method)
      case (UPSTREAM_SPLIT)
         do k=1,kmax-1
            do j=jjmin,jjmax
               do i=iimin,iimax
                  cu(i,j,k) = _ZERO_
                  if (az(i,j) .eq. 1) then
                     if (ww(i,j,k) .gt. _ZERO_) then
                        cu(i,j,k)=ww(i,j,k)*f(i,j,k)
                     else
                        cu(i,j,k)=ww(i,j,k)*f(i,j,k+1)
	             end if
	          end if
               end do
            end do
         end do
      case ((P2),(Superbee),(MUSCL),(P2_PDM))
         do k=1,kmax-1
            do j=jjmin,jjmax
               do i=iimin,iimax
                  cu(i,j,k) = _ZERO_
	          if (az(i,j) .eq. 1) then
                     if (ww(i,j,k) .gt. _ZERO_) then
                        c=ww(i,j,k)*dt/(0.5*(hi(i,j,k)+hi(i,j,k+1)))
		        if (k .gt. 1) then
		           fu=f(i,j,k-1)         ! upstream
			else
		           fu=f(i,j,k)
			end if
		        fc=f(i,j,k  )            ! central
		        fd=f(i,j,k+1)            ! downstream
		        if (abs(fd-fc) .gt. 1e-10) then
		           r=(fc-fu)/(fd-fc)     ! slope ratio
		        else
		           r=(fc-fu)*1.e10
		        end if
		     else
                        c=-ww(i,j,k)*dt/(0.5*(hi(i,j,k)+hi(i,j,k+1)))
		        if (k .lt. kmax-1) then
		           fu=f(i,j,k+2)         ! upstream
			else
		           fu=f(i,j,k+1)
			end if
		        fc=f(i,j,k+1)            ! central
		        fd=f(i,j,k  )            ! downstream
		        if (abs(fc-fd) .gt. 1e-10) then
		           r=(fu-fc)/(fc-fd)     ! slope ratio
		        else
		           r=(fu-fc)*1.e10
		        end if
		     end if
                     x = one6th*(1.-2.0*c)
                     Phi=(0.5+x)+(0.5-x)*r
		     select case (method)
		     case ((P2),(P2_PDM))
                        x = one6th*(1.-2.0*c)
                        Phi=(0.5+x)+(0.5-x)*r
			if (method.eq.P2) then
			   limit=Phi
			else
		           limit=max(_ZERO_,min(Phi,2./(1.-c),2.*r/(c+1.e-10)))
			end if
		     case (Superbee)
		        limit=max(_ZERO_, min(1.0, 2.0*r), min(r,2.0) )
		     case (MUSCL) 	
		        limit=max(_ZERO_,min(2.0,2.0*r,0.5*(1.0+r)))
                     case default
                        FATAL 'This is not so good - do_advection_3d()'
	                stop
		     end select
	             cu(i,j,k)=ww(i,j,k)*(fc+0.5*limit*(1-c)*(fd-fc))
		  end if
               end do
            end do
         end do
   end select

   do k=1,kmax   ! Doing a w-advection step
      do j=jjmin,jjmax
         do i=iimin,iimax
	    if (az(i,j) .eq. 1) then
               hio(i,j,k)=hi(i,j,k)
               hi(i,j,k)=hio(i,j,k)-splitfac*dt*(ww(i,j,k)-ww(i,j,k-1))
               f(i,j,k)=(f(i,j,k)*hio(i,j,k)-        &
                         splitfac*dt*(cu(i,j,k)-cu(i,j,k-1)))/hi(i,j,k)
            end if
         end do
      end do
   end do


#ifdef DEBUG
   write(debug,*) 'Leaving w_split_adv()'
   write(debug,*)
#endif
   return
   end subroutine w_split_adv

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  w_split_it_adv()
!
! !INTERFACE:
   subroutine w_split_it_adv(dt,f,ww,az,splitfac,method)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE , intent(in) :: ww(I3DFIELD),dt
   integer , intent(in)  :: az(E2DFIELD)
   REALTYPE, intent(in)	:: splitfac
   integer, intent(in)	:: method
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)	:: f(I3DFIELD)
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer	:: i,ii,j,jj,k,kk,it
   REALTYPE	:: c,alpha,beta,x,r,Phi,limit,fu,fc,fd,cmax
   LOGICAL      :: READY
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'w_split_it_adv() # ',Ncall
#endif

   cu = _ZERO_

! Calculating w-interface fluxes !

   select case (method)
      case (UPSTREAM_SPLIT)
         do j=jjmin,jjmax
            do i=iimin,iimax
               if (az(i,j) .eq. 1) then
                  do k=1,kmax-1
                     cu(i,j,k) = _ZERO_
                     if (ww(i,j,k) .gt. _ZERO_) then
                        cu(i,j,k)=ww(i,j,k)*f(i,j,k)
                     else
                        cu(i,j,k)=ww(i,j,k)*f(i,j,k+1)
	             end if
                  end do
                  do k=1,kmax   ! Doing a w-advection step
                     hio(i,j,k)=hi(i,j,k)
                     hi(i,j,k)=hio(i,j,k)-splitfac*dt*(ww(i,j,k)-ww(i,j,k-1))
                     f(i,j,k)=(f(i,j,k)*hio(i,j,k)-        &
                               splitfac*dt*(cu(i,j,k)-cu(i,j,k-1)))/hi(i,j,k)
                  end do
	       end if
            end do
         end do
      case ((P2),(Superbee),(MUSCL),(P2_PDM))
         do j=jjmin,jjmax
            do i=iimin,iimax
	       if (az(i,j) .eq. 1) then
	          cmax=0.
		  it=1
		  ready=.false.
111		  do ii=1,it
                  do k=1,kmax-1
                  cu(i,j,k) = _ZERO_
                     if (ww(i,j,k) .gt. _ZERO_) then
                        c=ww(i,j,k)/float(it)*dt/(0.5*(hi(i,j,k)+hi(i,j,k+1)))
			if (c.gt.cmax) cmax=c
		        if (k .gt. 1) then
		           fu=f(i,j,k-1)         ! upstream
			else
		           fu=f(i,j,k)
			end if
		        fc=f(i,j,k  )            ! central
		        fd=f(i,j,k+1)            ! downstream
		        if (abs(fd-fc) .gt. 1e-10) then
		           r=(fc-fu)/(fd-fc)     ! slope ratio
		        else
		           r=(fc-fu)*1.e10
		        end if
		     else
                        c=-ww(i,j,k)/float(it)*dt/(0.5*(hi(i,j,k)+hi(i,j,k+1)))
			if (c.gt.cmax) cmax=c
		        if (k .lt. kmax-1) then
		           fu=f(i,j,k+2)         ! upstream
			else
		           fu=f(i,j,k+1)
			end if
		        fc=f(i,j,k+1)            ! central
		        fd=f(i,j,k  )            ! downstream
		        if (abs(fc-fd) .gt. 1e-10) then
		           r=(fu-fc)/(fc-fd)     ! slope ratio
		        else
		           r=(fu-fc)*1.e10
		        end if
		     end if
                     x = one6th*(1.-2.0*c)
                     Phi=(0.5+x)+(0.5-x)*r
		     select case (method)
		     case ((P2),(P2_PDM))
                        x = one6th*(1.-2.0*c)
                        Phi=(0.5+x)+(0.5-x)*r
			if (method.eq.P2) then
			   limit=Phi
			else
		           limit=max(_ZERO_,min(Phi,2./(1.-c),2.*r/(c+1.e-10)))
			end if
		     case (Superbee)
		        limit=max(_ZERO_, min(1.0, 2.0*r), min(r,2.0) )
		     case (MUSCL) 	
		        limit=max(_ZERO_,min(2.0,2.0*r,0.5*(1.0+r)))
                     case default
                        FATAL 'This is not so good - do_advection_3d()'
	                stop
		     end select
	             cu(i,j,k)=ww(i,j,k)*(fc+0.5*limit*(1-c)*(fd-fc))
                  end do
		  if (.not.READY) then
		     it=min(200,int(cmax)+1)
		     if (it.gt.1) write(95,*) i,j,it,cmax
		     if (it.gt.1) write(*,*) i,j,it,cmax
		  end if
		  if ((it.gt.1).and.(.not.READY)) then
		     READY=.true.
		     goto 111
		  end if
                  do k=1,kmax   ! Doing a w-advection step
                     hio(i,j,k)=hi(i,j,k)
                     hi(i,j,k)=hio(i,j,k)                  &
		      -splitfac/float(it)*dt*(ww(i,j,k)-ww(i,j,k-1))
                     f(i,j,k)=(f(i,j,k)*hio(i,j,k)         &
                      -splitfac/float(it)*dt*(cu(i,j,k)-cu(i,j,k-1)))/hi(i,j,k)
                  end do
                  end do
	       end if
            end do
         end do
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving w_split_it_adv()'
   write(debug,*)
#endif
   return
   end subroutine w_split_it_adv

!-----------------------------------------------------------------------

   end module advection_3d

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------

