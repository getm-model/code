!$Id: uv_advect_3d.F90,v 1.2 2003-04-07 13:36:38 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: uv_advect_3d() - mumentum advection - horizontal.
!
! !INTERFACE:
   subroutine uv_advect_3d(hor_adv,ver_adv,strang)
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: iimin,iimax,jjmin,jjmax,kmax,az,au,av,ax
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dyc,arud1,dxx,dyx,arvd1,dxc
#else
   use domain, only: dx,dy,ard1
#endif
   use variables_3d, only: dt,kumin,kvmin,uu,vv,ww,hun,hvn,huo,hvo,uuEx,vvEx
#ifdef UV_TVD
   use variables_3d, only: uadv,vadv,wadv,huadv,hvadv,hoadv,hnadv
#endif
   use advection_3d, only: do_advection_3d
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: hor_adv,ver_adv,strang
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: uv_advect_3d.F90,v $
!  Revision 1.2  2003-04-07 13:36:38  kbk
!  parallel support, cleaned code + NO_3D, NO_BAROCLINIC
!
!  Revision 1.1.1.1  2002/05/02 14:00:57  gotm
!  recovering after CVS crash
!
!  Revision 1.7  2001/10/12 11:39:20  bbh
!  TVD moved out of ??_momentum_3d.F90 and into uv_advect_3d.F90
!
!  Revision 1.6  2001/08/01 08:31:22  bbh
!  CURVILINEAR now implemented
!
!  Revision 1.5  2001/07/26 13:47:18  bbh
!  Fixed some typos
!
!  Revision 1.4  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.3  2001/05/03 20:12:31  bbh
!  Use of variables_3d
!
!  Revision 1.2  2001/05/01 07:13:27  bbh
!  use: kmax from m3d to domain
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
   integer	:: i,j,k,rc
   REALTYPE	:: PP(iimin-1:iimax+1,jjmin-1:jjmax+1,1:kmax)
   REALTYPE	:: www(0:kmax)
#ifdef UV_TVD
   REALTYPE	:: dxuadv(I2DFIELD),dxvadv(I2DFIELD),area_inv(I2DFIELD)
   REALTYPE	:: dyuadv(I2DFIELD),dyvadv(I2DFIELD)
   integer	:: azadv(I2DFIELD),auadv(I2DFIELD),avadv(I2DFIELD)
   REALTYPE	:: AH=_ZERO_
#endif

!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'uv_advect_3d() # ',Ncall
#endif

#ifdef UV_TVD
! Here begins dimensional split advection for u-velocity
   do k=1,kmax
      do j=jjmin,jjmax
         do i=iimin,iimax
            uadv(i,j,k)=0.5*(uu(i+1,j,k)+uu(i,j,k))
            vadv(i,j,k)=0.5*(vv(i+1,j,k)+vv(i,j,k))
            wadv(i,j,k)=0.5*(ww(i+1,j,k)+ww(i,j,k))
            huadv(i,j,k)=0.5*(hun(i+1,j,k)+hun(i,j,k))
            hvadv(i,j,k)=0.5*(hvn(i+1,j,k)+hvn(i,j,k))
            hoadv(i,j,k)=huo(i,j,k)
            hnadv(i,j,k)=hun(i,j,k)
	 end do
      end do
   end do
   do j=jjmin,jjmax
      do i=iimin,iimax
         azadv(i,j)=au(i,j)
         auadv(i,j)=az(i,j)
         avadv(i,j)=ax(i,j)
#if defined(SPHERICAL) || defined(CURVILINEAR) 
         dxuadv(i,j)=dxc(i,j)
         dxvadv(i,j)=dxx(i,j)
         dyuadv(i,j)=dyc(i,j)
         dyvadv(i,j)=dyx(i,j)
	 area_inv(i,j)=arud1(i,j)
#else
         dxuadv(i,j)=dx
         dxvadv(i,j)=dx
         dyuadv(i,j)=dy
         dyvadv(i,j)=dy
	 area_inv(i,j)=1./(dx*dy)
#endif
      end do
   end do

   do k=1,kmax  ! uuEx is here the velocity to be transported.
      do j=jjmin,jjmax
         do i=iimin,iimax
            uuEx(i,j,k)=uu(i,j,k)/huo(i,j,k)
         end do
      end do
   end do
   call do_advection_3d(dt,uuEx,uadv,vadv,wadv,huadv,hvadv,hoadv,hnadv,	&
                        dxuadv,dxvadv,dyuadv,dyvadv,area_inv,		&
			azadv,auadv,avadv,hor_adv,ver_adv,strang,AH)
   uuEx=-(uuEx*hun-uu)/dt ! Here, uuEx is the advection term.  			

! Here begins dimensional split advection for v-velocity
   do k=1,kmax
      do j=jjmin,jjmax
         do i=iimin,iimax
            uadv(i,j,k)=0.5*(uu(i,j+1,k)+uu(i,j,k))
            vadv(i,j,k)=0.5*(vv(i,j+1,k)+vv(i,j,k))
            wadv(i,j,k)=0.5*(ww(i,j+1,k)+ww(i,j,k))
            huadv(i,j,k)=0.5*(hun(i,j+1,k)+hun(i,j,k))
            hvadv(i,j,k)=0.5*(hvn(i,j+1,k)+hvn(i,j,k))
            hoadv(i,j,k)=hvo(i,j,k)
            hnadv(i,j,k)=hvn(i,j,k)
         end do
      end do
   end do
   do j=jjmin,jjmax
      do i=iimin,iimax
         azadv(i,j)=av(i,j)
         auadv(i,j)=ax(i,j)
         avadv(i,j)=az(i,j)
#if defined(SPHERICAL) || defined(CURVILINEAR)
	 dxuadv(i,j)=dxx(i,j)
	 dxvadv(i,j)=dxc(i,j)
	 dyuadv(i,j)=dyx(i,j)
	 dyvadv(i,j)=dyc(i,j)
	 area_inv(i,j)=arvd1(i,j)
#else
	 dxuadv(i,j)=dx
	 dxvadv(i,j)=dx
	 dyuadv(i,j)=dy
	 dyvadv(i,j)=dy
	 area_inv(i,j)=_ONE_/(dx*dy)
#endif
      end do
   end do

   do k=1,kmax   ! vvEx is here the velocity to be transported.
      do j=jjmin,jjmax
         do i=iimin,iimax
            vvEx(i,j,k)=vv(i,j,k)/hvo(i,j,k)
         end do
      end do
   end do
   call do_advection_3d(dt,vvEx,uadv,vadv,wadv,huadv,hvadv,hoadv,hnadv,	&
                        dxuadv,dxvadv,dyuadv,dyvadv,area_inv,		&
                        azadv,auadv,avadv,hor_adv,ver_adv,strang,AH)
   vvEx=-(vvEx*hvn-vv)/dt ! Here, vvEx is the advection term.

#else  ! First-order upstream, one three-dimensional  step

! Upstream for dx(uu^2/hun)
   do k=1,kmax
      do i=iimin,iimax+1          ! PP defined on T-points
         do j=jjmin,jjmax+1
            if (az(i,j) .ge. 1) then
               if (k .ge. kumin(i,j)) then
                  PP(i,j,k)=0.5*(uu(i-1,j,k)+uu(i,j,k))
                  if (PP(i,j,k) .gt. _ZERO_) then
                     PP(i,j,k)=PP(i,j,k)*uu(i-1,j,k)/hun(i-1,j,k)*DYC
                  else
                     PP(i,j,k)=PP(i,j,k)*uu(i,j,k)/hun(i,j,k)*DYC
                  end if
               end if
            end if
         end do
      end do
   end do

   do k=1,kmax
      do j=jjmin,jjmax         ! uuEx defined on U-points
         do i=iimin,iimax
            if (au(i,j).ge.1) then
               if (k .ge. kumin(i,j)) then
                  uuEx(i,j,k)=(PP(i+1,j,k)-PP(i,j,k))*ARUD1
               end if
            end if
         end do
      end do
   end do

! Upstream for dy(uu*vv/hun)
   do k=1,kmax
      do j=jjmin-1,jjmax          ! PP defined on X-points
         do i=iimin-1,iimax
            if (au(i,j) .ge. 1 .or. au(i,j+1) .ge. 1) then
               if (k .ge. kumin(i,j)) then
                  PP(i,j,k)=0.5*(vv(i+1,j,k)+vv(i,j,k))
                  if (PP(i,j,k) .gt. _ZERO_) then
                     PP(i,j,k)=PP(i,j,k)*uu(i,j,k)/hun(i,j,k)*DXX
                  else
                     PP(i,j,k)=PP(i,j,k)*uu(i,j+1,k)/hun(i,j+1,k)*DXX
                  end if
               end if
            end if
         end do
      end do
   end do

   do k=1,kmax
      do j=jjmin,jjmax
         do i=iimin,iimax
            if (au(i,j) .ge. 1) then
               if (k .ge. kumin(i,j)) then
                  uuEx(i,j,k)=(uuEx(i,j,k)+(PP(i,j,k)-PP(i,j-1,k))*ARUD1)
               end if
            end if
         end do
      end do
   end do

! Upstream for dx(uu*vv/hvn)
   do k=1,kmax
      do j=jjmin-1,jjmax         !  PP defined on X-points
         do i=iimin-1,iimax
            if (av(i,j) .ge. 1 .or. av(i+1,j) .ge. 1) then
               if (k .ge. kvmin(i,j)) then
                  PP(i,j,k)=0.5*(uu(i,j,k)+uu(i,j+1,k))
                  if (PP(i,j,k) .gt. _ZERO_) then
                     PP(i,j,k)=PP(i,j,k)*vv(i,j,k)/hvn(i,j,k)*DYX
                  else
                     PP(i,j,k)=PP(i,j,k)*vv(i+1,j,k)/hvn(i+1,j,k)*DYX
                  end if
               end if
            end if
         end do
      end do
   end do

   do k=1,kmax
      do j=jjmin,jjmax          ! vvEx defined on V-points
         do i=iimin,iimax
            if (av(i,j) .ge. 1) then
               if (k .ge. kvmin(i,j)) then
                  vvEx(i,j,k)=(PP(i,j,k)-PP(i-1,j,k))*ARVD1
               end if
            end if
         end do
      end do
   end do

! Upstream for dy(vv^2/hvn)
   do k=1,kmax
      do j=jjmin,jjmax+1
         do i=iimin,iimax+1
            if (az(i,j) .ge. 1) then
               if (k .ge. kvmin(i,j)) then
                  PP(i,j,k)=0.5*(vv(i,j-1,k)+vv(i,j,k))
                  if (PP(i,j,k) .gt. _ZERO_) then
                     PP(i,j,k)=PP(i,j,k)*vv(i,j-1,k)/hvn(i,j-1,k)*DXC
                  else
                     PP(i,j,k)=PP(i,j,k)*vv(i,j,k)/hvn(i,j,k)*DXC
                  end if
               end if
            end if
         end do
      end do
   end do

   do k=1,kmax
      do j=jjmin,jjmax          ! vvEx defined on V-points
         do i=iimin,iimax
            if (av(i,j) .ge. 1) then
               if (k .ge. kvmin(i,j)) then
                  vvEx(i,j,k)=(vvEx(i,j,k)+(PP(i,j+1,k)-PP(i,j,k))*ARVD1)
               end if
            end if
         end do
      end do
   end do

! Upstream for (uu*ww)_k - (uu*ww)_{k-1}
   do j=jjmin,jjmax
      do i=iimin,iimax
         if (au(i,j) .eq. 1) then
            www(kumin(i,j)-1)= _ZERO_
            do k=kumin(i,j),kmax-1
               www(k)=0.5*(ww(i+1,j,k)+ww(i,j,k))
               if (www(k) .gt. _ZERO_) then
	          www(k)=www(k)*uu(i,j,k)/hun(i,j,k)
	       else
	          www(k)=www(k)*uu(i,j,k+1)/hun(i,j,k+1)
	       end if
            end do
            www(kmax)= _ZERO_
            do k=kumin(i,j),kmax
               uuEx(i,j,k)=uuEx(i,j,k)+www(k)-www(k-1)
            end do
	 end if
      end do
   end do

! Upstream for (vv*ww)_k - (vv*ww)_{k-1}
   do j=jjmin+1,jjmax-2
      do i=iimin,iimax
         if (av(i,j) .eq. 1) then
            www(kvmin(i,j)-1)= _ZERO_
            do k=kvmin(i,j),kmax-1
               www(k)=0.5*(ww(i,j+1,k)+ww(i,j,k))
	       if (www(k) .gt. _ZERO_) then
	          www(k)=www(k)*vv(i,j,k)/hvn(i,j,k)
	       else
	          www(k)=www(k)*vv(i,j,k+1)/hvn(i,j,k+1)
               end if
            end do
            www(kmax)= _ZERO_
            do k=kvmin(i,j),kmax
               vvEx(i,j,k)=vvEx(i,j,k)+www(k)-www(k-1)
            end do
	 end if
      end do
   end do

#endif  ! End of three-dimensional first-order upstream

#ifdef DEBUG
   write(debug,*) 'Leaving uv_advect_3d()'
   write(debug,*)
#endif
   return
   end subroutine uv_advect_3d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
