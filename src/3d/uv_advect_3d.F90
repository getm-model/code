!$Id: uv_advect_3d.F90,v 1.11 2006-02-10 22:41:56 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: uv_advect_3d - 3D momentum advection \label{sec-uv-advect-3d}
!
! !INTERFACE:
   subroutine uv_advect_3d(hor_adv,ver_adv,adv_split)
!
! !DESCRIPTION:
!
! For the discretisation of the momentum advection terms, two 
! conceptionally different methods have been implemented in GETM. 
! The first is the straight-forward multidimensional advection
! scheme, which is here realised as the first-order upwind scheme,
! see paragraph {\bf Multidimensional approach} on page 
! \pageref{uvadvect-multi}. 
! 
! In order to make use of the higher-order directional-split
! methods for tracers (see section \ref{sec-do-advection-3d}),
! an alternative method is implemented, in which the complete advection step
! is first made, and then the resulting advection terms,
! which are needed for the calculation of the slow terms, see equations
! (\ref{SxA}) and (\ref{SyA})) are calculated from this
! (see paragraph {\bf Directional-split approach} on page
! \pageref{uvadvect-direct}). 
!
! The choice which of the two methods to be used is made by means
! of the compiler option {\tt UV\_TVD} which has to be set in the
! {\tt Makefile} of the application in order to activate the
! more accurate but computationally more demanding high-order directional-split
! method. The effect of hih-order advaction can be impressively studied
! by means of the freshwater lense test case described in detail by
! \cite{BURCHARDea02}.
!
! When working with the option {\tt SLICE\_MODEL}, the calculation of
! all gradients in $y$-direction is suppressed.
! 
! \paragraph{Multidimensional approach}\label{uvadvect-multi} 
! 
! The advective terms in the momentum equation are discretised in
! a momentum-conservative form. This is carried out here for the 
! advective terms in the $u$-equation (\ref{uEqviCurvi}) and the 
! $v$-equation (\ref{vEqviCurvi}) (after multiplying these
! equations with $mn$). 
! 
! First advection term in (\ref{uEqviCurvi}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! \left(mn\,\partial_{\cal X}\left(\frac{u_kp_k}{n}\right)\right)_{i,j,k}\approx \\ \\
! \quad
! \displaystyle
! \frac{
! \frac12(p_{i+1,j,k}+p_{i,j,k})\tilde u_{i+1,j,k}\Delta y^c_{i+1,j}-
! \frac12(p_{i,j,k}+p_{i-1,j,k})\tilde u_{i,j,k}\Delta y^c_{i,j}
! }{\Delta x^u_{i,j}\Delta y^u_{i,j}}
! \end{array}
! \end{equation}
! 
! For an upwind scheme, the inter-facial velocities which are defined
! on T-points are here
! calculated as:
! 
! \begin{equation}
! \tilde u_{i,j,k}=
! \left\{
! \begin{array}{ll}
! u_{i-1,j,k} & \mbox{ for } \frac12(p_{i,j,k}+p_{i-1,j,k})>0\\ \\
! u_{i,j,k} & \mbox{ else. } 
! \end{array}
! \right.
! \end{equation}
! 
! Second advection term in (\ref{uEqviCurvi}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle 
! \left(mn\,\partial_{\cal Y}y\left(\frac{v_kp_k}{m}\right)\right)_{i,j,k}\approx \\ \\ 
! \displaystyle 
! \quad
! \frac{
! \frac12(q_{i+1,j,k}+q_{i,j,k})\tilde u_{i,j,k}\Delta x^+_{i,j}-
! \frac12(q_{i+1,j-1,k}+q_{i,j-1,k})\tilde u_{i,j-1,k}\Delta x^+_{i,j-1}
! }{\Delta x^u_{i,j}\Delta y^u_{i,j}}
! \end{array}
! \end{equation}
! 
! For an upwind scheme, the inter-facial velocities which are defined on
! X-points are here
! calculated as:
! 
! \begin{equation}
! \tilde u_{i,j,k}=
! \left\{
! \begin{array}{ll}
! u_{i,j,k} & \mbox{ for } \frac12(q_{i+1,j,k}+q_{i,j,k})>0\\ \\
! u_{i,j+1,k} & \mbox{ else. } 
! \end{array}
! \right.
! \end{equation}
! 
! First advection term in (\ref{vEqviCurvi}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle 
! \left(mn\,\partial_{\cal X}\left(\frac{v_kq_k}{n}\right)\right)_{i,j,k}\approx \\ \\ 
! \displaystyle 
! \quad
! \frac{
! \frac12(p_{i,j+1,k}+p_{i,j,k})\tilde v_{i,j,k}\Delta y^+_{i,j}-
! \frac12(p_{i-1,j+1,k}+p_{i-1,j,k})\tilde v_{i-1,j,k}\Delta y^+_{i-1,j}
! }{\Delta x^v_{i,j}\Delta y^v_{i,j}}
! \end{array}
! \end{equation}
! 
! For an upwind scheme, the interfacial velocities which are defined on
! X-points are here
! calculated as:
! 
! \begin{equation}
! \tilde v_{i,j,k}=
! \left\{
! \begin{array}{ll}
! v_{i,j,k} & \mbox{ for } \frac12(p_{i+1,j,k}+p_{i,j,k})>0\\ \\
! v_{i+1,j,k} & \mbox{ else. } 
! \end{array}
! \right.
! \end{equation}
! 
! Second advection term in (\ref{vEqviCurvi}):
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! \left(mn\,\partial_{\cal Y}\left(\frac{v_kq_k}{m}\right)\right)_{i,j,k}\approx \\ \\
! \quad
! \displaystyle
! \frac{
! \frac12(q_{i,j+1,k}+q_{i,j,k})\tilde v_{i,j+1,k}\Delta x^c_{i,j+1}-
! \frac12(q_{i,j,k}+q_{i,j-1,k})\tilde v_{i,j,k}\Delta x^c_{i,j}
! }{\Delta x^v_{i,j}\Delta y^v_{i,j}}
! \end{array}
! \end{equation}
! 
! For an upwind scheme, the interfacial velocities which are defined
! on T-points are here
! calculated as:
! 
! \begin{equation}
! \tilde v_{i,j,k}=
! \left\{
! \begin{array}{ll}
! v_{i,j-1,k} & \mbox{ for } \frac12(q_{i,j,k}+q_{i,j-1,k})>0\\ \\
! v_{i,j,k} & \mbox{ else. } 
! \end{array}
! \right.
! \end{equation}
! 
! The vertical advection terms in equations (\ref{uEqviCurvi})
! and (\ref{vEqviCurvi})
! can be discretised in an upstream scheme as well. 
! 
! Vertical advective flux in (\ref{uEqviCurvi}): 
! \begin{equation}
! \left(\bar w_{k} \tilde u_{k}\right)_{i,j}\approx
! \frac12(w_{i+1,j,k}+w_{i,j,k})\tilde u_{i,j,k} 
! \end{equation}
! 
! with 
! 
! \begin{equation}
! \tilde u_{i,j,k}=
! \left\{
! \begin{array}{ll}
! u_{i,j,k} & \mbox{ for } \frac12(w_{i+1,j,k}+w_{i,j,k}) > 0,\\ \\
! u_{i,j,k+1} & \mbox{ else.}
! \end{array}
! \right.
! \end{equation}
! 
! Vertical advective flux in (\ref{vEqviCurvi}): 
! \begin{equation}
! \left(\bar w_{k} \tilde v_{k}\right)_{i,j}\approx
! \frac12(w_{i,j+1,k}+w_{i,j,k})\tilde v_{i,j,k} 
! \end{equation}
! 
! with 
! 
! \begin{equation}
! \tilde v_{i,j,k}=
! \left\{
! \begin{array}{ll}
! v_{i,j,k} & \mbox{ for } \frac12(w_{i,j+1,k}+w_{i,j,k}) > 0,\\ \\
! v_{i,j,k+1} & \mbox{ else.}
! \end{array}
! \right.
! \end{equation}
! 
! 
! \paragraph{Directional-split approach}\label{uvadvect-direct} 
! 
! Multidimensional treatment of advective terms in three-dimensional
! models is often quite unhandy, especially when higher-order
! advection schemes are needed. On the other hand, directional-split
! methods (which update the advected fields
! in each directional step and then "forget" the advection terms)
! as discussed in section \ref{sec-do-advection-3d} on page 
! \pageref{sec-do-advection-3d}, 
! cannot directly be used 
! for momentum advection when the models are based on mode splitting as
! e.g.\ GETM. The reason for this is that the three-dimensional 
! advection terms are also needed for calculating the slow terms of the
! barotropic (external) mode, see equations (\ref{SxA}) and (\ref{SyA}). 
! 
! The procedure suggested here is as follows. First, the pure momentum advection
! equations are formally solved with the directional-split method
! described in section \ref{sec-do-advection-3d}:
! 
! \begin{equation}\label{uEqvi_adv}
! \begin{array}{l}
! \partial_t p_k
! +\partial_x(u_kp_k)+\partial_y(v_kp_k)
! +\bar w_k \tilde u_k -\bar w_{k-1} \tilde u_{k-1}
! = 0,
! \end{array}
! \end{equation}
! 
! \begin{equation}\label{vEqvi_adv}
! \partial_t q_k
! +\partial_x(u_kq_k)+\partial_y(v_kq_k)
! +\bar w_k  \tilde v_k -\bar w_{k-1}  \tilde v_{k-1}
! =0.
! \end{equation}
! 
! The new solutions $\hat p_{i,j,k}$ and $\hat q_{i,j,k}$ 
! are however not further used, but instead the resulting advective terms
! $-(\hat p_{i,j,k}-p_{i,j,k})/\Delta t$ and 
! $-(\hat q_{i,j,k}-q_{i,j,k})/\Delta t$
! are later applied to the momentum equations (together with other
! processes such as horizontal diffusion, pressure gradients, etc.)
! and also used for the calculation of the slow terms in
! (\ref{SxA}) and (\ref{SyA}). 
! 
! With this method, all higher-order directional-split advection schemes
! are now available for the momentum advection. The advective 
! fluxes needed for this have to be averaged from the conservative
! advective fluxes resulting from the continuity equation
! Continuity will
! still be retained due to the linearity of the continuity equation.    
! 
! !
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
   use halo_zones, only: update_3d_halo,wait_halo,U_TAG,V_TAG

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: hor_adv,ver_adv,adv_split
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: uv_advect_3d.F90,v $
!  Revision 1.11  2006-02-10 22:41:56  hb
!  Source code documentation extended
!
!  Revision 1.10  2006-02-04 11:47:26  hb
!  Source code documentation extended
!
!  Revision 1.9  2005-10-06 09:54:01  hb
!  added support for vertical slice model - via -DSLICE_MODEL
!
!  Revision 1.8  2005/05/25 10:32:13  kbk
!  merged from stabe branch v1_2_1
!
!  Revision 1.7.2.1  2005/05/25 08:41:38  kbk
!  fixed loop boundaries + update HALO's when -DUV_TVD
!
!  Revision 1.7  2003/08/28 15:20:37  kbk
!  use ax mask, always set PP
!
!  Revision 1.6  2003/08/14 13:00:40  kbk
!  do not use masks calculating adv. velocities
!
!  Revision 1.5  2003/06/28 10:40:41  kbk
!  changed loop order
!
!  Revision 1.4  2003/05/02 06:55:49  hb
!  momemtum advection only for mask=1
!
!  Revision 1.3  2003/04/23 12:16:34  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 13:36:38  kbk
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
   integer                   :: i,j,k,rc
   REALTYPE                  :: PP(iimin-1:iimax+1,jjmin-1:jjmax+1,1:kmax)
   REALTYPE                  :: www(0:kmax)
#ifdef UV_TVD
   integer                   :: azadv(I2DFIELD),auadv(I2DFIELD),avadv(I2DFIELD)
   REALTYPE                  :: dxuadv(I2DFIELD),dxvadv(I2DFIELD),area_inv(I2DFIELD)
   REALTYPE                  :: dyuadv(I2DFIELD),dyvadv(I2DFIELD)
   REALTYPE                  :: AH=_ZERO_
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
      do j=jjmin-HALO,jjmax+HALO
         do i=iimin-HALO,iimax+HALO-1
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

   do j=jjmin-HALO,jjmax+HALO
      do i=iimin-HALO,iimax+HALO
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

   call update_3d_halo(uuEx,uuEx,au,iimin,jjmin,iimax,jjmax,kmax,U_TAG)
   call wait_halo(U_TAG)

   call do_advection_3d(dt,uuEx,uadv,vadv,wadv,huadv,hvadv,hoadv,hnadv,&
                        dxuadv,dxvadv,dyuadv,dyvadv,area_inv,          &
                        azadv,auadv,avadv,hor_adv,ver_adv,adv_split,AH)

   uuEx=-(uuEx*hun-uu)/dt ! Here, uuEx is the advection term.

! Here begins dimensional split advection for v-velocity
   do k=1,kmax
      do j=jjmin-HALO,jjmax+HALO-1
         do i=iimin-HALO,iimax+HALO
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

   do j=jjmin-HALO,jjmax+HALO
      do i=iimin-HALO,iimax+HALO
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

   call update_3d_halo(vvEx,vvEx,av,iimin,jjmin,iimax,jjmax,kmax,V_TAG)
   call wait_halo(V_TAG)

   call do_advection_3d(dt,vvEx,uadv,vadv,wadv,huadv,hvadv,hoadv,hnadv,&
                        dxuadv,dxvadv,dyuadv,dyvadv,area_inv,          &
                        azadv,auadv,avadv,hor_adv,ver_adv,adv_split,AH)

   vvEx=-(vvEx*hvn-vv)/dt ! Here, vvEx is the advection term.

#else  ! First-order upstream, one three-dimensional  step

! Upstream for dx(uu^2/hun)
   do k=1,kmax
      do j=jjmin,jjmax             ! PP defined on T-points
         do i=iimin,iimax+1
            PP(i,j,k)=_ZERO_
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
            if (au(i,j) .eq. 1) then
               if (k .ge. kumin(i,j)) then
                  uuEx(i,j,k)=(PP(i+1,j,k)-PP(i,j,k))*ARUD1
               end if
            end if
         end do
      end do
   end do

#ifndef SLICE_MODEL
! Upstream for dy(uu*vv/hun)
   do k=1,kmax
      do j=jjmin-1,jjmax          ! PP defined on X-points
         do i=iimin,iimax
            PP(i,j,k)=_ZERO_
            if (ax(i,j) .ge. 1) then
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
            if (au(i,j) .eq. 1) then
               if (k .ge. kumin(i,j)) then
                  uuEx(i,j,k)=(uuEx(i,j,k)+(PP(i,j,k)-PP(i,j-1,k))*ARUD1)
               end if
            end if
         end do
      end do
   end do
#endif

! Upstream for dx(uu*vv/hvn)
   do k=1,kmax
      do j=jjmin-1,jjmax         !  PP defined on X-points
         do i=iimin-1,iimax
            PP(i,j,k)=_ZERO_
            if (ax(i,j) .ge. 1) then
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
            if (av(i,j) .eq. 1) then
               if (k .ge. kvmin(i,j)) then
                  vvEx(i,j,k)=(PP(i,j,k)-PP(i-1,j,k))*ARVD1
               end if
            end if
         end do
      end do
   end do

#ifndef SLICE_MODEL
! Upstream for dy(vv^2/hvn)
   do k=1,kmax
      do j=jjmin,jjmax+1
         do i=iimin,iimax
            PP(i,j,k)=_ZERO_
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
            if (av(i,j) .eq. 1) then
               if (k .ge. kvmin(i,j)) then
                  vvEx(i,j,k)=(vvEx(i,j,k)+(PP(i,j+1,k)-PP(i,j,k))*ARVD1)
               end if
            end if
         end do
      end do
   end do
#endif

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
   do j=jjmin,jjmax
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
