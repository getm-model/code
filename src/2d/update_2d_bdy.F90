!$Id: update_2d_bdy.F90,v 1.9 2008-12-09 00:31:57 kb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: update_2d_bdy - update 2D boundaries every time step.
!
! !INTERFACE:
   subroutine update_2d_bdy(loop,bdyramp)
!
! !DESCRIPTION:
!
! In this routine sea surface elevation boundary conditions are read
! in from a file, interpolated to the actual time step, and distributed
! to the open boundary grid boxes. 
! Only for a special test case ({\tt SYLT\_TEST}), ascii data reading is
! supported. For a few special simple cases, analytical calculation
! of boundary elevations is supported. The generic way is reading in
! boundary data from a netcdf file, which is managed in
! {\tt get\_2d\_bdy} via {\tt get\_2d\_bdy\_ncdf}.
!
! !USES:
   use domain, only: NWB,NNB,NEB,NSB,H,min_depth,imin,imax,jmin,jmax,az
   use domain, only: wi,wfj,wlj,nj,nfi,nli,ei,efj,elj,sj,sfi,sli
   use domain, only: bdy_index,nsbv
   use domain, only: bdy_2d_type
   use m2d, only: dtm,bdyfmt_2d,bdy_data,bdy_data_u,bdy_data_v
   use variables_2d, only: z,D,U,DU,V,DV
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxc,dyc
#else
   use domain, only: dx,dy
#endif

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: loop,bdyramp
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   logical, save             :: first=.true.
   REALTYPE, save            :: time_array(1000),zbo(1000),zbn(1000)
   REALTYPE, save            :: t,t1,t2
   REALTYPE                  :: a,amp,ratio,fac
   integer                   :: i,j,k,l,n
   REALTYPE, parameter       :: FOUR=4.*_ONE_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'update_2d_bdy() # ',Ncall
#endif

   zbo = _ZERO_
   zbn = _ZERO_

   t = loop*dtm

   select case (bdyfmt_2d)
      case (NO_DATA)
      case (ANALYTICAL)
#define OMEGA 2.*3.141592654/43200.
#ifdef COAST_TEST
         amp = 1.5
         bdy_data = amp*sin(OMEGA*t)
#endif
#ifdef WADDEN_SEA_TEST
         amp = 1.5
         bdy_data = amp*sin(OMEGA*t)
#endif
#ifdef FRESHWATER_LENSE_TEST
         bdy_data=_ZERO_
#endif
#undef OMEGA
#ifdef CHANNEL_TEST
         stop 'CHANNEL_TEST - update2dbdy'
         z(1,:)  =  0.01
         z(51,:) = -0.01
#endif
#undef KBK_TESTING
#ifdef CURVI_TEST
         k = 0
         do n=1,NWB
            i = wi(n)
            if (n .eq. 1) then
               a = 0.05
            else
               a = -0.05
            end if
            do j=wfj(1),wlj(1)
               k = k+1
               bdy_data(k) = a
            end do
         end do
         do n=1,NEB
            i = ei(n)
            do j=efj(1),elj(1)
               k = k+1
               bdy_data(k) = -0.05
            end do
         end do
#endif
      case (ASCII)
#ifdef SYLT_TEST
         if (first) then
            first = .false.

            if (bdyfmt_2d .eq. 1) then
               open(BDYDATA,file='bdy_data.dat')
               t1 = 0.
               do i=1,nsbv
                  read(BDYDATA,*) zbo(i)
               end do
               t2 = 44714./80
               do i=1,nsbv
                  read(BDYDATA,*) zbn(i)
               end do
            end if
         end if

!  Read in new boundary values
         if (t .ge. t2) then
            t1 = t2
            t2 = t2+44714./80
            zbo = zbn

            STDERR 'Reading new boundary values.... '
            do i=1,nsbv
               read(BDYDATA,*) zbn(i)
            end do
         end if
#endif
      case (NETCDF)
!        Read in get_2d_bdy() via get_2d_bdy_ncdf()
      case default
         FATAL 'A non valid communication method has been chosen'
         stop 'update2dbdy'
   end select

!  Data read - do time interpolation

   ratio = _ONE_
   fac = _ONE_
   if(bdyramp .gt. 1) fac=min( _ONE_ ,FOUR*loop/float(bdyramp))

! implicit Sommerfeldt
#if 0
eta = 1/(1+Cx)(eta_old + Cx*eta_new)
               Cx=dtm/dx*srt(g*D)
               z(i,j) = _ONE_/(_ONE_+Cx)*(z(i,j)+Cx*z(i,j-1))
#endif

   l = 0
   do n = 1,NWB
      l = l+1
      k = bdy_index(l)
      i = wi(n)
      do j = wfj(n),wlj(n)
         select case (bdy_2d_type(l))
            case (ZERO_GRADIENT)
               z(i,j) = z(i+1,j)
            case (SOMMERFELDT)
               z(i,j) = z(i,j) + dtm/DXC*sqrt(9.81*D(i,j))*(z(i+1,j)-z(i,j))
            case (CLAMPED)
               z(i,j) = max(fac*bdy_data(k),-H(i,j)+min_depth)
            case (FLATHER_ELEV)
               a = sqrt(DU(i,j)/9.81)*(U(i,j)/DU(i,j)-bdy_data_u(k))
               z(i,j) = max(fac*(bdy_data(k) - a),-H(i,j)+min_depth)
            case default
               FATAL 'Illegal NWB 2D boundary type selection'
               stop 'update_2d_bdy()'
         end select
         k = k+1
      end do
   end do

   do n = 1,NNB
      l = l+1
      k = bdy_index(l)
      j = nj(n)
      do i = nfi(n),nli(n)
         select case (bdy_2d_type(l))
            case (ZERO_GRADIENT)
               z(i,j) = z(i,j-1)
            case (SOMMERFELDT)
               z(i,j) = z(i,j) - dtm/DXC*sqrt(9.81*D(i,j))*(z(i,j)-z(i,j-1))
            case (CLAMPED)
               z(i,j) = max(fac*bdy_data(k),-H(i,j)+min_depth)
            case (FLATHER_ELEV)
               a = sqrt(DV(i,j)/9.81)*(V(i,j-1)/DV(i,j-1)+bdy_data_v(k))
               z(i,j) = max(fac*(bdy_data(k) + a),-H(i,j)+min_depth)
            case default
               FATAL 'Illegal NNB 2D boundary type selection'
               stop 'update_2d_bdy()'
         end select
         k = k+1
      end do
   end do

   do n = 1,NEB
      l = l+1
      k = bdy_index(l)
      i = ei(n)
      do j = efj(n),elj(n)
         select case (bdy_2d_type(l))
            case (ZERO_GRADIENT)
               z(i,j) = z(i-1,j)
            case (SOMMERFELDT)
               z(i,j) = z(i,j) - dtm/DXC*sqrt(9.81*D(i,j))*(z(i,j)-z(i-1,j))
            case (CLAMPED)
               z(i,j) = max(fac*bdy_data(k),-H(i,j)+min_depth)
            case (FLATHER_ELEV)
               a = sqrt(DU(i,j)/9.81)*(U(i-1,j)/DU(i-1,j)-bdy_data_u(k))
               z(i,j) = max(fac*(bdy_data(k) + a),-H(i,j)+min_depth)
            case default
               FATAL 'Illegal NEB 2D boundary type selection'
               stop 'update_2d_bdy()'
         end select
         k = k+1
      end do
   end do

   do n = 1,NSB
      l = l+1
      k = bdy_index(l)
      j = sj(n)
      do i = sfi(n),sli(n)
         select case (bdy_2d_type(l))
            case (ZERO_GRADIENT)
               z(i,j) = z(i,j+1)
            case (SOMMERFELDT)
               z(i,j) = z(i,j) + dtm/DXC*sqrt(9.81*D(i,j))*(z(i,j+1)-z(i,j))
            case (CLAMPED)
               z(i,j) = max(fac*bdy_data(k),-H(i,j)+min_depth)
            case (FLATHER_ELEV)
               a = sqrt(DV(i,j)/9.81)*(V(i,j)/DV(i,j)-bdy_data_v(k))
               z(i,j) = max(fac*(bdy_data(k) - a),-H(i,j)+min_depth)
            case default
               FATAL 'Illegal NSB 2D boundary type selection'
               stop 'update_2d_bdy()'
         end select
         k = k+1
      end do
   end do

#ifdef WADDEN_SEA_TEST
   i=imin
   do j=1,90
      if (az(i+2,j).eq.2) then
         a=z(i+1,j)+z(i+1,j)-z(i+2,j)
      else
         a=z(i+1,j)
      end if
      if (az(i,j).eq.2) z(i,j)= max(a,-H(i,j)+min_depth)
   end do
   i=imax
   do j=1,127
      if (az(i-2,j).eq.2) then
         a=z(i-1,j)+z(i-1,j)-z(i-2,j)
      else
         a=z(i-1,j)
      end if
      if (az(i,j).eq.2) z(i,j)= max(a,-H(i,j)+min_depth)
   end do
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving update_2d_bdy()'
   write(debug,*)
#endif
   return
   end subroutine update_2d_bdy
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
