!$Id: gotm.F90,v 1.8 2004-08-06 15:14:35 hb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gotm() - a wrapper to call GOTM.
!
! !INTERFACE:
   subroutine gotm()
!
! !DESCRIPTION:
!
! !USES:
   use halo_zones, only: update_3d_halo,wait_halo,H_TAG
   use domain, only: iimin,iimax,jjmin,jjmax,kmax,az,min_depth,crit_depth
   use variables_2d, only: D,zub,zvb,z
   use variables_3d, only: dt,kmin,ho,hn,tke,eps,SS,num,taus,taub
#ifndef NO_BAROCLINIC
   use variables_3d, only: NN,nuh
#endif
   use variables_3d, only: avmback,avhback
   use turbulence, only: do_turbulence,cde
   use turbulence, only: tke1d => tke, eps1d => eps, L1d => L
   use turbulence, only: num1d => num, nuh1d => nuh
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: gotm.F90,v $
!  Revision 1.8  2004-08-06 15:14:35  hb
!  num and nuh now properly initialised and no gotm call for CONSTANT_VISCOSITY
!
!  Revision 1.7  2003/12/17 10:22:41  kbk
!  now compiles with -DNO_BAROCLINIC and -DNO_3D
!
!  Revision 1.6  2003/12/16 15:58:54  kbk
!  back ground viscosity and diffusivity (manuel)
!
!  Revision 1.5  2003/05/05 15:51:57  kbk
!  proper update of halo zones for num and nuh
!
!  Revision 1.4  2003/04/23 12:16:34  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.3  2003/04/07 13:36:38  kbk
!  parallel support, cleaned code + NO_3D, NO_BAROCLINIC
!
!  Revision 1.1.1.1  2002/05/02 14:00:54  gotm
!  recovering after CVS crash
!
!  Revision 1.11  2001/10/23 07:06:43  bbh
!  PARABOLIC VISCOSITY --> PARABOLIC_VISCOSITY
!
!  Revision 1.10  2001/10/23 07:05:11  bbh
!  Parabolic viscosity - -DPARABOLIC_VISCOSITY
!
!  Revision 1.8  2001/09/09 20:11:46  bbh
!  Buoyancy production has now right sign
!
!  Revision 1.7  2001/08/27 11:50:17  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.6  2001/07/26 12:54:11  bbh
!  Testing advection schems - using ifdef HAIDVOGEL_TEST
!
!  Revision 1.5  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.4  2001/05/21 13:07:19  bbh
!  dt and cnpar is in variables_3d.F90
!
!  Revision 1.3  2001/05/20 07:49:22  bbh
!  Also diffusivities + partial fix of z0b
!
!  Revision 1.2  2001/05/03 20:12:31  bbh
!  Use of variables_3d
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   REALTYPE                  :: u_taus,u_taub,z0s,z0b
   REALTYPE                  :: h(0:kmax),dry,zz
   REALTYPE                  :: P(0:kmax),B(0:kmax)
   REALTYPE                  :: NN1d(0:kmax),SS1d(0:kmax)
   logical, save             :: first=.true.
   integer, save             :: n = 0
   integer                   :: kk
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'gotm() # ',Ncall
#endif

   do j=jjmin,jjmax
      do i=iimin,iimax

         if (az(i,j) .eq. 1 ) then

            u_taus = sqrt(taus(i,j))
            u_taub = sqrt(taub(i,j))

            h = hn(i,j,:)
            SS1d = SS(i,j,:)
#ifndef NO_BAROCLINIC
            NN1d = NN(i,j,:)
#endif

            P  = num(i,j,:)*SS1d
#ifndef NO_BAROCLINIC
            B  = -nuh(i,j,:)*NN1d
#endif

            tke1d=tke(i,j,:)
            eps1d=eps(i,j,:)
            L1d=cde*tke1d**1.5/eps1d
            num1d=num(i,j,:)
#ifndef NO_BAROCLINIC
            nuh1d=nuh(i,j,:)
#endif

            z0s = 0.1
            z0b=0.5*(max(zub(i-1,j),zub(i,j))+max(zvb(i,j-1),zvb(i,j)))
            if (z0s .gt. D(i,j)/10.) z0s=D(i,j)/10.

#ifdef PARABOLIC_VISCOSITY
            zz = _ZERO_
            do k=1,kmax-1
               zz=zz+hn(i,j,k)
               tke1d(k)=max(1.e-10,3.333333*u_taub**2*(1.-zz/D(i,j)))
               L1d(k)=0.4*(zz+z0b)*sqrt(1.-zz/D(i,j))
               eps1d(k)=0.16431677*tke1d(k)**1.5/L1d(k)
               num1d(k)=0.09*tke1d(k)**2/eps1d(k)
#ifndef NO_BAROCLINIC
               nuh1d(k)=num1d(k)
#endif
            end do
#else
            call do_turbulence(kmax,dt,D(i,j),u_taus,u_taub,z0s,z0b,h, &
                               NN1d,SS1d,P,B)
#endif

            tke(i,j,:) = tke1d
            eps(i,j,:) = eps1d
            num(i,j,:) = num1d + avmback
#ifndef NO_BAROCLINIC
            nuh(i,j,:) = nuh1d + avhback
#endif
         end if
      end do
   end do
   
   call update_3d_halo(num,num,az,iimin,jjmin,iimax,jjmax,kmax,H_TAG)
   call wait_halo(H_TAG)
#ifndef NO_BAROCLINIC
   call update_3d_halo(nuh,nuh,az,iimin,jjmin,iimax,jjmax,kmax,H_TAG)
   call wait_halo(H_TAG)
#endif


#ifdef DEBUG
   write(debug,*) 'Leaving gotm()'
   write(debug,*)
#endif
   return
   end subroutine gotm
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
