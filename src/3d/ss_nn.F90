!$Id: ss_nn.F90,v 1.3 2003-04-23 12:16:34 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ss_nn() - calculates the Prandtl and Brunt-Vaisalla freq.
!
! !INTERFACE:
   subroutine ss_nn()
!
! !DESCRIPTION:
!  Brunt Vaisalla frequency: Calculation has to be done
!  as weighted average over the central point (NNc) and
!  the four neighbouring positions (NNe,NNw,NNn,NNs).
!  Otherwise the turbulence model would create too
!  much noise. Later this stencil should be improved
!  by including isopycnal slopes.
!
! !USES:
   use domain, only: iimin,iimax,jjmin,jjmax,kmax,au,av,az
   use variables_3d, only: kmin,kumin,uu,hun,kvmin,vv,hvn,SS
#ifndef NO_BAROCLINIC
   use variables_3d, only: hn,NN,rho,num
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: ss_nn.F90,v $
!  Revision 1.3  2003-04-23 12:16:34  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 13:36:38  kbk
!  parallel support, cleaned code + NO_3D, NO_BAROCLINIC
!
!  Revision 1.1.1.1  2002/05/02 14:00:55  gotm
!  recovering after CVS crash
!
!  Revision 1.7  2001/10/12 09:22:10  bbh
!  Removed some bugs
!
!  Revision 1.6  2001/09/09 20:12:57  bbh
!  Weighted averaging for shear production coded
!
!  Revision 1.5  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.4  2001/05/22 08:28:27  bbh
!  Fixed calculation of buoyancy frequency - averaging
!
!  Revision 1.3  2001/05/20 07:56:09  bbh
!  Included bouyancy frequency
!
!  Revision 1.2  2001/05/03 20:12:31  bbh
!  Use of variables_3d
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,nb
   REALTYPE                  :: dz,NNc,NNe,NNw,NNn,NNs
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'ss_nn() # ',Ncall
#endif

#undef NEW_SS
#define NEW_SS

!  Prandtl frequency
   do j=jjmin,jjmax
      do i=iimin,iimax
         if (az(i,j) .eq. 1 ) then
            do k=1,kmax-1
! This is an older version which we should keep here.
#ifndef NEW_SS
              SS(i,j,k)=0.5* (                                               &
                   ( (uu(i,j,k+1)/hun(i,j,k+1)-uu(i,j,k)/hun(i,j,k))         &
                   /(0.5*(hun(i,j,k+1)+hun(i,j,k))) )**2                     &
                +  ( (uu(i-1,j,k+1)/hun(i-1,j,k+1)-uu(i-1,j,k)/hun(i-1,j,k)) &
                   /(0.5*(hun(i-1,j,k+1)+hun(i-1,j,k))) )**2                 &
                +  ( (vv(i,j,k+1)/hvn(i,j,k+1)-vv(i,j,k)/hvn(i,j,k))         &
                   /(0.5*(hvn(i,j,k+1)+hvn(i,j,k))) )**2                     &
                +  ( (vv(i,j-1,k+1)/hvn(i,j-1,k+1)-vv(i,j-1,k)/hvn(i,j-1,k)) &
                   /(0.5*(hvn(i,j-1,k+1)+hvn(i,j-1,k))) )**2                 &
                            )
#else
! This version should better conserve energy.
              SS(i,j,k)=0.5* (                                                &
                   (uu(i,j,k+1)/hun(i,j,k+1)-uu(i,j,k)/hun(i,j,k))**2         &
                   /(0.5*(hun(i,j,k+1)+hun(i,j,k)))                           &
                    *0.5*(num(i,j,k)+num(i+1,j,k))                            &
               +  (uu(i-1,j,k+1)/hun(i-1,j,k+1)-uu(i-1,j,k)/hun(i-1,j,k))**2 &
                  /(0.5*(hun(i-1,j,k+1)+hun(i-1,j,k)))                       &
                   *0.5*(num(i-1,j,k)+num(i,j,k))                            &
                +  (vv(i,j,k+1)/hvn(i,j,k+1)-vv(i,j,k)/hvn(i,j,k))**2         &
                   /(0.5*(hvn(i,j,k+1)+hvn(i,j,k)))                           &
                    *0.5*(num(i,j,k)+num(i,j+1,k))                            &
                +  (vv(i,j-1,k+1)/hvn(i,j-1,k+1)-vv(i,j-1,k)/hvn(i,j-1,k))**2 &
                   /(0.5*(hvn(i,j-1,k+1)+hvn(i,j-1,k)))                       &
                    *0.5*(num(i,j-1,k)+num(i,j,k))                            &
                            )/(0.5*(hn(i,j,k)+hn(i,j,k+1)))/num(i,j,k)
#endif
            end do
         end if
      end do
   end do

#ifndef NO_BAROCLINIC
#define NEW_NN
#undef NEW_NN
   do j=jjmin,jjmax
      do i=iimin,iimax
         if (az(i,j) .eq. 1 ) then
            do k=kmax-1,1,-1
               dz=0.5*(hn(i,j,k+1)+hn(i,j,k))
               NNc =(rho(i,j,k+1)-rho(i,j,k))/dz
#ifndef NEW_NN
               if (az(i+1,j) .eq. 1) then
                  dz=0.5*(hn(i+1,j,k+1)+hn(i+1,j,k))
                  NNe=(rho(i+1,j,k+1)-rho(i+1,j,k))/dz
               else
                  NNe=NNc
               end if 
               if (az(i-1,j) .eq. 1) then
                  dz=0.5*(hn(i-1,j,k+1)+hn(i-1,j,k))
                  NNw=(rho(i-1,j,k+1)-rho(i-1,j,k))/dz
               else
                  NNw=NNc
               end if
               if (az(i,j+1) .eq. 1) then
                  dz=0.5*(hn(i,j+1,k+1)+hn(i,j+1,k))
                  NNn=(rho(i,j+1,k+1)-rho(i,j+1,k))/dz
               else
                  NNn=NNc
               end if
               if (az(i,j-1) .eq. 1) then
                  dz=0.5*(hn(i,j-1,k+1)+hn(i,j-1,k))
                  NNs=(rho(i,j-1,k+1)-rho(i,j-1,k))/dz
               else
                  NNs=NNc
               end if
               NN(i,j,k)=0.3333333*NNc+0.1666666*(NNe+NNw+NNn+NNs)
#else
               NN(i,j,k)=NNc
#endif
            end do
         end if
      end do
   end do
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving ss_nn()'
   write(debug,*)
#endif
   return
   end subroutine ss_nn
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
