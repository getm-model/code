!$Id: bottom_friction.F90,v 1.2 2003-03-20 15:48:12 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: bottom_friction() - calculates the 2D bottom friction.
!
! !INTERFACE:
   subroutine bottom_friction(runtype)
!
! !DESCRIPTION:
!
! !USES:
   use parameters, only: kappa,avmmol
   use domain,     only: imin,imax,jmin,jmax,au,av,min_depth
   use variables_2d
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)	:: runtype
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: bottom_friction.F90,v $
!  Revision 1.2  2003-03-20 15:48:12  gotm
!  fixed small bug in calc. of rvv + cleaning the code
!
!  Revision 1.1.1.1  2002/05/02 14:00:41  gotm
!  recovering after CVS crash
!
!  Revision 1.6  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.5  2001/05/20 07:46:38  bbh
!  Fixed masks
!
!  Revision 1.4  2001/05/18 13:03:34  bbh
!  Optimize for speed + masks for update_2d_halo() - CHECK
!
!  Revision 1.3  2001/05/03 19:35:01  bbh
!  Use of variables_2d
!
!  Revision 1.2  2001/04/24 08:24:58  bbh
!  Use runtype instead of macro
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
!  !LOCAL VARIABLES:
   integer 	:: i,j
   REALTYPE 	:: uloc(E2DFIELD),vloc(E2DFIELD)
   REALTYPE 	:: HH(E2DFIELD),fricvel(E2DFIELD)
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'bottom_friction() # ',Ncall
#endif

   do j=jmin,jmax
      do i=imin,imax
         if (au(i,j) .gt. 0) then
            vloc(i,j)=0.25* ( V(i  ,j  )/DV(i  ,j  )	&
                             +V(i+1,j  )/DV(i+1,j  )	&
                             +V(i  ,j-1)/DV(i  ,j-1)	&
                             +V(i+1,j-1)/DV(i+1,j-1) )
         else
            vloc(i,j) = _ZERO_
         end if
      end do
   end do

!  The x-direction

#ifndef DEBUG
   where (au .gt. 0)
      uloc=U/DU
      HH=max(min_depth,DU)
      ruu=(kappa/log((zub+0.5*HH)/zub))**2
   end where
#else
   uloc=U/DU
   HH=max(min_depth,DU)
   ruu=(zub+0.5*HH)/zub

   do j=jmin,jmax
      do i=imin,imax
         if (ruu(i,j) .le. _ONE_) then
            STDERR 'bottom_friction friction coefficient infinite.'
            STDERR 'i = ',i,' j = ',j
            STDERR 'min_depth = ',min_depth,' DU = ',DU(i,j)
            stop
         end if
      end do
   end do
   where (au .gt. 0)
      ruu=(kappa/log(ruu))**2
   end where
#endif

   if (runtype .eq. 1) then
      where (au .gt. 0)
         fricvel=sqrt(ruu*(uloc**2+vloc**2))
         zub=min(HH,zub0+0.1*avmmol/max(avmmol,fricvel))
         ruu=(zub+0.5*HH)/zub
         ruu=(kappa/log(ruu))**2
      end where
   end if

   where (au .gt. 0)
      ru=ruu*sqrt(uloc**2+vloc**2)
   end where

!  The y-direction
   do j=jmin,jmax
      do i=imin,imax
         if (av(i,j) .gt. 0) then
            uloc(i,j)=0.25* ( U(i  ,j  )/DU(i  ,j  )	&
                             +U(i-1,j  )/DU(i-1,j  )	&
                             +U(i  ,j+1)/DU(i  ,j+1)	&
                             +U(i-1,j+1)/DU(i-1,j+1) )
         else
            uloc(i,j) = _ZERO_
         end if
      end do
   end do

#ifndef DEBUG
   where (av .gt. 0)
      vloc=V/DV
      HH=max(min_depth,DV)
      rvv=(kappa/log((zvb+0.5*HH)/zvb))**2
   end where
#else
   vloc=V/DV
   HH=max(min_depth,DV)
   rvv=(zvb+0.5*HH)/zvb

   do j=jmin,jmax
      do i=imin,imax
         if (rvv(i,j) .le. _ONE_) then
            STDERR 'bottom_friction friction coefficient infinite.'
            STDERR 'i = ',i,' j = ',j
            STDERR 'min_depth = ',min_depth,' DV = ',DU(i,j)
            stop
         end if
      end do
   end do

   where (av .gt. 0)
      rvv=(kappa/log(rvv))**2
   end where
#endif

   if (runtype .eq. 1) then
      where (av .gt. 0)
         fricvel=sqrt(rvv*(uloc**2+vloc**2))
         zvb=min(HH,zvb0+0.1*avmmol/max(avmmol,fricvel))
         rvv=(zvb+0.5*HH)/zvb
         rvv=(kappa/log(rvv))**2
      end where
   end if

   where (av .gt. 0)
      rv=rvv*sqrt(uloc**2+vloc**2)
   end where

#ifdef DEBUG
   write(debug,*) 'Leaving bottom_friction()'
   write(debug,*)
#endif
   return
   end subroutine bottom_friction
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
