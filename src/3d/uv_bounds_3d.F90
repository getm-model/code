!$Id: uv_bounds_3d.F90,v 1.1.1.1 2002-05-02 14:01:00 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: uv_bounds_3d() - boundaries for 3D u and v.
!
! !INTERFACE:
   subroutine uv_bounds_3d
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: iimin,iimax,jjmin,jjmax,kmax,HU,HV
   use m2d,    only: Uint,Vint
   use variables_3d,    only: kumin,kvmin,uu,vv,hun,hvn,ssun,ssvn
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
!  $Log: uv_bounds_3d.F90,v $
!  Revision 1.1.1.1  2002-05-02 14:01:00  gotm
!  recovering after CVS crash
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
!
! !LOCAL VARIABLES:
   integer	:: i,j,k
   REALTYPE 	:: x
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'uv_bounds_3d() # ',Ncall
#endif

   do j=jjmin,jjmax

      i=iimin-1              ! uu at western boundary
      x = Uint(i,j)/(ssun(i,j)+HU(i,j))
      do k=kumin(i,j),kmax
!kbk         uu(i,j,k)=Uint(i,j)/(ssun(i,j)+HU(i,j))*hun(i,j,k)
         uu(i,j,k)=x*hun(i,j,k)
      end do

      i=iimax                ! uu at eastern boundary
      x = Uint(i,j)/(ssun(i,j)+HU(i,j))
      do k=kumin(i,j),kmax
!kbk         uu(i,j,k)=Uint(i,j)/(ssun(i,j)+HU(i,j))*hun(i,j,k)
         uu(i,j,k)=x*hun(i,j,k)
      end do

   end do

   do i=iimin,iimax

      j=jjmin-1              ! vv at southern boundary
      x = Vint(i,j)/(ssvn(i,j)+HV(i,j))
      do k=kvmin(i,j),kmax
!kbk         vv(i,j,k)=Vint(i,j)/(ssvn(i,j)+HV(i,j))*hvn(i,j,k)
         vv(i,j,k)=x*hvn(i,j,k)
      end do

      j=jjmax                ! vv at northern boundary
      x = Vint(i,j)/(ssvn(i,j)+HV(i,j))
      do k=kvmin(i,j),kmax
!kbk         vv(i,j,k)=Vint(i,j)/(ssvn(i,j)+HV(i,j))*hvn(i,j,k)
         vv(i,j,k)=x*hvn(i,j,k)
      end do

   end do

#ifdef DEBUG
   write(debug,*) 'Leaving uv_bounds_3d()'
   write(debug,*)
#endif
   return
   end subroutine uv_bounds_3d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
