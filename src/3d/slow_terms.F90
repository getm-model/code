!$Id: slow_terms.F90,v 1.1 2002-05-02 14:00:55 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: slow_terms() - .....
!
! !INTERFACE:
   subroutine slow_terms
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: iimin,iimax,jjmin,jjmax,kmax,HU,HV,au,av
   use m2d,    only: Uint,Vint,UEx,VEx,Slru,Slrv,SlUx,SlVx,ru,rv
   use variables_3d,    only: kumin,kvmin,uu,vv,huo,hun,hvo,hvn,idpdx,idpdy
   use variables_3d,    only: ssuo,ssun,ssvo,ssvn,uuEx,vvEx,rru,rrv
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
!  $Log: slow_terms.F90,v $
!  Revision 1.1  2002-05-02 14:00:55  gotm
!  Initial revision
!
!  Revision 1.6  2001/08/27 11:50:17  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.5  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.4  2001/05/20 07:51:40  bbh
!  Internal pressure included
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
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'slow_terms() # ',Ncall
#endif

   if (kmax.gt.1) then

      do i=iimin,iimax
         do j=jjmin,jjmax
	    if (au(i,j).ge.1) then
               SlUx(i,j)=-UEx(i,j)
	    end if
         end do
      end do

      do k=1,kmax
         do j=jjmin,jjmax
            do i=iimin,iimax
	       if (au(i,j).ge.1) then
                  if (k.ge.kumin(i,j)) then
                     SlUx(i,j)=SlUx(i,j)+uuEx(i,j,k)-idpdx(i,j,k)
                  end if
               end if
            end do
         end do
      end do

      do i=iimin,iimax
         do j=jjmin,jjmax
	    if (av(i,j).ge.1) then
               SlVx(i,j)=-VEx(i,j)
	    end if
         end do
      end do

      do k=1,kmax
         do j=jjmin,jjmax
            do i=iimin,iimax
	       if (av(i,j).ge.1) then
                  if (k.ge.kvmin(i,j)) then
                     SlVx(i,j)=SlVx(i,j)+vvEx(i,j,k)-idpdy(i,j,k)
                  end if
	       end if
            end do
         end do
      end do

   else

      do k=1,kmax
         do j=jjmin,jjmax
            do i=iimin,iimax
	       if (au(i,j).ge.1) then
                  if (k.ge.kumin(i,j)) then
                     SlUx(i,j)=-idpdx(i,j,k)
                  end if
               end if
            end do
         end do
      end do

      do k=1,kmax
         do j=jjmin,jjmax
            do i=iimin,iimax
	       if (av(i,j).ge.1) then
                  if (k.ge.kvmin(i,j)) then
                     SlVx(i,j)=-idpdy(i,j,k)
                  end if
	       end if
            end do
         end do
      end do
   endif

   do j=jjmin,jjmax
      do i=iimin,iimax
      if (au(i,j).ge.1) then
         k=kumin(i,j)
	 if (kmax.gt.1) then
#ifdef NO_SLR
            STDERR 'NO_SLR U'
            Slru(i,j)= _ZERO_
#else
            Slru(i,j)=-Uint(i,j)/(0.5*(ssuo(i,j)+ssun(i,j))	&
                                    +HU(i,j))*ru(i,j)		&
                  +uu(i,j,k)/(0.5*(huo(i,j,k)+hun(i,j,k)))*rru(i,j)
#endif
         else
            Slru(i,j)= _ZERO_
         end if
      end if
      end do
   end do

   do j=jjmin,jjmax
      do i=iimin,iimax
         if (av(i,j).ge.1) then
            k=kvmin(i,j)
	    if (kmax.gt.1) then
#ifdef NO_SLR
               STDERR 'NO_SLR V'
               Slrv(i,j)=0.
#else
               Slrv(i,j)=-Vint(i,j)/(0.5*(ssvo(i,j)+ssvn(i,j))	&
                                    +HV(i,j))*rv(i,j)		&
                  +vv(i,j,k)/(0.5*(hvo(i,j,k)+hvn(i,j,k)))*rrv(i,j)
#endif
            else
               Slrv(i,j)=0.
            end if  	
         end if
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'Leaving slow_terms()'
   write(debug,*)
#endif
   return
   end subroutine slow_terms
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
