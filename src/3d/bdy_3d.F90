!$Id: bdy_3d.F90,v 1.1 2002-05-02 14:00:59 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  bdy_3d - boundary conditions 3D.
!
! !INTERFACE:
   module bdy_3d
!
! !DESCRIPTION:
!  Description still missing
!
! !USES:
   use domain, only: imin,jmin,imax,jmax,H,az,au,av
   use domain, only: iimin,jjmin,iimax,jjmax,kmax
   use domain, only: nsbv,NWB,NNB,NEB,NSB
   use domain, only: wi,wfj,wlj,nj,nfi,nli,ei,efj,elj,sj,sfi,sli
   use variables_3d
   IMPLICIT NONE
!
   private
!
! !PUBLIC DATA MEMBERS:
   public init_bdy_3d, do_bdy_3d
   REALTYPE, public, allocatable	:: S_bdy(:,:),T_bdy(:,:)
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: bdy_3d.F90,v $
!  Revision 1.1  2002-05-02 14:00:59  gotm
!  Initial revision
!
!  Revision 1.1  2001/08/30 08:58:19  bbh
!  Initial import
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_bdy_3d
!
! !INTERFACE:
   subroutine init_bdy_3d()
!
! !DESCRIPTION:
!  Description still missing
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer	:: rc,i,j,k,n
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_bdy_3d() # ',Ncall
#endif

   LEVEL2 'init_bdy_3d()'
   allocate(S_bdy(nsbv,0:kmax),stat=rc)
   if (rc /= 0) stop 'init_init_bdy_3d: Error allocating memory (S_bdy)'

   allocate(T_bdy(nsbv,0:kmax),stat=rc)
   if (rc /= 0) stop 'init_init_bdy_3d: Error allocating memory (T_bdy)'

#ifdef DEBUG
   write(debug,*) 'Leaving init_bdy_3d()'
   write(debug,*)
#endif
   return
   end subroutine init_bdy_3d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_bdy_3d()
!
! !INTERFACE:
   subroutine do_bdy_3d(tag,field)
!
! !DESCRIPTION:
!  Description still missing
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)	:: tag
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)	:: field(I3DFIELD)
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer	:: i,j,k,n,ii,jj
   REALTYPE	:: sp(1:4),rat
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_bdy_3d() # ',Ncall
#endif

#if 0
   select case (tag)
      case (1)
      ! Lateral zero-gradient boundary condition (north & south)
         do k=1,kmax
            do i=iimin,iimax
               if (au(i,jjmin) .eq. 3) field(i,jjmin,k)=field(i,jjmin+1,k)
               if (au(i,jjmax) .eq. 3) field(i,jjmax,k)=field(i,jjmax-1,k)
            end do
         end do
      case (2)
      ! Lateral zero-gradient boundary conditions (west & east)
         do k=1,kmax
            do j=jjmin,jjmax
               if (av(iimin,j) .eq. 3) field(iimin,j,k)=field(iimin+1,j,k)
               if (av(iimax,j) .eq. 3) field(iimax,j,k)=field(iimax-1,j,k)
            end do
         end do
      case default
         FATAL 'Non valid tag'
	 stop 'do_bdy_3d'
   end select
#endif

! Sponge layer factors according to Martinsen and Engedahl, 1987.
   sp(1)=1.0
   sp(2)=0.5625
   sp(3)=0.25
   sp(4)=0.0625

   k = 0
   do n=1,NWB
      i = wi(n)
      do j=wfj(1),wlj(1)
         k = k+1
         do ii=1,4
	    if (az(i-1+ii,j).gt.0) then
               S(i-1+ii,j,:) = sp(ii)*S_bdy(k,:)+(1.-sp(ii))*S(i-1+ii,j,:)
               T(i-1+ii,j,:) = sp(ii)*T_bdy(k,:)+(1.-sp(ii))*T(i-1+ii,j,:)
	    end if
         end do
      end do
   end do
   do n = 1,NNB
      j = nj(n)
      do i = nfi(n),nli(n)
         k = k+1
         do jj=1,4
	    if (az(i,j+1-jj).gt.0) then
               S(i,j+1-jj,:) = sp(jj)*S_bdy(k,:)+(1.-sp(jj))*S(i,j+1-jj,:)
               T(i,j+1-jj,:) = sp(jj)*T_bdy(k,:)+(1.-sp(jj))*T(i,j+1-jj,:)
	    end if
         end do
      end do
   end do
   do n=1,NEB
      i = ei(n)
      do j=efj(1),elj(1)
         k = k+1
         do ii=1,4
	    if (az(i+1-ii,j).gt.0) then
               S(i+1-ii,j,:) = sp(ii)*S_bdy(k,:)+(1.-sp(ii))*S(i+1-ii,j,:)
               T(i+1-ii,j,:) = sp(ii)*T_bdy(k,:)+(1.-sp(ii))*T(i+1-ii,j,:)
	    end if
         end do
      end do
   end do
   do n = 1,NSB
      j = sj(n)
      do i = sfi(n),sli(n)
         k = k+1
         do jj=1,4
	    if (az(i,j-1+jj).gt.0) then
               S(i,j-1+jj,:) = sp(jj)*S_bdy(k,:)+(1.-sp(jj))*S(i,j-1+jj,:)
               T(i,j-1+jj,:) = sp(jj)*T_bdy(k,:)+(1.-sp(jj))*T(i,j-1+jj,:)
	    end if
         end do
      end do
   end do

#ifdef DEBUG
   write(debug,*) 'leaving do_bdy_3d()'
   write(debug,*)
#endif
   return
   end subroutine do_bdy_3d
!EOC

!-----------------------------------------------------------------------

   end module bdy_3d

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Karsten Bolding and Hans Burchard               !
!-----------------------------------------------------------------------
