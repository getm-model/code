!$Id: halo_zones.F90,v 1.1 2002-05-02 14:01:29 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: halo_zones - update halo zones in 'getm'
!
! !INTERFACE:
   module halo_zones

! !DESCRIPTION:
!  This module is included only to supply myid and nprocs used in various
!  places in 'getm'. From version 1.4 real use of MPI will be implemented.
!
! !USES:
   use halo_mpi
   IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_halo_zones,update_2d_halo,update_3d_halo,wait_halo
!
! !PUBLIC DATA MEMBERS:
   integer, parameter	:: H_TAG=10,HU_TAG=11,HV_TAG=12
   integer, parameter	:: D_TAG=20,DU_TAG=21,DV_TAG=22
   integer, parameter   :: z_TAG=30,U_TAG=31,V_TAG=32
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: halo_zones.F90,v $
!  Revision 1.1  2002-05-02 14:01:29  gotm
!  Initial revision
!
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_halo_zones - initialize MPI environment
!
! !INTERFACE:
   subroutine init_halo_zones(parallel,iextr,jextr,imin,imax,jmin,jmax,kmax)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Initialize Parallel environment
!
! !INPUT PARAMETERS:
   logical, intent(in)	:: parallel
   integer, intent(in), optional	:: iextr,jextr,kmax
#ifdef STATIC
   integer, intent(in), optional	:: imin,imax,jmin,jmax
#endif
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
#ifndef STATIC
   integer, intent(out), optional	:: imin,imax,jmin,jmax
#endif
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_halo_zones'
#endif

   if ( .not. parallel ) then
      myid   = -1
      nprocs =  1
   else
      call init_mpi(iextr,jextr,imin,imax,jmin,jmax,kmax)
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving init_halo_zones()'
   write(debug,*)
#endif
   return
   end subroutine init_halo_zones
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: update_2d_halo - updates the halo zones for 2D fields.
!
! !INTERFACE:
   subroutine update_2d_halo(f1,f2,mask,imin,jmin,imax,jmax,tag)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Print information on the MPI environment
!
! !INPUT PARAMETERS:
   integer, intent(in)		:: imin,jmin,imax,jmax
   integer, intent(in)		:: tag
   integer, intent(in) 		:: mask(-HALO+1:,-HALO+1:)
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout) 	:: f1(E2DFIELD),f2(E2DFIELD)
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer	:: i,j,k
   integer	:: ilow,jlow,ihigh,jhigh
!EOP
!-------------------------------------------------------------------------
!BOC
   select case (tag)
      case(HU_TAG, U_TAG , DU_TAG) ! for variables defined on u-grid
         ilow=imin;ihigh=imax-1;jlow=jmin;jhigh=jmax
      case(HV_TAG, V_TAG , DV_TAG) ! for variables defined on v-grid
         ilow=imin;ihigh=imax;jlow=jmin;jhigh=jmax-1
      case default                 ! for variables defined on scalar-grid
         ilow=imin;ihigh=imax;jlow=jmin;jhigh=jmax
   end select
   if (nprocs .eq. 1) then
      STDERR 'one process'
#if 1
      forall(i=ilow:ilow,j=jlow:jhigh, mask(i,j) .ne. 0) &
            f1(i-1,j) = f2(i,j)
      forall(i=ihigh:ihigh,j=jlow:jhigh, mask(i,j) .ne. 0) &
            f1(i+1,j) = f2(i,j)
      forall(i=ilow:ihigh,j=jlow:jlow, mask(i,j) .ne. 0) &
            f1(i,j-1) = f2(i,j)
      forall(i=ilow:ihigh,j=jhigh:jhigh, mask(i,j) .ne. 0) &
            f1(i,j+1) = f2(i,j)
#else
      f1(ilow -1, : )  = f2(ilow,  :  )
      f1(ihigh+1, : )  = f2(ihigh, :  )
      f1( :, jlow -1 ) = f2( :, jlow  )
      f1( :, jhigh+1 ) = f2( :, jhigh )
#endif
   else
      call update_2d_halo_mpi(f1,f2,imin,jmin,imax,jmax,tag)
   end if
   return
   end subroutine update_2d_halo
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: update_3d_halo - updates the halo zones for 3D fields.
!
! !INTERFACE:
   subroutine update_3d_halo(f1,f2,mask,iimin,jjmin,iimax,jjmax,kmax,tag)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Print information on the MPI environment
!
! !INPUT PARAMETERS:
   integer, intent(in)	:: iimin,jjmin,iimax,jjmax,kmax
   integer, intent(in)	:: tag
   integer, intent(in)	:: mask(-HALO+1:,-HALO+1:)
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout):: f1(I3DFIELD),f2(I3DFIELD)
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer	:: i,j,k
   integer	:: ilow,jlow,ihigh,jhigh
!EOP
!-------------------------------------------------------------------------
!BOC
   select case (tag)
      case(HU_TAG, U_TAG , DU_TAG) ! for variables defined on u-grid
         ilow=iimin;ihigh=iimax-1;jlow=jjmin+1;jhigh=jjmax-1
         ilow=iimin;ihigh=iimax-1;jlow=jjmin;jhigh=jjmax
      case(HV_TAG, V_TAG , DV_TAG) ! for variables defined on v-grid
         ilow=iimin+1;ihigh=iimax-1;jlow=jjmin;jhigh=jjmax-1
         ilow=iimin;ihigh=iimax;jlow=jjmin;jhigh=jjmax-1
      case default                 ! for variables defined on scalar-grid
         ilow=iimin;ihigh=iimax;jlow=jjmin;jhigh=jjmax
   end select
   if (nprocs .eq. 1) then
#if 1
      forall(i=ilow:ilow,j=jlow:jhigh, mask(i,j) .ne. 0) &
         f1(i-1,j,0:kmax) = f2(i,j,0:kmax)
      forall(i=ihigh:ihigh,j=jlow:jhigh, mask(i,j) .ne. 0) &
         f1(i+1,j,0:kmax) = f2(i,j,0:kmax)
      forall(i=ilow:ihigh,j=jlow:jlow, mask(i,j) .ne. 0) &
         f1(i,j-1,0:kmax) = f2(i,j,0:kmax)
      forall(i=ilow:ihigh,j=jhigh:jhigh, mask(i,j) .ne. 0) &
         f1(i,j+1,0:kmax) = f2(i,j,0:kmax)
#else
      f1(ilow -1, : , : )  = f2(ilow , : , :  )
      f1(ihigh+1, : , : )  = f2(ihigh, : , :  )
      f1( : , jlow -1, : ) = f2( : , jlow , : )
      f1( : , jhigh+1, : ) = f2( : , jhigh, : )
#endif
   else
      call update_3d_halo_mpi(f1,f2,iimin,jjmin,iimax,jjmax,kmax,tag)
   end if
   return
   end subroutine update_3d_halo
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: wait_halo - waits for any un-finished communications
!
! !INTERFACE:
   subroutine wait_halo(tag)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Print information on the MPI environment
!
! !INPUT PARAMETERS:
   integer, intent(in)	:: tag
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
!EOP
!-------------------------------------------------------------------------
!BOC
   if (nprocs .gt. 1) then
      call wait_halo_mpi(tag)
   end if
   return
   end subroutine wait_halo
!EOC

!-----------------------------------------------------------------------

   end module halo_zones

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
