!$Id: parallel.F90,v 1.3 2003-04-01 15:31:11 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: commhalo - mpi interface to 'getm'
!
! !INTERFACE:
   module commhalo

! !DESCRIPTION:
!  This module is included only to supply myid and nprocs used in various
!  places in 'getm'. From version 1.4 real use of MPI will be implemented.
!
! !USES:
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   integer		:: myid, nprocs
   integer		:: comm_hd=-1
   integer, parameter	:: H_TAG=10,HU_TAG=11,HV_TAG=12
   integer, parameter	:: D_TAG=20,DU_TAG=21,DV_TAG=22
   integer, parameter   :: z_TAG=30,U_TAG=31,V_TAG=32

!  Different mesh specification methods
   integer, parameter	:: ONE_CELL=-1
!  Methods of communication
   integer, parameter	:: ONE_PROCESS=-1
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: parallel.F90,v $
!  Revision 1.3  2003-04-01 15:31:11  gotm
!  removed dead code
!
!  Revision 1.2  2003/03/24 14:21:11  gotm
!  corrected boundary indices
!
!  Revision 1.1.1.1  2002/05/02 14:01:29  gotm
!  recovering after CVS crash
!
!  Revision 1.4  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.3  2001/05/18 12:53:08  bbh
!  Prepared for mask in update_2d_halo - but not used yet
!
!  Revision 1.2  2001/05/18 10:03:44  bbh
!  Added mask in parameter list to update_3d_halo()
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
   integer, private	:: active_comm=-1
   integer, private	:: comm_method=ONE_PROCESS
!EOP
!-----------------------------------------------------------------------
!BOC

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_mpi - initialize MPI environment
!
! !INTERFACE:
   subroutine init_mpi()
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Initialize Parallel environment
!
! !REVISION HISTORY:
!
!  22Apr99   Karsten Bolding & Hans Burchard  Initial code.
!
! !LOCAL VARIABLES:
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_mpi'
#endif

   myid   = -1
   nprocs =  1

#ifdef DEBUG
   write(debug,*) 'Leaving init_mpi()'
   write(debug,*)
#endif

   return
   end subroutine init_mpi

!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: print_mpi_info - write various MPI related info.
!
! !INTERFACE:
   subroutine print_mpi_info()
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Print information on the MPI environment
!
! !REVISION HISTORY:
!
!  22Apr99   Karsten Bolding & Hans Burchard  Initial code.
!
! !LOCAL VARIABLES:
!
!EOP
!-------------------------------------------------------------------------
!BOC
   LEVEL1 'OK - you specified a parallel run'
   LEVEL1 'prallel execution is not available in this version'
   LEVEL1 'setting myid and nprocs to useable values'

   end subroutine print_mpi_info
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
! !INPUT PARAMETERS:
   integer, intent(in)		:: imin,jmin,imax,jmax
   integer, intent(in)		:: tag
   REALTYPE, intent(in) 	:: f2(E2DFIELD)
   integer, intent(in) 		:: mask(E2DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)	:: f1(E2DFIELD)
!
! !DESCRIPTION:
!  Print information on the MPI environment
!
! !REVISION HISTORY:
!
!  22Apr99   Karsten Bolding & Hans Burchard  Initial code.
!
! !LOCAL VARIABLES:
   integer	:: i,j,k
   integer	:: ilow,jlow,ihigh,jhigh
!EOP
!-------------------------------------------------------------------------
!BOC
#if 0
   select case (tag)
      case(HU_TAG, U_TAG , DU_TAG) ! for variables defined on u-grid
         ilow=imin;ihigh=imax-1;jlow=jmin;jhigh=jmax
      case(HV_TAG, V_TAG , DV_TAG) ! for variables defined on v-grid
         ilow=imin;ihigh=imax;jlow=jmin;jhigh=jmax-1
      case default                 ! for variables defined on scalar-grid
         ilow=imin;ihigh=imax;jlow=jmin;jhigh=jmax
   end select
#endif

   ilow=imin;ihigh=imax;jlow=jmin;jhigh=jmax

   select case (comm_method)
      case(ONE_PROCESS)
         f1(ilow -1, : )  = f2(ilow,  :  )
         f1(ihigh+1, : )  = f2(ihigh, :  )
         f1( :, jlow -1 ) = f2( :, jlow  )
         f1( :, jhigh+1 ) = f2( :, jhigh )
      case default
         FATAL 'A non valid communication method has been chosen'
         stop 'update_2d_halo'
   end select
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
! !INPUT PARAMETERS:
   integer, intent(in)	:: iimin,jjmin,iimax,jjmax,kmax
   integer, intent(in)	:: tag
   integer, intent(in)	:: mask(-HALO+1:,-HALO+1:)
   REALTYPE, intent(in)	:: f2(I3DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out):: f1(I3DFIELD)
!
! !DESCRIPTION:
!  Print information on the MPI environment
!
! !REVISION HISTORY:
!
!  22Apr99   Karsten Bolding & Hans Burchard  Initial code.
!
! !LOCAL VARIABLES:
   integer	:: i,j,k
   integer	:: ilow,jlow,ihigh,jhigh
!EOP
!-------------------------------------------------------------------------
!BOC

#if 0
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
#endif

   ilow=iimin;ihigh=iimax;jlow=jjmin;jhigh=jjmax

   select case (comm_method)
      case(ONE_PROCESS)
         f1(ilow -1, : , : )  = f2(ilow , : , :  )
         f1(ihigh+1, : , : )  = f2(ihigh, : , :  )
         f1( : , jlow -1, : ) = f2( : , jlow , : )
         f1( : , jhigh+1, : ) = f2( : , jhigh, : )
      case default
         FATAL 'A non valid communication method has been chosen'
         stop 'update_3d_halo'
   end select
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
!  Waits for any un-finished communications. In the case of serial model run,
!  only 1 process or blocking communication this routine does not do anything.
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
!   if (nprocs .gt. 1) then
!      call wait_halo_mpi(tag)
!   end if
   return
   end subroutine wait_halo
!EOC

!-----------------------------------------------------------------------

   end module commhalo

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
