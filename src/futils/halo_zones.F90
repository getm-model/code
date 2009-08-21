!$Id: halo_zones.F90,v 1.7 2009-08-21 08:56:34 bjb Exp $
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
#ifdef GETM_PARALLEL
   use halo_mpi
#endif
   IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_halo_zones,update_2d_halo,update_3d_halo,wait_halo
!
! !PUBLIC DATA MEMBERS:
#ifndef GETM_PARALLEL
   integer, parameter                  :: H_TAG=10,HU_TAG=11,HV_TAG=12
   integer, parameter                  :: D_TAG=20,DU_TAG=21,DV_TAG=22
   integer, parameter                  :: z_TAG=30,U_TAG=31,V_TAG=32
#endif
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: halo_zones.F90,v $
!  Revision 1.7  2009-08-21 08:56:34  bjb
!  Fix name clash on PARALLEL with OpenMP key word
!
!  Revision 1.6  2008-09-16 10:03:24  kb
!  added Holtermanns emergency break algorithm
!
!  Revision 1.5  2007-06-07 10:25:19  kbk
!  iimin,iimax,jjmin,jjmax -> imin,imax,jmin,jmax
!
!  Revision 1.4  2003-08-03 08:49:51  kbk
!  cleaned code
!
!  Revision 1.3  2003/05/09 11:52:08  kbk
!  do not mirror coordinate info + use mask for inverse area calculation
!
!  Revision 1.2  2003/04/23 12:02:43  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.1  2003/04/07 12:05:42  kbk
!  new parallel related files
!
!  Revision 1.1.1.1  2002/05/02 14:01:29  gotm
!  recovering after CVS crash
!
! !LOCAL VARIABLES:
#ifndef GETM_PARALLEL
   integer, parameter        :: nprocs=1
#endif
!EOP
!-----------------------------------------------------------------------
!BOC

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_halo_zones - 
!
! !INTERFACE:
   subroutine init_halo_zones()
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Initialize Parallel environment
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
! !LOCAL VARIABLES:
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_halo_zones'
#endif

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
   subroutine update_2d_halo(f1,f2,mask,imin,jmin,imax,jmax,tag,mirror)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Print information on the MPI environment
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: imin,jmin,imax,jmax
   integer, intent(in)                 :: mask(-HALO+1:,-HALO+1:)
   integer, intent(in)                 :: tag
   logical, optional, intent(in)       :: mirror
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout)             :: f1(E2DFIELD),f2(E2DFIELD)
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   integer                   :: il,jl,ih,jh
   logical                   :: do_mirror=.true.
!EOP
!-------------------------------------------------------------------------
!BOC
   il=imin;ih=imax;jl=jmin;jh=jmax

   if ( present(mirror) ) do_mirror = mirror

   if (nprocs .eq. 1) then
      if ( do_mirror ) then
      f1(il-1, : )  = f2(il,  :  )
      f1(ih+1, : )  = f2(ih, :  )
      f1( :, jl-1 ) = f2( :, jl  )
      f1( :, jh+1 ) = f2( :, jh )
      f1(il-1,jh+1) = f2(il,jh)
      f1(ih+1,jh+1) = f2(ih,jh)
      f1(ih+1,jl-1) = f2(ih,jl)
      f1(il-1,jl-1) = f2(il,jl)
      end if
   else
#ifdef GETM_PARALLEL
      if( present(mirror) ) then
         call update_2d_halo_mpi(f1,f2,imin,jmin,imax,jmax,tag, &
                                 mirror=mirror)
      else
         call update_2d_halo_mpi(f1,f2,imin,jmin,imax,jmax,tag)
      endif
#endif
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
   subroutine update_3d_halo(f1,f2,mask,imin,jmin,imax,jmax,kmax,tag)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Print information on the MPI environment
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: imin,jmin,imax,jmax,kmax
   integer, intent(in)                 :: tag
   integer, intent(in)                 :: mask(-HALO+1:,-HALO+1:)
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
   integer                   :: i,j,k
   integer                   :: il,jl,ih,jh
!EOP
!-------------------------------------------------------------------------
!BOC
   il=imin;ih=imax;jl=jmin;jh=jmax

   if (nprocs .eq. 1) then
      f1(il-1, : , : )  = f2(il, : , :  )
      f1(ih+1, : , : )  = f2(ih, : , :  )
      f1( : , jl-1, : ) = f2( : , jl, : )
      f1( : , jh+1, : ) = f2( : , jh, : )
   else
#ifdef GETM_PARALLEL
      call update_3d_halo_mpi(f1,f2,imin,jmin,imax,jmax,kmax,tag)
#endif
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
   integer, intent(in)                 :: tag
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
#ifdef GETM_PARALLEL
   if (nprocs .gt. 1) then
      call wait_halo_mpi(tag)
   end if
#endif
   return
   end subroutine wait_halo
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_flag - set a flag across sub-domains
!
! !INTERFACE:
   subroutine set_flag(n,flag,flags)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Call the parallel library (MPI) routine for setting an array of flags
!  - or - in case of a serial run just set the local value. 
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: n,flag
!
! !OUTPUT PARAMETERS:
   integer, intent(out)                :: flags(n)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL VARIABLES:
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef GETM_PARALLEL
   if (nprocs .gt. 1) then
      call set_flag_mpi(n,flag,flags)
   end if
#else
   flags(n) = flag
#endif
   return
   end subroutine set_flag
!EOC

!-----------------------------------------------------------------------

   end module halo_zones

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
