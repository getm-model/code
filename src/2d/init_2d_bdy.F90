!$Id: init_2d_bdy.F90,v 1.1.1.1 2002-05-02 14:00:44 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_2d_bdy() - checks whether this node has boundaries.
!
! !INTERFACE:
   subroutine init_2d_bdy(bdyfmt,FName)
!
! !DESCRIPTION:
!
! !USES:
   use domain,     only: iextr,jextr,NWB,NNB,NEB,NSB
   use m2d,        only: EWBdy,ENBdy,EEBdy,ESBdy
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
   integer, intent(in)	:: bdyfmt
   character(len=*), intent(in)	:: FName
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: init_2d_bdy.F90,v $
!  Revision 1.1.1.1  2002-05-02 14:00:44  gotm
!  recovering after CVS crash
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
   integer	:: rc
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_2d_bdy() # ',Ncall
#endif

#ifndef STATIC
   if (NWB .ge. 1) then
      allocate(EWBdy(jextr),stat=rc)
      if (rc /= 0) stop 'init_2d_bdy: Error allocating memory (EWBdy)'
   end if
   if (NNB .ge. 1) then
      allocate(ENBdy(iextr),stat=rc)
      if (rc /= 0) stop 'init_2d_bdy: Error allocating memory (ENBdy)'
   end if
   if (NEB .ge. 1) then
      allocate(EEBdy(jextr),stat=rc)
      if (rc /= 0) stop 'init_2d_bdy: Error allocating memory (EEBdy)'
   end if
   if (NSB .ge. 1) then
      allocate(ESBdy(iextr),stat=rc)
      if (rc /= 0) stop 'init_2d_bdy: Error allocating memory (ESBdy)'
   end if
#endif

   select case (bdyfmt)
      case (ASCII)
      case (NETCDF)
!kbk         call ncdf_init_2d_bdy(fname)
      case default
        STDERR 'Fatal error: init_2d_bdy'
        stop 'init_2d_bdy'
   end select


#ifdef DEBUG
   write(debug,*) 'Leaving init_2d_bdy()'
   write(debug,*)
#endif
   return
   end subroutine init_2d_bdy
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
