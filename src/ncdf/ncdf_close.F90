!$Id: ncdf_close.F90,v 1.4 2006-02-02 17:51:37 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:  ncdf_close() - closes the specified NetCDF file.
!
! !INTERFACE:
   subroutine ncdf_close()
!
! !DESCRIPTION:
!
! !USES:
   use ncdf_2d, only: nc2d => ncid
#ifndef NO_3D
   use ncdf_3d, only: nc3d => ncid
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
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ncdf_close.F90,v $
!  Revision 1.4  2006-02-02 17:51:37  kbk
!  do not try and save to un-opened NetCDF files
!
!  Revision 1.3  2003/04/23 11:54:03  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 12:46:06  kbk
!  NO_3D
!
!  Revision 1.1.1.1  2002/05/02 14:01:47  gotm
!  recovering after CVS crash
!
!  Revision 1.2  2001/09/13 14:56:58  bbh
!  Also updated
!
!  Revision 1.1.1.1  2001/04/17 08:43:07  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
   integer                   :: err
   REALTYPE                  :: dummy=-_ONE_
!EOP
!-----------------------------------------------------------------------
!BOC
   include "netcdf.inc"
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'ncdf_close() # ',Ncall
#endif

   if (nc2d .ge. 0) then
      call save_2d_ncdf(dummy)
      err = nf_close(nc2d)
         call save_2d_ncdf(dummy)
      if (err .NE. NF_NOERR) go to 10
   end if
#ifndef NO_3D
   if (nc3d .ge. 0) then
      err = nf_close(nc3d)
      if (err .NE. NF_NOERR) go to 10
   end if
#endif
   return

10 FATAL 'ncdf_close: ',nf_strerror(err)
   stop

#ifdef DEBUG
   write(debug,*) 'Leaving ncdf_close()'
   write(debug,*)
#endif
   return
   end subroutine ncdf_close
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
