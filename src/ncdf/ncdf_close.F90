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
   use netcdf
   use ncdf_2d, only: nc2d => ncid
#ifndef NO_3D
   use ncdf_3d, only: nc3d => ncid
#endif
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: err
   REALTYPE                  :: dummy=-_ONE_
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'ncdf_close() # ',Ncall
#endif

   if (nc2d .ge. 0) then
      call save_2d_ncdf(dummy)
      err = nf90_close(nc2d)
      if (err .NE. NF90_NOERR) go to 10
   end if
#ifndef NO_3D
   if (nc3d .ge. 0) then
      err = nf90_close(nc3d)
      if (err .NE. NF90_NOERR) go to 10
   end if
#endif
   return

10 FATAL 'ncdf_close: ',nf90_strerror(err)
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
