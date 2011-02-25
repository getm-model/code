!$Id: get_field_ncdf.F90,v 1.3 2007-06-07 10:25:19 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: get_field_ncdf -
!
! !INTERFACE:
   subroutine get_field_ncdf(fname,var,f)
!
! !DESCRIPTION:
!  From a NetCDF files - fname - read the variable - var - into the field - f.
!
! !USES:
   use netcdf
   use domain, only: imin,imax,jmin,jmax,kmax
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fname,var
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: f(I3DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: get_field_ncdf.F90,v $
!  Revision 1.3  2007-06-07 10:25:19  kbk
!  iimin,iimax,jjmin,jjmax -> imin,imax,jmin,jmax
!
!  Revision 1.2  2003-04-23 11:54:03  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.1.1.1  2002/05/02 14:01:48  gotm
!  recovering after CVS crash
!
!  Revision 1.3  2001/10/22 08:10:43  bbh
!  De-allocate wrk #ifdef FORTRAN90
!
!  Revision 1.1  2001/05/10 11:38:29  bbh
!  Added get_field_ncdf() + various small bug fixes
!
! !LOCAL VARIABLES:
   integer                   :: rc,err,ncid,var_id,i,j,k,size,indx
   integer                   :: start(3),edges(3)
   REAL_4B, allocatable      :: wrk(:)
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   write(debug,*) 'get_field_ncdf (NetCDF)'
   write(debug,*) 'Reading from: ',trim(fname)
#endif

   LEVEL3 'get_field_ncdf'

   err = nf90_open(fname,NF90_NOWRITE,ncid)
   if (err .NE. NF90_NOERR) go to 10
#if 0
   err = nf90_inq_unlimdim(ncid, rec_id)
   if (err .ne. NF90_NOERR) go to 10

   err = nf90_inq_dimlen(ncid, rec_id, nsets)
   if (err .ne. NF90_NOERR) go to 10

   err = nf90_inq_dimid(ncid, 'nbdyp', bdy_id)
   if (err .ne. NF90_NOERR)  go to 10

   err = nf90_inq_dimlen(ncid, bdy_id, bdy_len)
   if (err .ne. NF90_NOERR) go to 10
#endif
   err = nf90_inq_varid(ncid,trim(var),var_id)
   if (err .NE. NF90_NOERR) go to 10

   size = (imax-imin+1)*(jmax-jmin+1)*(kmax+1)
   allocate(wrk(size),stat=rc)
   if (rc /= 0) stop 'get_field_ncdf: Error allocating work-space'

   start(1) = imin ; start(2) = jmin; start(3) = 1
   edges(1) = imax-imin+1 ; edges(2) = jmax-jmin+1; edges(3) = kmax
   err = nf90_get_var(ncid,var_id,wrk,start,edges)
   if (err .NE. NF90_NOERR) go to 10

   indx = 1
   do k=1,kmax
      do j=jmin,jmax
         do i=imin,imax
            f(i,j,k) = wrk(indx)
            indx = indx+1
         end do
      end do
   end do

   err = nf90_close(ncid)
   if (err .NE. NF90_NOERR) go to 10

#ifdef DEBUG
   write(debug,*) 'Leaving get_field_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'get_field_ncdf: ',nf90_strerror(err)
   stop
   end subroutine get_field_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
