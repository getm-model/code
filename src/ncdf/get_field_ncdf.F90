!$Id: get_field_ncdf.F90,v 1.1.1.1 2002-05-02 14:01:48 gotm Exp $
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
   use domain, only: iimin,iimax,jjmin,jjmax,kmax
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)	:: fname,var
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)	:: f(I3DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: get_field_ncdf.F90,v $
!  Revision 1.1.1.1  2002-05-02 14:01:48  gotm
!  recovering after CVS crash
!
!  Revision 1.3  2001/10/22 08:10:43  bbh
!  De-allocate wrk #ifdef FORTRAN90
!
!  Revision 1.1  2001/05/10 11:38:29  bbh
!  Added get_field_ncdf() + various small bug fixes
!
! !LOCAL VARIABLES:
   integer	:: rc,err,ncid,var_id,i,j,k,size,indx
   integer	:: start(3),edges(3)
   REAL_4B, allocatable	:: wrk(:)
!EOP
!-------------------------------------------------------------------------
!BOC
   include "netcdf.inc"
#ifdef DEBUG
   write(debug,*) 'get_field_ncdf (NetCDF)'
   write(debug,*) 'Reading from: ',trim(fname)
#endif

   LEVEL3 'get_field_ncdf'

   err = nf_open(fname,NCNOWRIT,ncid)
   if (err .NE. NF_NOERR) go to 10
#if 0
   err = nf_inq_unlimdim(ncid, rec_id)
   if (err .ne. NF_NOERR) go to 10

   err = nf_inq_dimlen(ncid, rec_id, nsets)
   if (err .ne. NF_NOERR) go to 10

   err = nf_inq_dimid(ncid, 'nbdyp', bdy_id)
   if (err .ne. NF_NOERR)  go to 10

   err = nf_inq_dimlen(ncid, bdy_id, bdy_len)
   if (err .ne. NF_NOERR) go to 10
#endif
   err = nf_inq_varid(ncid,trim(var),var_id)
   if (err .NE. NF_NOERR) go to 10

   size = (iimax-iimin+1)*(jjmax-jjmin+1)*(kmax+1)
   allocate(wrk(size),stat=rc)
   if (rc /= 0) stop 'get_field_ncdf: Error allocating work-space'

   start(1) = iimin ; start(2) = jjmin; start(3) = 1
   edges(1) = iimax-iimin+1 ; edges(2) = jjmax-jjmin+1; edges(3) = kmax
   err = nf_get_vara_real(ncid,var_id,start,edges,wrk)
   if (err .NE. NF_NOERR) go to 10

   indx = 1
   do k=1,kmax
      do j=jjmin,jjmax
         do i=iimin,iimax
            f(i,j,k) = wrk(indx)
            indx = indx+1
         end do
      end do
   end do

   err = nf_close(ncid)
   if (err .NE. NF_NOERR) go to 10

#ifdef DEBUG
   write(debug,*) 'Leaving get_field_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'get_field_ncdf: ',nf_strerror(err)
   stop
   end subroutine get_field_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
