!$Id: read_field_ncdf.F90,v 1.2 2003-04-07 12:39:59 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: read_field_ncdf -
!
! !INTERFACE:
   subroutine read_field_ncdf(fname,var,n,f)
!
! !DESCRIPTION:
!  From a NetCDF files - fname - read the variable - var - into the field - f.
!
! !USES:
   use domain, only: imin,jmin,imax,jmax,ioff,joff
   use domain, only: iimin,jjmin,iimax,jjmax,kmax
   use domain, only: H,az
   use variables_3d, only: hn
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)	:: fname,var
   integer			:: n
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(inout)	:: f(I3DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: read_field_ncdf.F90,v $
!  Revision 1.2  2003-04-07 12:39:59  kbk
!  parallel support
!
!  Revision 1.1.1.1  2002/05/02 14:01:47  gotm
!  recovering after CVS crash
!
!
! !LOCAL VARIABLES:
   integer	:: ih,jh,kh,nh
   integer	:: rc,err,ncid,var_id,i,j,k
   integer	:: start(4),edges(4)
   integer	:: ndims,rec_id
   REAL_4B, allocatable		:: zax(:), wrk(:,:,:)
   REALTYPE, allocatable	:: zax_2d(:), wrk_2d(:,:,:)
!EOP
!-------------------------------------------------------------------------
!BOC
   include "netcdf.inc"
#ifdef DEBUG
   write(debug,*) 'read_field_ncdf (NetCDF)'
   write(debug,*) 'Reading from: ',trim(fname)
#endif

   LEVEL3 'read_field_ncdf'

   err = nf_open(fname,NCNOWRIT,ncid)
   if (err .NE. NF_NOERR) go to 10

   err = nf_inq_ndims(ncid, ndims)
   if (err .ne. NF_NOERR) go to 10

!  We assume dimension order as follows 1->i, 2->j, 3->k, 4->n

   err = nf_inq_dimlen(ncid, 1, ih)
   if (err .ne. NF_NOERR) go to 10
   err = nf_inq_dimlen(ncid, 2, jh)
   if (err .ne. NF_NOERR) go to 10
   err = nf_inq_dimlen(ncid, 3, kh)
   if (err .ne. NF_NOERR) go to 10
   err = nf_inq_dimlen(ncid, 4, nh)
   if (err .ne. NF_NOERR) go to 10

   err = nf_inq_varid(ncid,"zax",var_id)
   if (err .NE. NF_NOERR) go to 10
   allocate(zax(kh),stat=rc)
   if (rc /= 0) stop 'read_field_ncdf: Error allocating zax'
   err = nf_get_var_real(ncid,var_id,zax)
   if (err .NE. NF_NOERR) go to 10

   err = nf_inq_varid(ncid,trim(var),var_id)
   if (err .NE. NF_NOERR) go to 10

   allocate(wrk(ih,jh,kh),stat=rc)
   if (rc /= 0) stop 'read_field_ncdf: Error allocating wrk'

   start(1) = iimin+ioff ; start(2) = jjmin+joff; 
   start(3) = 1; start(4) = n;
   edges(1) = iimax-iimin+1 ; edges(2) = jjmax-jjmin+1; 
   edges(3) = kh; edges(4) = 1

   err = nf_get_vara_real(ncid,var_id,start,edges,wrk)
   if (err .NE. NF_NOERR) go to 10

   allocate(zax_2d(kh),stat=rc)
   if (rc /= 0) stop 'read_field_ncdf: Error allocating zax_2d'
   do k=1,kh
      zax_2d(k) = -zax(kh-k+1)
   end do
   allocate(wrk_2d(ih,jh,kh),stat=rc)
   if (rc /= 0) stop 'read_field_ncdf: Error allocating wrk_2d'
   wrk_2d = wrk

   call kbk_interpol(kh,zax_2d,wrk_2d,imin,jmin,imax,jmax,az,H, &
                     iimin,jjmin,iimax,jjmax,kmax,hn,f)

   err = nf_close(ncid)
   if (err .NE. NF_NOERR) go to 10

#ifdef FORTRAN90
   deallocate(zax,stat=rc)
   if (rc /= 0) stop 'read_field_ncdf: Error de-allocating memory (zax)'
   deallocate(zax_2d,stat=rc)
   if (rc /= 0) stop 'read_field_ncdf: Error de-allocating memory (zax_2d)'
   deallocate(wrk,stat=rc)
   if (rc /= 0) stop 'read_field_ncdf: Error de-allocating memory (wrk)'
   deallocate(wrk_2d,stat=rc)
   if (rc /= 0) stop 'read_field_ncdf: Error de-allocating memory (wrk_2d)'
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving read_field_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'read_field_ncdf: ',nf_strerror(err)
   stop
   end subroutine read_field_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
