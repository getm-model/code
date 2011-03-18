#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: read_field_ncdf -
!
! !INTERFACE:
   subroutine read_field_ncdf(fname,var,nf,f)
!
! !DESCRIPTION:
!  From a NetCDF files - fname - read the variable - var - into the field - f.
!
! !USES:
   use netcdf
   use domain, only: imin,jmin,imax,jmax,kmax,iextr,jextr,ioff,joff
   use domain, only: H,az
   use variables_3d, only: hn
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fname,var
   integer, intent(in)                 :: nf
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(inout)             :: f(I3DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: il,jl,iloc,jloc,indx
   integer                   :: ih,jh,kh,nh
   integer                   :: rc,err,ncid,var_id,i,j,k,n
   integer                   :: start(4),edges(4)
   integer                   :: ndims
   integer                   :: xax_id=-1,yax_id=-1,zax_id=-1,time_id=-1
   character(len=256)        :: dimname
   REAL_4B, allocatable      :: zax(:), tax(:), wrk(:)
   REALTYPE, allocatable     :: zax_2d(:), wrk_2d(:,:,:)
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   write(debug,*) 'read_field_ncdf (NetCDF)'
   write(debug,*) 'Reading from: ',trim(fname)
#endif

   LEVEL3 'read_field_ncdf'

   err = nf90_open(fname,NF90_NOWRITE,ncid)
   if (err .NE. NF90_NOERR) go to 10

   err = nf90_inquire(ncid,nDimensions=ndims)
   if (err .ne. NF90_NOERR) go to 10

   if(ndims .ne. 4) then
      FATAL 'we assume to read from a 4D field - and ndims = ',ndims
      stop 'read_field_ncdf()'
   end if

   do n=1,ndims
      err = nf90_inquire_dimension(ncid,n,name=dimname)
      if (err .ne. NF90_NOERR) go to 10

#if 0
      if( dimname .eq. 'lon' ) then
         xax_id = n
         err = nf90_inq_dimlen(ncid,xax_id,iextr)
         if (err .ne. NF90_NOERR) go to 10
!         LEVEL4 'xax_id  --> ',xax_id,', len = ',iextr
      end if
      if( dimname .eq. 'lat' ) then
         yax_id = n
         err = nf90_inq_dimlen(ncid,yax_id,jextr)
         if (err .ne. NF90_NOERR) go to 10
!         LEVEL4 'yax_id  --> ',yax_id,', len = ',jextr
      end if
#endif

      if( dimname .eq. 'zax' ) then
         zax_id = n
         err = nf90_inquire_dimension(ncid,zax_id,len = kh)
         if (err .ne. NF90_NOERR) go to 10
!         LEVEL4 'zax_id  --> ',zax_id,', len = ',kh
      end if
      if( dimname .eq. 'time' ) then
         time_id = n
         err = nf90_inquire_dimension(ncid,time_id,len = nh)
         if (err .ne. NF90_NOERR) go to 10
!         LEVEL4 'time_id --> ',time_id,', len = ',nh
      end if
   end do

#if 0
   if(xax_id .eq. -1) then
      FATAL 'could not find x-axis information in ',trim(fname)
      stop 'read_field_ncdf()'
   end if
   if(yax_id .eq. -1) then
      FATAL 'could not find y-axis information in ',trim(fname)
      stop 'read_field_ncdf()'
   end if
#endif

   if(zax_id .eq. -1) then
      FATAL 'could not find z-axis information in ',trim(fname)
      stop 'read_field_ncdf()'
   end if
   if(time_id .eq. -1) then
      FATAL 'could not find time-axis information in ',trim(fname)
      stop 'read_field_ncdf()'
   end if

   if (kh .gt. 1) then
      err = nf90_inq_varid(ncid,"zax",var_id)
      if (err .NE. NF90_NOERR) go to 10
      allocate(zax(kh),stat=rc)
      if (rc /= 0) stop 'read_field_ncdf: Error allocating zax'
      err = nf90_get_var(ncid,var_id,zax)
      if (err .NE. NF90_NOERR) go to 10
   end if

#if 0
   err = nf90_inq_varid(ncid,"time",var_id)
   if (err .NE. NF90_NOERR) go to 10
   allocate(tax(nh),stat=rc)
   if (rc /= 0) stop 'read_field_ncdf: Error allocating tax'
   err = nf90_get_var(ncid,var_id,tax)
   if (err .NE. NF90_NOERR) go to 10
#endif

   il = max(imin+ioff,1); ih = min(imax+ioff,iextr)
   jl = max(jmin+joff,1); jh = min(jmax+joff,jextr)
   iloc = max(imin-ioff,1); jloc = max(jmin-joff,1)

   start(1) = il ; start(2) = jl
   start(3) = 1; 
   edges(1) = ih-il+1 ; edges(2) = jh-jl+1
   edges(3) = kh; edges(4) = 1

   start(4) = -1
   do n=1,nh
      if (n .eq. nf) then
         start(4) = n
         EXIT
      end if
   end do

   if(start(4) .lt. 0) then
      FATAL 'could not find requested field: ',nf
      stop 'read_field_ncdf()' 
   end if

   allocate(wrk((imax-imin+1)*(jmax-jmin+1)*kh),stat=rc)
   if (rc /= 0) stop 'read_field_ncdf: Error allocating wrk'

   err = nf90_inq_varid(ncid,trim(var),var_id)
   if (err .NE. NF90_NOERR) go to 10

   err = nf90_get_var(ncid,var_id,wrk,start,edges)
   if (err .NE. NF90_NOERR) go to 10

   if (kh .gt. 1) then
      allocate(zax_2d(kh),stat=rc)
      if (rc /= 0) stop 'read_field_ncdf: Error allocating zax_2d'
      do k=1,kh
         zax_2d(k) = -zax(kh-k+1)
      end do

      allocate(wrk_2d(imin:imax,jmin:jmax,kh),stat=rc)
      if (rc /= 0) stop 'read_field_ncdf: Error allocating wrk_2d'
      indx = 1
      do k=1,kh
         do j=jl,jh
            do i=il,ih
               wrk_2d(i-il+iloc,j-jl+jloc,k) = wrk(indx)
               indx = indx+1
            end do
         end do
      end do

      call kbk_interpol(kh,zax_2d,wrk_2d,imin,jmin,imax,jmax,kmax, &
                        az,H,hn,f)
   else
      indx = 1
      do j=jl,jh
         do i=il,ih
            f(i-il+iloc,j-jl+jloc,:) = wrk(indx)
            indx = indx+1
         end do
      end do
   end if

   err = nf90_close(ncid)
   if (err .NE. NF90_NOERR) go to 10

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
10 FATAL 'read_field_ncdf: ',nf90_strerror(err)
   stop
   end subroutine read_field_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
