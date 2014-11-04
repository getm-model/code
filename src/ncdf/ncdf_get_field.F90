#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdf_get_field()
!
! !INTERFACE:
   module ncdf_get_field
!
! !DESCRIPTION:
!  Provides 2 subroutines for reading 2D and 3D fields from NetCDF files.
!  Vertical interpolation to the model grid is done for 3D fields.
!
! !USES:
   use netcdf
   use exceptions
   IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
   public get_2d_field_ncdf, get_3d_field_ncdf
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: get_2d_field_ncdf()
!
! !INTERFACE:
   subroutine get_2d_field_ncdf(fn,varname,il,ih,jl,jh,break_on_missing,field)
! !USES:
    IMPLICIT NONE
!
! !DESCRIPTION:
!  A two-dimensional netCDF variable with specified global range
! {\tt il < i < ih} and {\tt jl < j < jh} is read into {\tt field}.
! It is checked if the sizes of the fields correspond exactly.
! When calling this funtions, remember that  FORTRAN netCDF variables
! start with index 1.
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn,varname
   integer,          intent(in)        :: il,ih,jl,jh
   logical, intent(in)                 :: break_on_missing
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: field(:,:)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding, Lars Umlauf
!
! !LOCAL VARIABLES:
   integer, dimension(2)               :: start
   integer, dimension(2)               :: edges
   integer, dimension(2)               :: ubounds
   integer                             :: status,ncid,varid
!EOP
!-------------------------------------------------------------------------
   LEVEL3 'get_2d_field_ncdf()'

   start(1) = il
   start(2) = jl
   edges(1) = ih-il+1
   edges(2) = jh-jl+1

   ubounds =  ubound(field)

   if ((ubounds(1) .ne. edges(1)) .or. ubounds(2) .ne. edges(2) ) then
      call getm_error("get_2d_field_ncdf()", &
           "Array bounds inconsistent.")
   endif

   status = nf90_open(trim(fn),NF90_NOWRITE,ncid)
   if (status .NE. NF90_NOERR) then
      call netcdf_error(status,"get_2d_field_ncdf()", &
           "Error opening file "//trim(fn))
   end if

   status = nf90_inq_varid(ncid,trim(varname),varid)
   if (status .ne. NF90_NOERR) then
      if (break_on_missing) then
         call netcdf_error(status,"get_2d_field_ncdf()", &
              "Error inquiring "//trim(varname))
      else
         LEVEL4 trim(varname),': does not exist - continuing'
         return
      end if
   endif

   status = nf90_get_var(ncid,varid,field,start,edges)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"get_2d_field_ncdf()", &
           "Error reading "//trim(varname))
   endif

  status = nf90_close(ncid)
   if (status .ne. NF90_NOERR) then
      call netcdf_error(status,"get_2d_field_ncdf()", &
           "Error closing file")
   endif

   return
   end subroutine get_2d_field_ncdf
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_3d_field_ncdf -
!
! !INTERFACE:
   subroutine get_3d_field_ncdf(fname,var,nf,break_on_missing,f)
!
! !DESCRIPTION:
!  From a NetCDF files - fname - read the variable - var - into the field - f.
!
! !USES:
!   use netcdf
   use domain, only: imin,jmin,imax,jmax,kmax,iextr,jextr,ioff,joff
   use domain, only: il_domain=>il,ih_domain=>ih,jl_domain=>jl,jh_domain=>jh
   use domain, only: H,az
#ifndef NO_3D
   use variables_3d, only: hn
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fname,var
   integer, intent(in)                 :: nf
   logical, intent(in)                 :: break_on_missing
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(inout)             :: f(I3DFIELD)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
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

   LEVEL3 'get_3d_field_ncdf'

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

   il = max(imin+ioff,il_domain); ih = min(imax+ioff,ih_domain)
   jl = max(jmin+joff,jl_domain); jh = min(jmax+joff,jh_domain)
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
   if (break_on_missing) then
      if (err .NE. NF90_NOERR) go to 10
   else
      if (err .NE. NF90_NOERR) then
         LEVEL4 trim(var),': does not exist - continuing'
         return
      end if
   end if

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

#ifndef NO_3D
      call kbk_interpol(kh,zax_2d,wrk_2d,imin,jmin,imax,jmax,kmax, &
                        az,H,hn,f)
#else
      FATAL 'vertical interpolation impossible for NO_3D'
      stop 'get_3d_field_ncdf'
#endif

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
   end subroutine get_3d_field_ncdf
!EOC

!-----------------------------------------------------------------------

   end module ncdf_get_field

!-----------------------------------------------------------------------
! Copyright (C) 2013 - Hans Burchard and Karsten Bolding (BB)          !
!-----------------------------------------------------------------------
