!$Id: ncdf_topo.F90,v 1.4 2003-05-02 08:19:14 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: get_dimensions - reads domain dimensions from NetCDF-file
! !INTERFACE:
#ifdef STATIC
   subroutine get_dimensions(fname,rc)
#else
   subroutine get_dimensions(fname,iextr,jextr,rc)
#endif
!
! !DESCRIPTION:
!  To be written
!
! !USES:
   use ncdfin
   implicit none
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fname
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
#ifdef STATIC
   integer,intent(out)                 :: rc
#else
   integer,intent(out)                 :: iextr,jextr,rc
#endif
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ncdf_topo.F90,v $
!  Revision 1.4  2003-05-02 08:19:14  kbk
!  fixed calculationof lonx andlatx if not in file
!
!  Revision 1.3  2003/04/23 11:54:03  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 12:39:59  kbk
!  parallel support
!
!  Revision 1.1.1.1  2002/05/02 14:01:47  gotm
!  recovering after CVS crash
!
!  Revision 1.5  2001/10/22 12:07:05  bbh
!  Reading curvi-linear specific fields
!
!  Revision 1.4  2001/09/26 10:01:41  bbh
!  lat and lon maps now read in ncdf_topo.F90
!
!  Revision 1.3  2001/09/19 11:20:32  bbh
!  Explicit de-allocates memory when -DFORTRAN90
!
!  Revision 1.2  2001/08/01 08:48:17  bbh
!  CURVILINEAR now implemented - reading additional variables
!
!  Revision 1.1.1.1  2001/04/17 08:43:07  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
   integer                   :: err
!
!EOP
!-------------------------------------------------------------------------
!BOC
   include "netcdf.inc"
#ifdef DEBUG
   write(debug,*) 'get_dimensions (NetCDF)'
   write(debug,*) 'Reading from: ',TRIM(fname)
#endif
   err = nf_open(fname,nf_nowrite,ncbathy)
   if (err .NE. NF_NOERR) go to 10

! get bathymetry and dimension information
#ifndef CURVILINEAR
   err = nf_inq_varid(ncbathy,'bathymetry',h_id)
#else
   err = nf_inq_varid(ncbathy,'h',h_id)
#endif
   if (err .ne. NF_NOERR) go to 10
   err = nf_inq_varndims(ncbathy, h_id, dims)
   if (err .ne. NF_NOERR) go to 10
!kbk      if (dims.ne.2) call fatal(E_TOPO_DIMS)
   err = nf_inq_vardimid(ncbathy,h_id,dimids)
   if (err .ne. NF_NOERR) go to 10

#ifndef STATIC
! get and check x dimension
   err = nf_inq_dimlen(ncbathy,dimids(1),iextr)
   if (err .ne. NF_NOERR) go to 10

! get and check y dimension
   err = nf_inq_dimlen(ncbathy,dimids(2),jextr)
   if (err .ne. NF_NOERR) go to 10
#endif

   rc = 0

#ifdef DEBUG
   write(debug,*) 'Leaving get_dimensions()'
   write(debug,*)
#endif
   return
10 FATAL 'get_dimensions: ',nf_strerror(err)
   stop
   return
   end subroutine get_dimensions
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: get_bathymetry - reads bathymetry from NetCDF file
!
! !INTERFACE:
   subroutine get_bathymetry(H,Hland,iextr,jextr,ioff,joff, &
                             imin,imax,jmin,jmax,rc)
!
! !USES:
   use ncdfin
   use domain, only: lonc,latc
#if ! ( defined(SPHERICAL) || defined(CURVILINEAR) ) 
   use domain, only: dx,dy,conv
#else
   use domain, only: lonx,latx,lonu,latu,lonv,latv
   use domain, only: dxdyc,dydxc,angle,conv
#if defined(CURVILINEAR)
   use domain, only: xx,yx,xc,yc,xu,yu,xv,yv
#endif
#endif
   IMPLICIT NONE
!
! !DESCRIPTION:

! !INPUT PARAMETERS:
   integer, intent(in)                 :: iextr,jextr,ioff,joff
   integer, intent(in)                 :: imin,imax,jmin,jmax
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: H(E2DFIELD)
   REALTYPE, intent(out)               :: Hland
   integer, intent(out)                :: rc
!
! !REVISION HISTORY:
!
!  22Apr99   Karsten Bolding & Hans Burchard  Initial code.
!
! !LOCAL VARIABLES:
   integer                   :: size,i,j,indx
   integer                   :: id,err
   integer                   :: il,ih,jl,jh,iloc,jloc
   REAL_4B, allocatable      :: wrk(:)
   REALTYPE, parameter       :: pi=3.1415927,deg2rad=pi/180.,rad2deg=180./pi
   REALTYPE                  :: x
!EOP
!-------------------------------------------------------------------------
#include"netcdf.inc"
#ifdef DEBUG
   write(debug,*) 'get_bathymetry (NetCDF)'
   write(debug,*) imin,imax,jmin,jmax
#endif
   size = (imax-imin+1)*(jmax-jmin+1)
   allocate(wrk(size),stat=rc)
   if (rc /= 0) stop 'get_bathymetry: Error allocating work-space'

#if ! ( defined(SPHERICAL) || defined(CURVILINEAR) )
   err = nf_get_att_real(ncbathy,h_id,'missing_value',wrk(1))
   if (err .ne. NF_NOERR) go to 10
#else
   wrk(1) = -10.
#endif
   Hland = wrk(1)

   H=Hland

   LEVEL3 'reading bathymetry'

   il = max(imin+ioff,1); ih = min(imax+ioff,iextr)
   jl = max(jmin+joff,1); jh = min(jmax+joff,jextr)
   iloc = max(imin-ioff,1); jloc = max(jmin-joff,1)

   call read_2d_field("bathymetry",imin,imax,jmin,jmax,H, &
                      il,ih,jl,jh,iloc,jloc)
   where ( H .gt. 20000.)
      H = Hland
   end where

! get grid information
#if ! ( defined(SPHERICAL) || defined(CURVILINEAR) ) 
   LEVEL3 'reading cartesian grid information'
   err = nf_inq_varid(ncbathy,'dx',dx_id)
   if (err .ne. NF_NOERR) go to 10
   err = nf_get_var_real(ncbathy,dx_id,wrk(1))
   if (err .ne. NF_NOERR) go to 10
   dx = wrk(1)

   err = nf_inq_varid(ncbathy,'dy',dy_id)
   if (err .ne. NF_NOERR) go to 10
   err = nf_get_var_real(ncbathy,dy_id,wrk(1))
   if (err .ne. NF_NOERR) go to 10
   dy = wrk(1)

   call read_2d_field("lon",imin,imax,jmin,jmax,lonc,il,ih,jl,jh,iloc,jloc)
   call read_2d_field("lat",imin,imax,jmin,jmax,latc,il,ih,jl,jh,iloc,jloc)
   conv = _ZERO_
   call read_2d_field("conv",imin,imax,jmin,jmax,conv,il,ih,jl,jh,iloc,jloc)
#endif

#ifdef SPHERICAL
   LEVEL3 'reading spherical grid information'

!  Reading longitude for C and X points
   il = max(imin+ioff,1); ih = min(imax+ioff,iextr)
   start(1) = il ; edges(1) = ih-il+1
   err = nf_inq_varid(ncbathy,"lon",id)
   if (err .ne. NF_NOERR) go to 10
   err = nf_get_vara_real(ncbathy,id,start,edges,wrk)
   if (err .ne. NF_NOERR) go to 10
   indx=1
   do i=il,ih
     lonc(iloc+i-il,:) = wrk(indx)
     indx=indx+1
   end do
   x = lonc(iloc+1,1) - lonc(iloc,1)
   lonc(iloc-1,:) = lonc(iloc,:) - x
   x = lonc(iloc+ih-il,1) - lonc(iloc+ih-il-1,1)
   lonc(iloc+ih-il+1,:) = lonc(iloc+ih-il,:) + x

   start(1) = il
   edges(1) = ih-il+2
   err = nf_inq_varid(ncbathy,"lonx",id)
   if (err .ne. NF_NOERR) then
      LEVEL4 'Can not read lonx - generating from lonc'
      do i=il-1,ih
         x = (lonc(iloc+i-il,1) + lonc(iloc+i-il+1,1))/2. 
         lonx(i-il+iloc,:) = x
      end do
   else
      err = nf_get_vara_real(ncbathy,id,start,edges,wrk)
      if (err .ne. NF_NOERR) go to 10
      indx=1
      do i=il-1,ih
        lonx(iloc+i-il,:) = wrk(indx)
        indx=indx+1
      end do
   end if
   x = lonc(iloc+ih-il+1,1)-lonc(iloc+ih-il,1) 
   lonx(iloc+ih-il+1,:) = lonc(iloc+ih-il+1,:) + x/2.

!  Reading latitude for C and X points
   jl = max(jmin+joff,1); jh = min(jmax+joff,jextr)
   start(1) = jl ; edges(1) = jh-jl+1
   err = nf_inq_varid(ncbathy,"lat",id)
   if (err .ne. NF_NOERR) go to 10
   err = nf_get_vara_real(ncbathy,id,start,edges,wrk)
   if (err .ne. NF_NOERR) go to 10
   indx=1
   do j=jl,jh
     latc(:,j-jl+jloc) = wrk(indx)
     indx=indx+1
   end do
   x = latc(1,jloc+1) - latc(1,jloc)
   latc(:,jloc-1) = latc(:,jloc) - x
   x = latc(1,jh-jl+jloc) - latc(1,jh-jl-1+jloc)
   latc(:,jh-jl+jloc+1) = latc(:,jh-jl+jloc) + x

   start(1) = jl
   edges(1) = jh+1
   err = nf_inq_varid(ncbathy,"latx",id)
   if (err .ne. NF_NOERR) then
      LEVEL4 'Can not read latx - generating from latc'
      do j=jl-1,jh
         x = (latc(1,jloc+j-jl) + latc(1,jloc+j-jl+1))/2. 
         latx(:,jloc+j-jl) = x
      end do
   else
      err = nf_get_vara_real(ncbathy,id,start,edges,wrk)
      if (err .ne. NF_NOERR) go to 10
      indx=1
      do j=jl-1,jh
        latx(:,j-jl+jloc) = wrk(indx)
        indx=indx+1
      end do
   endif
   x = latc(1,jloc+jh-jl+1)-latc(1,jloc+jh-jl) 
   latx(:,jloc+jh-jl+1) = latc(:,jloc+jh-jl+1) + x/2.
#endif

#ifdef CURVILINEAR
   LEVEL3 'reading curvi-linear grid information'
!  need to read xx,yx,xc,yc,xu,yu,xv,yv
!  need to read lonx,latx,lonc,latc,lonu,lonv,lonv,latv
!  also read dxdyc,dydxc,angle

!  0:imax,0:jmax
   il = imin+ioff-1; ih = imax+ioff
   jl = jmin+joff-1; jh = jmax+joff
   call read_2d_field("xx",imin,imax,jmin,jmax,xx,il,ih,jl,jh)
   call read_2d_field("yx",imin,imax,jmin,jmax,yx,il,ih,jl,jh)
!  1:imax,1:jmax
   il = imin+ioff;   ih = imax+ioff
   jl = jmin+joff;   jh = jmax+joff
   call read_2d_field("xc",imin,imax,jmin,jmax,xc,il,ih,jl,jh)
   call read_2d_field("yc",imin,imax,jmin,jmax,yc,il,ih,jl,jh)
!  0:imax,1:jmax
   il = imin+ioff-1; ih = imax+ioff
   jl = jmin+joff;   jh = jmax+joff
   call read_2d_field("xu",imin,imax,jmin,jmax,xu,il,ih,jl,jh)
   call read_2d_field("yu",imin,imax,jmin,jmax,yu,il,ih,jl,jh)
!  1:imax,0:jmax
   il = imin+ioff;   ih = imax+ioff
   jl = jmin+joff-1; jh = jmax+joff
   call read_2d_field("xv",imin,imax,jmin,jmax,xv,il,ih,jl,jh)
   call read_2d_field("yv",imin,imax,jmin,jmax,yv,il,ih,jl,jh)
!  0:imax,0:jmax
   il = imin+ioff-1; ih = imax+ioff
   jl = jmin+joff-1; jh = jmax+joff
   call read_2d_field("lonx",imin,imax,jmin,jmax,lonx,il,ih,jl,jh)
   call read_2d_field("latx",imin,imax,jmin,jmax,latx,il,ih,jl,jh)
!  1:imax,1:jmax
   il = imin+ioff;   ih = imax+ioff
   jl = jmin+joff;   jh = jmax+joff
   call read_2d_field("lonc",imin,imax,jmin,jmax,lonc,il,ih,jl,jh)
   call read_2d_field("latc",imin,imax,jmin,jmax,latc,il,ih,jl,jh)
!  0:imax,1:jmax
   il = imin+ioff-1; ih = imax+ioff
   jl = jmin+joff;   jh = jmax+joff
   call read_2d_field("lonu",imin,imax,jmin,jmax,lonu,il,ih,jl,jh)
   call read_2d_field("latu",imin,imax,jmin,jmax,latu,il,ih,jl,jh)
!  1:imax,0:jmax
   il = imin+ioff;   ih = imax+ioff
   jl = jmin+joff-1; jh = jmax+joff
   call read_2d_field("lonv",imin,imax,jmin,jmax,lonv,il,ih,jl,jh)
   call read_2d_field("latv",imin,imax,jmin,jmax,latv,il,ih,jl,jh)
!  1:imax,1:jmax
   il = imin+ioff;   ih = imax+ioff
   jl = jmin+joff;   jh = jmax+joff
   call read_2d_field("dxdyc",imin,imax,jmin,jmax,dxdyc,il,ih,jl,jh)
   call read_2d_field("dydxc",imin,imax,jmin,jmax,dydxc,il,ih,jl,jh)
   call read_2d_field("angle",imin,imax,jmin,jmax,angle,il,ih,jl,jh)

#ifdef DK_COARSE_CURV_TEST
!  All coordinates are here divided by 1.85185 since the curvilinear grid
!  of this test case has this factor as an error. The more accurate
!  correction suggested by Carsten Hansen needs to be further checked.
!  KBK 20020223 - one version divides by 1.07 one multiplies
!   xx=xx/(cos(latx*deg2rad)*1.07)
!   xc=xc/(cos(latc*deg2rad)*1.07)
!   xu=xu/(cos(latu*deg2rad)*1.07)
!   xv=xv/(cos(latv*deg2rad)*1.07)
!   yx=yx/(cos(latx*deg2rad)*1.07)
!   yc=yc/(cos(latc*deg2rad)*1.07)
!   yu=yu/(cos(latu*deg2rad)*1.07)
!   yv=yv/(cos(latv*deg2rad)*1.07)
   xx=xx/1.85185
   xc=xc/1.85185
   xu=xu/1.85185
   xv=xv/1.85185
   yx=yx/1.85185
   yc=yc/1.85185
   yu=yu/1.85185
   yv=yv/1.85185
#endif

   conv=-angle*rad2deg

#endif 

#ifdef DEBUG
   write(debug,*) 'Leaving get_bathymetry()'
   write(debug,*)
#endif
   return
10 FATAL 'get_bathymetry: ',nf_strerror(err)
   stop
   end subroutine get_bathymetry
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: read_2d_field -
!
! !INTERFACE:
   subroutine read_2d_field(name,imin,imax,jmin,jmax,field, &
                            il,ih,jl,jh,iloc,jloc)
!
! !USES:
   use ncdfin
   IMPLICIT NONE
!
! !DESCRIPTION:

! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: name
   integer, intent(in)                 :: imin,imax,jmin,jmax
   integer, intent(in)                 :: il,ih,jl,jh
   integer, intent(in)                 :: iloc,jloc
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: field(E2DFIELD)
!
! !REVISION HISTORY:
!
!  22Apr99   Karsten Bolding & Hans Burchard  Initial code.
!
! !LOCAL VARIABLES:
   integer                   :: id,size,i,j,indx
   integer                   :: err,rc
   REAL_4B, allocatable      :: wrk(:)
!EOP
!-------------------------------------------------------------------------
#include"netcdf.inc"

   size = (ih-il+1)*(jh-jl+1)
   allocate(wrk(size),stat=rc)
   if (rc /= 0) stop 'read_2d_field: Error allocating work-space'

   start(1) = il ; start(2) = jl
   edges(1) = ih-il+1 ; edges(2) = jh-jl+1

   err = nf_inq_varid(ncbathy,name,id)
   if (err .ne. NF_NOERR) go to 10

   err = nf_get_vara_real(ncbathy,id,start,edges,wrk)
   if (err .ne. NF_NOERR) go to 10
   indx = 1
   do j=jl,jh
     do i=il,ih
       field(i-il+iloc,j-jl+jloc) = wrk(indx)
       indx = indx+1
     end do
   end do

#ifdef FORTRAN90
   deallocate(wrk,stat=rc)
   if (rc /= 0) stop 'read_2d_field: Error de-allocating work-space'
#endif
   return
10 STDERR 'read_2d_field: error reading ',trim(name)
   STDERR 'this might cause problems later'
   return
   end subroutine read_2d_field
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
