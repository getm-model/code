!$Id: ncdf_rivers.F90,v 1.1 2002-05-02 14:01:48 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdf_river -
!
! !INTERFACE:
   module ncdf_river
!
! !DESCRIPTION:
!
! !USES:
   use time, only: string_to_julsecs,TimeDiff,add_secs,in_interval
   use time, only: jul0,secs0,julianday,secondsofday,timestep
   use rivers, only: nriver,river_data,river_name,river_flux,river_factor,ok
   IMPLICIT NONE
!
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_river_input_ncdf,get_river_data_ncdf
!
! !PRIVATE DATA MEMBERS:
   REALTYPE	:: offset
   integer 	:: ncid,ndims,dims(2),unlimdimid,textr
   integer 	:: start(1),edges(1)
   integer	:: timedim,time_id
   integer, allocatable	:: r_ids(:)
   REAL_4B, allocatable	:: river_times(:)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ncdf_rivers.F90,v $
!  Revision 1.1  2002-05-02 14:01:48  gotm
!  Initial revision
!
!  Revision 1.1  2001/10/07 14:50:22  bbh
!  Reading river data implemented - NetCFD
!
!
! !TO DO:
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_river_input_ncdf -
!
! !INTERFACE:
   subroutine init_river_input_ncdf(fn,nstart)
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)	:: fn
   integer, intent(in)		:: nstart
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   integer      :: i,j,n
   integer      :: err
   integer      :: j1,s1,j2,s2
   character(len=256)	:: time_units
!EOP
!-------------------------------------------------------------------------
   include "netcdf.inc"
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_river_input_ncdf() # ',Ncall
#endif

   LEVEL3 'init_river_input_ncdf'

   allocate(r_ids(nriver),stat=err)
   if (err /= 0) stop 'ncdf_river: Error allocating memory (r_ids)'

   err = nf_open(fn,NCNOWRIT,ncid)
   if (err .ne. NF_NOERR) go to 10

   err = nf_inq_unlimdim(ncid,unlimdimid)
   if (err .NE. NF_NOERR) go to 10

   err = nf_inq_dimlen(ncid,unlimdimid,textr)
   if (err .ne. NF_NOERR) go to 10

   err = nf_inq_varid(ncid,"time",time_id)
   if (err .ne. NF_NOERR) go to 10

   do n=1,nriver
      err = nf_inq_varid(ncid,river_name(n),r_ids(n))
      if (err .ne. NF_NOERR) go to 10
   end do

   allocate(river_times(textr),stat=err)
   if (err /= 0) stop  &
      'init_river_input_ncdf: Error allocating memory (river_times)'

   err =  nf_get_att_text(ncid,time_id,'units',time_units)
   if (err .ne. NF_NOERR) go to 10
   call string_to_julsecs(time_units,j1,s1)
   err = nf_get_var_real(ncid,time_id,river_times)
   if (err .ne. NF_NOERR) go to 10

   offset = TimeDiff(jul0,secs0,j1,s1)
   if( offset .lt. _ZERO_ ) then    !HB Karsten check, I changed gt to lt
      FATAL 'Model simulation starts before available river data'
      stop 'init_river_input_ncdf'
   endif

   call add_secs(j1,s1,nint(river_times(textr)),j2,s2)
!kbkSTDERR TimeDiff(j1,s1,j2,s2)
!   if( TimeDiff(j1,s1,j2,s2) .lt. _ZERO_ ) then
!      FATAL 'Not sufficient river data available'
!      stop 'init_river_input_ncdf'
!   endif

   call get_river_data_ncdf(nstart)

#ifdef DEBUG
   write(debug,*) 'Leaving init_river_input_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'init_river_input_ncdf: ',nf_strerror(err)
   stop
   end subroutine init_river_input_ncdf
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_river_data_ncdf - .
!
! !INTERFACE:
   subroutine get_river_data_ncdf(loop)
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   integer, intent(in)	:: loop
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   integer      :: i,n,indx,err
   REALTYPE	:: t
   REAL_4B	:: x(1)
   logical, save 	:: first=.true.
   integer, save 	:: save_n=1,last_indx=-1
   REALTYPE, save	:: t_1,t_2
!EOP
!-------------------------------------------------------------------------
   include "netcdf.inc"
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'get_river_data_ncdf() # ',Ncall
#endif
   edges(1) = 1
#define NO_INTERPOL
#ifdef NO_INTERPOL
   t = (loop-1)*timestep
   do indx=save_n,textr
      if (river_times(indx) .ge. real(t + offset)) EXIT
   end do
   if (indx .gt. last_indx) then
      LEVEL3 'reading river data - indx = ',indx
      last_indx = indx
      start(1) = indx
      do n =1,nriver
         if (ok(n) .gt. 0) then
            err = nf_get_vara_real(ncid,r_ids(n),start,edges,x)
            if (err .ne. NF_NOERR) go to 10
            river_flux(n) = river_factor*x(1)
         end if
      end do
   end if
#else
   t = loop*timestep
   do indx=save_n,textr
      if (river_times(indx) .gt. real(t + offset)) EXIT
   end do
   ! First time through we have to initialize t_1
   if (first) then
      LEVEL3 'reading first river data - indx = ',indx
      first = .false.
      if (indx .gt. 1) then
         indx = indx-1
      end if
      save_n = indx
      start(1) = indx
      t_1 = river_times(indx) - offset
      t_2 = t_1

      do n =1,nriver
         if (ok(n) .gt. 0) then
            err = nf_get_vara_real(ncid,r_ids(n),start,edges,x)
            if (err .ne. NF_NOERR) go to 10
            river_flux(n) = x(1)
            STDERR x(1),river_flux(n)
         end if
      end do
   else
      if (indx .gt. save_n) then
         LEVEL3 'reading new river data - indx = ',indx
         save_n = indx
         t_1 = t_2
         t_2 = river_times(indx) - offset
      end if
   end if
#endif
#ifdef DEBUG
   write(debug,*) 'Leaving get_river_data_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'get_river_data_ncdf: ',nf_strerror(err)
   stop
   end subroutine get_river_data_ncdf
!EOC

!-----------------------------------------------------------------------

   end module ncdf_river

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
