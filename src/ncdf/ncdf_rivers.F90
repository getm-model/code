!$Id: ncdf_rivers.F90,v 1.9 2010-03-30 10:06:51 kb Exp $
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
   use time, only: string_to_julsecs,time_diff,add_secs
   use time, only: julianday,secondsofday,juln,secsn,timestep
   use time, only: write_time_string,timestr
   use rivers, only: nriver,river_data,river_name,river_flow,river_factor
   use rivers, only: ok,rriver,real_river_name,river_split
   use rivers, only: temp_missing,salt_missing
   use rivers, only: use_river_temp,use_river_salt,river_temp,river_salt
#ifdef GETM_BIO
   use bio, only: bio_calc
   use bio_var, only: numc,var_names
   use rivers, only: river_bio
#endif
   IMPLICIT NONE
!
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_river_input_ncdf,get_river_data_ncdf
!
! !PRIVATE DATA MEMBERS:
   REALTYPE                            :: offset
   integer                             :: ncid,ndims,dims(2),unlimdimid,textr
   integer                             :: start(1),edges(1)
   integer                             :: timedim,time_id
   integer, allocatable                :: r_ids(:)
   integer, allocatable                :: salt_id(:)
   integer, allocatable                :: temp_id(:)
   integer, allocatable                :: r_salt(:)
   integer, allocatable                :: r_temp(:)
   REAL_4B, allocatable                :: river_times(:)
#ifdef GETM_BIO
   integer, allocatable                :: bio_id(:,:)
   integer, allocatable                :: r_bio(:,:)
#endif
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ncdf_rivers.F90,v $
!  Revision 1.9  2010-03-30 10:06:51  kb
!  removed temporarely diagnostic output
!
!  Revision 1.8  2007-09-30 13:00:43  kbk
!  prints real time as part of progessoutput
!
!  Revision 1.7  2005-09-23 11:27:43  kbk
!  support fo nutrient loading in rivers
!
!  Revision 1.6  2005/05/04 11:45:29  kbk
!  adding model time stamp on IO
!
!  Revision 1.5  2005/01/13 09:20:47  kbk
!  support for T and S specifications in rivers - Stips
!
!  Revision 1.4  2004/04/06 16:32:29  kbk
!  TimeDiff --> time_diff
!
!  Revision 1.3  2003/10/14 10:05:54  kbk
!  checks if indices are in subdomain + cleaning
!
!  Revision 1.2  2003/04/23 11:54:03  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.1.1.1  2002/05/02 14:01:48  gotm
!  recovering after CVS crash
!
!  Revision 1.1  2001/10/07 14:50:22  bbh
!  Reading river data implemented - NetCFD
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
   character(len=*), intent(in)        :: fn
   integer, intent(in)                 :: nstart
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
   integer                   :: i,j,m,n
   integer                   :: err
   character(len=19)         :: tbuf
   integer                   :: j1,s1,j2,s2
   character(len=256)        :: time_units
   character(len=256)        :: bio_name
!EOP
!-------------------------------------------------------------------------
   include "netcdf.inc"
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_river_input_ncdf() # ',Ncall
#endif

   LEVEL3 'init_river_input_ncdf'

   allocate(r_ids(rriver),stat=err)
   if (err /= 0) stop 'ncdf_river: Error allocating memory (r_ids)'

   allocate(r_salt(rriver),stat=err)
   if (err /= 0) stop 'ncdf_river: Error allocating memory (r_salt)'
   allocate(r_temp(rriver),stat=err)
   if (err /= 0) stop 'ncdf_river: Error allocating memory (r_temp)'

   allocate(salt_id(rriver),stat=err)
   if (err /= 0) stop 'ncdf_river: Error allocating memory (salt_id)'
   allocate(temp_id(rriver),stat=err)
   if (err /= 0) stop 'ncdf_river: Error allocating memory (temp_id)'

#ifdef GETM_BIO
   allocate(r_bio(rriver,numc),stat=err)
   if (err /= 0) stop 'ncdf_river: Error allocating memory (r_bio)'
   allocate(bio_id(rriver,numc),stat=err)
   if (err /= 0) stop 'ncdf_river: Error allocating memory (bio_id)'
   bio_id = -1
#endif

   err = nf_open(fn,NCNOWRIT,ncid)
   if (err .ne. NF_NOERR) go to 10

   err = nf_inq_unlimdim(ncid,unlimdimid)
   if (err .NE. NF_NOERR) go to 10

   err = nf_inq_dimlen(ncid,unlimdimid,textr)
   if (err .ne. NF_NOERR) go to 10

   err = nf_inq_varid(ncid,"time",time_id)
   if (err .ne. NF_NOERR) go to 10

   do n=1,rriver
      err = nf_inq_varid(ncid,real_river_name(n),r_ids(n))
      if (err .ne. NF_NOERR) go to 10
      r_salt(n) = 0
      if ( use_river_salt ) then
         err =  nf_get_att_int(ncid,r_ids(n),'salt',r_salt(n))
!        if (err .ne. NF_NOERR) go to 10
         if (r_salt(n) .eq. 1 ) then
            LEVEL3 'river salinity:    ',trim(real_river_name(n))//trim('_salt')
            err =  nf_inq_varid(ncid,trim(real_river_name(n))//trim('_salt'),salt_id(n))
            if (err .ne. NF_NOERR) go to 10
         end if
      end if
      r_temp(n) = 0
      if ( use_river_temp ) then
         err =  nf_get_att_int(ncid,r_ids(n),'temp',r_temp(n))
!        if (err .ne. NF_NOERR) go to 10
         if (r_temp(n) .eq. 1 ) then
            LEVEL3 'river temperature: ',trim(real_river_name(n))//trim('_temp')
            err =  nf_inq_varid(ncid,trim(real_river_name(n))//trim('_temp'),temp_id(n))
            if (err .ne. NF_NOERR) go to 10
         end if
      end if

#ifdef GETM_BIO
      do m=1,numc
         bio_name=trim(real_river_name(n))//'_'//trim(var_names(m))
         err =  nf_inq_varid(ncid,trim(bio_name),bio_id(n,m))
         if (err .ne. NF_NOERR) then
            bio_id(n,m) = -1
         end if
         if ( bio_id(n,m) .ne. -1 ) then
            LEVEL4 trim(real_river_name(n)),': ',trim(var_names(m))
         end if
      end do
#endif
   end do

   allocate(river_times(textr),stat=err)
   if (err /= 0) stop  &
      'init_river_input_ncdf: Error allocating memory (river_times)'

   err =  nf_get_att_text(ncid,time_id,'units',time_units)
   if (err .ne. NF_NOERR) go to 10
   call string_to_julsecs(time_units,j1,s1)
   err = nf_get_var_real(ncid,time_id,river_times)
   if (err .ne. NF_NOERR) go to 10

   offset = time_diff(julianday,secondsofday,j1,s1)
   if( offset .lt. river_times(1) ) then
      FATAL 'Model simulation starts before available river data'
      call write_time_string(julianday,secondsofday,tbuf)
      FATAL 'Simulation starts: ',tbuf
      call add_secs(j1,s1,nint(river_times(1)),j2,s2)
      call write_time_string(j2,s2,tbuf)
      FATAL 'River file starts: ',tbuf
      stop 'init_river_input_ncdf'
   else
      LEVEL3 'River offset time ',offset
   endif

!  check if the bdy data file is long enough
   if( time_diff(juln,secsn,j1,s1) .gt. river_times(textr) ) then
      FATAL 'Not sufficient river data available'
      call write_time_string(juln,secsn,tbuf)
      FATAL 'Simulation ends: ',tbuf
      call add_secs(j1,s1,nint(river_times(textr)),j2,s2)
      call write_time_string(j2,s2,tbuf)
      FATAL 'River file ends: ',tbuf
      stop 'init_river_input_ncdf'
   endif

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
   integer, intent(in)                 :: loop
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
   integer                   :: i,j,n,nn,ni,m,indx,err
   REALTYPE                  :: t
   REAL_4B                   :: x(1)
   logical, save             :: first=.true.
   integer, save             :: save_n=1,last_indx=-1
   REALTYPE, save            :: t_1,t_2
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
      call write_time_string()
      LEVEL3 timestr, ': reading river data .... ',indx
      last_indx = indx
      start(1) = indx
      nn = 1
      ni = 1
      do n =1,nriver
         if (ni .le. nriver) then
            if (ok(ni) .ne. 0) then
               err = nf_get_vara_real(ncid,r_ids(nn),start,edges,x)
               if (err .ne. NF_NOERR) go to 10
               do m=1,river_split(ni)
                  river_flow(ni+m-1) = river_factor*x(1)
                  river_salt(ni+m-1) = salt_missing
                  river_temp(ni+m-1) = temp_missing
               end do
               if ( r_salt(nn) .eq. 1 ) then
                  err = nf_get_vara_real(ncid,salt_id(nn),start,edges,x)
                  if (err .ne. NF_NOERR) go to 10
                  do m=1,river_split(ni)
                     river_salt(ni+m-1) = x(1)
                  end do
               end if
               if ( r_temp(nn) .eq. 1 ) then
                  err = nf_get_vara_real(ncid,temp_id(nn),start,edges,x)
                  if (err .ne. NF_NOERR) go to 10
                  do m=1,river_split(ni)
                     river_temp(ni+m-1) = x(1)
                  end do
               end if
#ifdef GETM_BIO
               do j=1,numc
                  if (bio_id(nn,j) .gt. 0) then
                     err = nf_get_vara_real(ncid,bio_id(nn,j),start,edges,x)
                     if (err .ne. NF_NOERR) go to 10
                     do m=1,river_split(ni)
                        river_bio(ni+m-1,j) = x(1)
                     end do
                  end if
               end do
#endif
            end if
            nn = nn + 1 
            ni = ni + river_split(ni) 
         end if
      end do
   end if

#else
!AS this is not to be used !
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
         if (ok(n) .ne. 0) then
            err = nf_get_vara_real(ncid,r_ids(n),start,edges,x)
            if (err .ne. NF_NOERR) go to 10
            river_flow(n) = x(1)
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
