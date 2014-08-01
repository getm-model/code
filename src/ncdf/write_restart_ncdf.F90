#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Writes variables to a GETM NetCDF hotstart file
!
! !INTERFACE:
   subroutine write_restart_ncdf(runtype,secs,loop,julianday,secondsofday)
!
! !DESCRIPTION:
!  Writes to a NetCDF file previously created using the
!  create\_restart\_ncdf() subroutine all variables necessary to make
!  a correct GETM hotstart. The Fortran variables are written directly
!  into the corresponding NetCDF variable.
!
! !USES:
   use netcdf
   use ncdf_restart
   use domain, only: grid_type
   use domain, only: xc,yc,lonc,latc
   use domain, only: imin,imax,jmin,jmax,kmax
   use variables_2d
#ifndef NO_3D
   use variables_3d
#ifdef GETM_BIO
   use bio, only: bio_calc
   use bio_var, only: numc
#endif
#ifdef _FABM_
   use getm_fabm, only: fabm_pel,fabm_ben
#endif
#endif
#ifdef SPM
   use suspended_matter
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)       :: runtype
   REALTYPE, intent(in)      :: secs ! not used now
   integer, intent(in)       :: loop,julianday,secondsofday
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL VARIABLES:
   integer         :: k,n, rc
   REALTYPE, allocatable :: zax(:)
!EOP
!-----------------------------------------------------------------------
!BOC
   status = nf90_put_var(ncid,loop_id,loop)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_put_var(ncid,julianday_id,julianday)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_put_var(ncid,secondsofday_id,secondsofday)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_put_var(ncid,timestep_id,timestep)
   if (status .NE. NF90_NOERR) go to 10

   allocate(zax(0:kmax),stat=status)
   if (rc /= 0) &
      stop 'write_restart_ncdf(): Error allocating memory (zax)'
   do k=0,kmax
      zax(k)=k
   end do

#ifdef _WRITE_HOT_HALOS_
   LEVEL3 'include HALOs in NetCDF hotstart files'
   start(1) = 1; edges(1) = (imax+HALO)-(imin-HALO)+1
   start(2) = 1; edges(2) = (jmax+HALO)-(jmin-HALO)+1
#else
   start(1) = 1; edges(1) = imax-imin+1
   start(2) = 1; edges(2) = jmax-jmin+1
#endif
   start(3) = 1; edges(3) = kmax+1

   select case (grid_type)
      case (1) ! cartesian
         status = nf90_put_var(ncid,xax_id,xc(_IRANGE_,1))
         if (status .NE. NF90_NOERR) go to 10
         status = nf90_put_var(ncid,yax_id,yc(1,_JRANGE_))
         if (status .NE. NF90_NOERR) go to 10
      case (2) ! spherical
         status = nf90_put_var(ncid,xax_id,lonc(_IRANGE_,1))
         if (status .NE. NF90_NOERR) go to 10
         status = nf90_put_var(ncid,yax_id,latc(1,_JRANGE_))
         if (status .NE. NF90_NOERR) go to 10
   end select

   status = nf90_put_var(ncid,zax_id,zax(0:kmax))
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,z_id,z(_2D_W_HOT_),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,zo_id,zo(_2D_W_HOT_),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,U_id,U(_2D_W_HOT_),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,SlUx_id,SlUx(_2D_W_HOT_),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,Slru_id,Slru(_2D_W_HOT_),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,V_id,V(_2D_W_HOT_),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,SlVx_id,SlVx(_2D_W_HOT_),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,Slrv_id,Slrv(_2D_W_HOT_),start,edges)
   if (status .NE. NF90_NOERR) go to 10

#ifndef NO_3D
   if (runtype .ge. 2)  then
      status = &
      nf90_put_var(ncid,ssen_id,ssen(_2D_W_HOT_),start,edges)
      if (status .NE. NF90_NOERR) go to 10

      status = &
      nf90_put_var(ncid,ssun_id,ssun(_2D_W_HOT_),start,edges)
      if (status .NE. NF90_NOERR) go to 10

      status = &
      nf90_put_var(ncid,ssvn_id,ssvn(_2D_W_HOT_),start,edges)
      if (status .NE. NF90_NOERR) go to 10

      status = &
      nf90_put_var(ncid,sseo_id,sseo(_2D_W_HOT_),start,edges)
      if (status .NE. NF90_NOERR) go to 10

      status = &
      nf90_put_var(ncid,ssuo_id,ssuo(_2D_W_HOT_),start,edges)
      if (status .NE. NF90_NOERR) go to 10

      status = &
      nf90_put_var(ncid,ssvo_id,ssvo(_2D_W_HOT_),start,edges)
      if (status .NE. NF90_NOERR) go to 10

      status = &
      nf90_put_var(ncid,Uinto_id,Uinto(_2D_W_HOT_),start,edges)
      if (status .NE. NF90_NOERR) go to 10

      status = &
      nf90_put_var(ncid,Vinto_id,Vinto(_2D_W_HOT_),start,edges)
      if (status .NE. NF90_NOERR) go to 10

      status = &
      nf90_put_var(ncid,uu_id,uu(_3D_W_HOT_),start,edges)
      if (status .NE. NF90_NOERR) go to 10

      status = &
      nf90_put_var(ncid,vv_id,vv(_3D_W_HOT_),start,edges)
      if (status .NE. NF90_NOERR) go to 10

      status = &
      nf90_put_var(ncid,ww_id,ww(_3D_W_HOT_),start,edges)
      if (status .NE. NF90_NOERR) go to 10

      status = &
      nf90_put_var(ncid,uuEx_id,uuEx(_3D_W_HOT_),start,edges)
      if (status .NE. NF90_NOERR) go to 10

      status = &
      nf90_put_var(ncid,vvEx_id,vvEx(_3D_W_HOT_),start,edges)
      if (status .NE. NF90_NOERR) go to 10

      status = &
      nf90_put_var(ncid,tke_id,tke(_3D_W_HOT_),start,edges)
      if (status .NE. NF90_NOERR) go to 10

      status = &
      nf90_put_var(ncid,eps_id,eps(_3D_W_HOT_),start,edges)
      if (status .NE. NF90_NOERR) go to 10

      status = &
      nf90_put_var(ncid,num_id,num(_3D_W_HOT_),start,edges)
      if (status .NE. NF90_NOERR) go to 10

      status = &
      nf90_put_var(ncid,nuh_id,nuh(_3D_W_HOT_),start,edges)
      if (status .NE. NF90_NOERR) go to 10

      status = &
      nf90_put_var(ncid,ho_id,ho(_3D_W_HOT_),start,edges)
      if (status .NE. NF90_NOERR) go to 10

      status = &
      nf90_put_var(ncid,hn_id,hn(_3D_W_HOT_),start,edges)
      if (status .NE. NF90_NOERR) go to 10

#ifndef NO_BAROCLINIC
      if (runtype .ge. 3) then
         status = &
         nf90_put_var(ncid,T_id,T(_3D_W_HOT_),start,edges)
         if (status .NE. NF90_NOERR) go to 10

         status = &
         nf90_put_var(ncid,S_id,S(_3D_W_HOT_),start,edges)
         if (status .NE. NF90_NOERR) go to 10
      end if
#endif
#ifdef SPM
      if (spm_calc) then
         status = &
         nf90_put_var(ncid,spm_id,spm(_3D_W_HOT_),start,edges)
         if (status .NE. NF90_NOERR) go to 10

         status = &
         nf90_put_var(ncid,spmpool_id,spm_pool(_3D_W_HOT_)
         if (status .NE. NF90_NOERR) go to 10
      end if
#endif
#ifdef _FABM_
      if (allocated(fabm_pel)) then
         start(4) = 1; edges(4) = size(fabm_pel,4)
         status = &
         nf90_put_var(ncid,fabm_pel_id,fabm_pel(_3D_W_HOT_,:),start,edges)
         if  (status .NE. NF90_NOERR) go to 10

         start(3) = 1; edges(3) = size(fabm_ben,3)
         if (edges(3).gt.0) then
            status = &
            nf90_put_var(ncid,fabm_ben_id,fabm_ben(_2D_W_HOT_,:),start,edges)
            if  (status .NE. NF90_NOERR) go to 10
         end if
      end if
#endif
#ifdef GETM_BIO
      if (bio_calc) then

         start(1) = 1; edges(1) = numc
#ifdef _WRITE_HOT_HALOS_
         start(2) = 1; edges(2) = (imax+HALO)-(imin-HALO)+1
         start(3) = 1; edges(3) = (jmax+HALO)-(jmin-HALO)+1
#else
         start(2) = 1; edges(2) = imax-imin+1
         start(3) = 1; edges(3) = jmax-jmin+1
#endif
         start(4) = 1; edges(4) = kmax+1
#if 0
         status = &
         nf90_put_var(ncid,bio_id,cc3d(1:numc,_3D_W_HOT_), &
                      start,edges)
         if (status .NE. NF90_NOERR) go to 10
#else
         do n=1,numc
            start(1) = n; edges(1) = 1
            status = &
            nf90_put_var(ncid,bio_id,cc3d(n,_3D_W_HOT_), &
                         start,edges)
            if (status .NE. NF90_NOERR) go to 10
         end do
#endif
      end if
#endif
   end if
#endif

   status = nf90_sync(ncid)
   if (status .NE. NF90_NOERR) go to 10

   return

10 FATAL 'write_restart_ncdf: ',nf90_strerror(status)
   stop
   return

   end subroutine write_restart_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Karsten Bolding (BBH)                           !
!-----------------------------------------------------------------------
