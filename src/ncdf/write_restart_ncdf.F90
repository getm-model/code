!$Id: write_restart_ncdf.F90,v 1.1 2007-09-21 13:03:42 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Writes variables to a GETM NetCDF hotstart file
!
! !INTERFACE:
   subroutine write_restart_ncdf(secs,loop,julianday,secondsofday)
!
! !DESCRIPTION:
!  Writes to a NetCDF file previously created using the 
!  create_restart_ncdf() subroutine all variables necessary to make
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
   use variables_3d
#ifdef GETM_BIO
   use bio_var, only: numc
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in) :: secs
   integer, intent(in)       :: loop,julianday,secondsofday
!
! !DEFINED PARAMTERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  $Log: write_restart_ncdf.F90,v $
!  Revision 1.1  2007-09-21 13:03:42  kbk
!  added drop-in NetCDF replacement for binary hotstart file (default is binary)
!
!
! !LOCAL VARIABLES:
   integer         :: k
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

   start(1) = 1; edges(1) = imax-imin+1
   start(2) = 1; edges(2) = jmax-jmin+1
   start(3) = 1; edges(3) = kmax+1

   select case (grid_type)
      case (1) ! cartesian
         status = nf90_put_var(ncid,xax_id,xc(1:imax,1))
         if (status .NE. NF90_NOERR) go to 10
         status = nf90_put_var(ncid,yax_id,yc(1,1:jmax))
         if (status .NE. NF90_NOERR) go to 10
      case (2) ! spherical
         status = nf90_put_var(ncid,xax_id,lonc(1:imax,1))
         if (status .NE. NF90_NOERR) go to 10
         status = nf90_put_var(ncid,yax_id,latc(1,1:jmax))
         if (status .NE. NF90_NOERR) go to 10
   end select

   status = nf90_put_var(ncid,zax_id,zax(0:kmax))
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,z_id,z(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,zo_id,zo(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,U_id,U(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,zu_id,zu(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,zub_id,zub(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,SlUx_id,SlUx(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,Slru_id,Slru(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,V_id,V(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,zv_id,zv(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,zvb_id,zvb(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,SlVx_id,SlVx(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,Slrv_id,Slrv(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

#ifndef NO_3D
   status = &
   nf90_put_var(ncid,ssen_id,ssen(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,ssun_id,ssun(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,ssvn_id,ssvn(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,sseo_id,sseo(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,ssuo_id,ssuo(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,ssvo_id,ssvo(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,Uinto_id,Uinto(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,Vinto_id,Vinto(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,uu_id,uu(imin:imax,jmin:jmax,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,vv_id,vv(imin:imax,jmin:jmax,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,ww_id,ww(imin:imax,jmin:jmax,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,uuEx_id,uuEx(imin:imax,jmin:jmax,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,vvEx_id,vvEx(imin:imax,jmin:jmax,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,tke_id,tke(imin:imax,jmin:jmax,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,eps_id,eps(imin:imax,jmin:jmax,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,num_id,num(imin:imax,jmin:jmax,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,nuh_id,nuh(imin:imax,jmin:jmax,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

#ifndef NO_BAROCLINIC
   status = &
   nf90_put_var(ncid,T_id,T(imin:imax,jmin:jmax,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,S_id,S(imin:imax,jmin:jmax,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10
#endif
#ifdef SPM
   status = &
   nf90_put_var(ncid,spm_id,spm(imin:imax,jmin:jmax,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10

   status = &
   nf90_put_var(ncid,spmpool_id,spmpool(imin:imax,jmin:jmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10
#endif
#ifdef GETM_BIO
   start(4) = 1; edges(4) = numc
   status = &
   nf90_put_var(ncid,bio_id,cc3d(imin:imax,jmin:jmax,0:kmax,numc), &
                start,edges)
   if (status .NE. NF90_NOERR) go to 10
#endif
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
