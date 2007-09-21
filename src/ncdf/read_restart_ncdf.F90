!$Id: read_restart_ncdf.F90,v 1.1 2007-09-21 13:03:42 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Read variables from a GETM NetCDF hotstart file
!
! !INTERFACE:
   subroutine read_restart_ncdf(loop,julianday,secondsofday,tstep)
!
! !DESCRIPTION:
!  Reads from a NetCDF files (with handler ncid) opened with
!  open_restart_ncdf(). All variable id's are initialised. The variables
!  can be read from hotstart files with the same dimensions as given by
!  imin:imax,jmin:jmax - or - from a hotstart file with the same dimensions
!  as topo.nc (and on the same grid). This allows to use 'ncmerge' to
!  combine a number of hotstart files in to one - make a new sub-domain
!  decomposition and use the newly created hotstart file. It might be 
!  necessary to use 'ncks' to cut the file to be have the same dimensions 
!  as topo.nc. Allowing for the file naming scheme in GETM links for each 
!  sub-domain should be made - e.g. ln -s restart.in restart.000.in; ln -s 
!  restart.in restart.001.in etc.\newline
!  Halo-zones are updated using calls to update_2d_halo() and 
!  update_3d_halo().
!
! !USES:
   use netcdf
   use ncdf_restart
   use domain, only: iextr,jextr,ioff,joff
   use domain, only: az,au,av
   use halo_zones, only: update_2d_halo,update_3d_halo,wait_halo
   use halo_zones, only: H_TAG,U_TAG,V_TAG
   use variables_2d
   use variables_3d
#ifdef GETM_BIO
   use bio_var, only: numc
#endif
   IMPLICIT NONE
!
! !OUTPUT PARAMETERS:
   integer, intent(out)      :: loop,julianday,secondsofday
   REALTYPE, intent(out)     :: tstep
!
! !DEFINED PARAMTERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  $Log: read_restart_ncdf.F90,v $
!  Revision 1.1  2007-09-21 13:03:42  kbk
!  added drop-in NetCDF replacement for binary hotstart file (default is binary)
!
!
! !LOCAL VARIABLES:
   integer         :: il,ih,iloc,ilen
   integer         :: jl,jh,jloc,jlen
!EOP
!-----------------------------------------------------------------------
!BOC
   status = nf90_get_var(ncid,loop_id,loop)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_get_var(ncid,julianday_id,julianday)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_get_var(ncid,secondsofday_id,secondsofday)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_get_var(ncid,timestep_id,tstep)
   if (status .NE. NF90_NOERR) go to 10

!  allows reading from subdomain or topo.nc sized files
!  i.e. ncmerged files cut to the same size as topo.nc
   if (xlen .eq. (imax-imin+1) .and. ylen .eq. (jmax-jmin+1) ) then
      il = 1 ; ih = xlen
      jl = 1 ; jh = ylen
      iloc = 1 ; jloc = 1
   else
      il = max(imin+ioff,1); ih = min(imax+ioff,iextr)
      jl = max(jmin+joff,1); jh = min(jmax+joff,jextr)
      iloc = max(imin-ioff,1); jloc = max(jmin-joff,1)
   end if
   ilen = ih-il+1 
   jlen = jh-jl+1
   start(1) = il ; edges(1) = ih-il+1
   start(2) = jl ; edges(2) = jh-jl+1
   start(3) = 1  ; edges(3) = kmax+1

   status = &
   nf90_get_var(ncid,z_id,z(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(z,z,az,imin,jmin,imax,jmax,H_TAG)
   call wait_halo(H_TAG)

   status = &
   nf90_get_var(ncid,zo_id,zo(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(zo,zo,az,imin,jmin,imax,jmax,H_TAG)
   call wait_halo(H_TAG)

   status = &
   nf90_get_var(ncid,U_id,U(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(U,U,au,imin,jmin,imax,jmax,U_TAG)
   call wait_halo(U_TAG)

   status = &
   nf90_get_var(ncid,zu_id,zu(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(zu,zu,au,imin,jmin,imax,jmax,U_TAG)
   call wait_halo(U_TAG)

   status = &
   nf90_get_var(ncid,zub_id,zub(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(zub,zub,au,imin,jmin,imax,jmax,U_TAG)
   call wait_halo(U_TAG)

   status = &
   nf90_get_var(ncid,SlUx_id,SlUx(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(SlUx,SlUx,au,imin,jmin,imax,jmax,U_TAG)
   call wait_halo(U_TAG)

   status = &
   nf90_get_var(ncid,Slru_id,Slru(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(Slru,Slru,au,imin,jmin,imax,jmax,U_TAG)
   call wait_halo(U_TAG)

   status = &
   nf90_get_var(ncid,V_id,V(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(V,V,av,imin,jmin,imax,jmax,V_TAG)
   call wait_halo(V_TAG)

   status = &
   nf90_get_var(ncid,zv_id,zv(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(zv,zv,av,imin,jmin,imax,jmax,V_TAG)
   call wait_halo(V_TAG)

   status = &
   nf90_get_var(ncid,zvb_id,zvb(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(zvb,zvb,av,imin,jmin,imax,jmax,V_TAG)
   call wait_halo(V_TAG)

   status = &
   nf90_get_var(ncid,SlVx_id,SlVx(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(SlVx,SlVx,av,imin,jmin,imax,jmax,V_TAG)
   call wait_halo(V_TAG)

   status = &
   nf90_get_var(ncid,Slrv_id,Slrv(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(Slrv,Slrv,av,imin,jmin,imax,jmax,V_TAG)
   call wait_halo(V_TAG)

#ifndef NO_3D
   status = &
   nf90_get_var(ncid,ssen_id,ssen(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(ssen,ssen,az,imin,jmin,imax,jmax,H_TAG)
   call wait_halo(H_TAG)

   status = &
   nf90_get_var(ncid,ssun_id,ssun(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(ssun,ssun,au,imin,jmin,imax,jmax,U_TAG)
   call wait_halo(U_TAG)

   status = &
   nf90_get_var(ncid,ssvn_id,ssvn(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(ssvn,ssvn,av,imin,jmin,imax,jmax,V_TAG)
   call wait_halo(V_TAG)

   status = &
   nf90_get_var(ncid,sseo_id,sseo(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(sseo,sseo,az,imin,jmin,imax,jmax,H_TAG)
   call wait_halo(H_TAG)

   status = &
   nf90_get_var(ncid,ssuo_id,ssuo(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(ssuo,ssuo,au,imin,jmin,imax,jmax,U_TAG)
   call wait_halo(U_TAG)

   status = &
   nf90_get_var(ncid,ssvo_id,ssvo(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(ssvo,ssvo,av,imin,jmin,imax,jmax,V_TAG)
   call wait_halo(V_TAG)

   status = &
   nf90_get_var(ncid,Uinto_id,Uinto(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(Uinto,Uinto,au,imin,jmin,imax,jmax,U_TAG)
   call wait_halo(U_TAG)

   status = &
   nf90_get_var(ncid,Vinto_id,Vinto(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(Vinto,Vinto,av,imin,jmin,imax,jmax,V_TAG)
   call wait_halo(V_TAG)

   status = &
   nf90_get_var(ncid,uu_id,uu(iloc:ilen,jloc:jlen,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_3d_halo(uu,uu,au,imin,jmin,imax,jmax,kmax,U_TAG)
   call wait_halo(U_TAG)

   status = &
   nf90_get_var(ncid,vv_id,vv(iloc:ilen,jloc:jlen,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_3d_halo(vv,vv,av,imin,jmin,imax,jmax,kmax,V_TAG)
   call wait_halo(V_TAG)

   status = &
   nf90_get_var(ncid,ww_id,ww(iloc:ilen,jloc:jlen,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_3d_halo(ww,ww,az,imin,jmin,imax,jmax,kmax,H_TAG)
   call wait_halo(H_TAG)

   status = &
   nf90_get_var(ncid,uuEx_id,uuEx(iloc:ilen,jloc:jlen,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_3d_halo(uuEx,uuEx,au,imin,jmin,imax,jmax,kmax,U_TAG)
   call wait_halo(U_TAG)

   status = &
   nf90_get_var(ncid,vvEx_id,vvEx(iloc:ilen,jloc:jlen,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_3d_halo(vvEx,vvEx,av,imin,jmin,imax,jmax,kmax,V_TAG)
   call wait_halo(V_TAG)

   status = &
   nf90_get_var(ncid,tke_id,tke(iloc:ilen,jloc:jlen,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_3d_halo(tke,tke,az,imin,jmin,imax,jmax,kmax,H_TAG)
   call wait_halo(H_TAG)

   status = &
   nf90_get_var(ncid,eps_id,eps(iloc:ilen,jloc:jlen,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_3d_halo(eps,eps,az,imin,jmin,imax,jmax,kmax,H_TAG)
   call wait_halo(H_TAG)

   status = &
   nf90_get_var(ncid,num_id,num(iloc:ilen,jloc:jlen,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_3d_halo(num,num,az,imin,jmin,imax,jmax,kmax,H_TAG)
   call wait_halo(H_TAG)

   status = &
   nf90_get_var(ncid,nuh_id,nuh(iloc:ilen,jloc:jlen,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_3d_halo(nuh,nuh,az,imin,jmin,imax,jmax,kmax,H_TAG)
   call wait_halo(H_TAG)

#ifndef NO_BAROCLINIC
   status = &
   nf90_get_var(ncid,T_id,T(iloc:ilen,jloc:jlen,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_3d_halo(T,T,az,imin,jmin,imax,jmax,kmax,H_TAG)
   call wait_halo(H_TAG)

   status = &
   nf90_get_var(ncid,S_id,S(iloc:ilen,jloc:jlen,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_3d_halo(S,S,az,imin,jmin,imax,jmax,kmax,H_TAG)
   call wait_halo(H_TAG)
#endif
#ifdef SPM
   status = &
   nf90_get_var(ncid,spm_id,spm(iloc:ilen,jloc:jlen,0:kmax),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_3d_halo(spm,spm,az,imin,jmin,imax,jmax,kmax,H_TAG)
   call wait_halo(H_TAG)

   status = &
   nf90_get_var(ncid,spmpool_id,spmpool(iloc:ilen,jloc:jlen),start,edges)
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(spmpool,spmpool,az,imin,jmin,imax,jmax,H_TAG)
   call wait_halo(H_TAG)

#endif
#ifdef GETM_BIO
   start(4) = 1; edges(4) = numc
   status = &
   nf90_get_var(ncid,bio_id,cc3d(iloc:ilen,jloc:jlen,0:kmax,numc), &
                start,edges)
   if (status .NE. NF90_NOERR) go to 10
#endif
#endif

   status = nf90_close(ncid)
   if (status .NE. NF90_NOERR) go to 10

   return

10 FATAL 'read_restart_ncdf: ',nf90_strerror(status)
   stop
   return

   end subroutine read_restart_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Karsten Bolding (BBH)                           !
!-----------------------------------------------------------------------
