#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Read variables from a GETM NetCDF hotstart file
!
! !INTERFACE:
   subroutine read_restart_ncdf(runtype,loop,julianday,secondsofday,tstep)
!
! !DESCRIPTION:
!  Reads from a NetCDF files (with handler ncid) opened with
!  open\_restart\_ncdf(). All variable id's are initialised. The variables
!  can be read from hotstart files with the same dimensions as given by
!  imin:imax,jmin:jmax - or - from a hotstart file with the same dimensions
!  as topo.nc (and on the same grid). This allows to use 'ncmerge' to
!  combine a number of hotstart files in to one - make a new sub-domain
!  decomposition and use the newly created hotstart file. It might be
!  necessary to use 'ncks' to cut the file to be have the same dimensions
!  as topo.nc. Allowing for the file naming scheme in GETM links for each
!  sub-domain should be made - e.g. ln -s restart.in restart.000.in; ln -s
!  restart.in restart.001.in etc.\newline
!  Halo-zones are updated using calls to update\_2d\_halo() and
!  update\_3d\_halo().
!
! !USES:
   use netcdf
   use ncdf_restart
   use domain, only: iextr,jextr,ioff,joff
   use domain, only: az,au,av
   use halo_zones, only: update_2d_halo,update_3d_halo,wait_halo
   use halo_zones, only: H_TAG,U_TAG,V_TAG
   use variables_2d
   use exceptions, only: getm_error
#ifndef NO_3D
   use variables_3d
#ifdef GETM_BIO
   use bio, only: bio_calc
   use bio_var, only: numc
   use getm_bio, only: bio_init_method
#endif
#ifdef _FABM_
   use getm_fabm, only: fabm_init_method
   use getm_fabm, only: fabm_pel,fabm_ben
#endif
#endif
#ifdef SPM
   use suspended_matter
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)      :: runtype
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
! !LOCAL VARIABLES:
   integer         :: il,ih,iloc,ilen,i,istart,istop
   integer         :: jl,jh,jloc,jlen,j,jstart,jstop
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
#ifdef __READ_HOT_HALOS_
   if (xlen .eq. ((imax+HALO)-(imin-HALO)+1) .and.  &
       ylen .eq. ((jmax+HALO)-(jmin-HALO)+1) ) then
      LEVEL3 'hotstart file(s) include HALO-zones'
      il = 1 ; ih = xlen - 2*HALO
      jl = 1 ; jh = ylen - 2*HALO
      start(1) = il + HALO
      start(2) = jl + HALO
#else
   if (xlen .eq. (imax-imin+1) .and. ylen .eq. (jmax-jmin+1) ) then
      LEVEL3 'hotstart file(s) do NOT include HALO-zones'
      il = 1 ; ih = xlen
      jl = 1 ; jh = ylen
      start(1) = il
      start(2) = jl
#endif
      iloc = 1 ; jloc = 1
      istart=1 ; jstart=1
   else
      if ( xlen .ne. iextr .or. ylen .ne. jextr ) then
         FATAL 'x- or y-dimemsions of hotstart different from iextr or jextr'
         call getm_error('read_restart_ncdf()','x- or y-dimemsions of hotstart different from iextr or jextr')
         stop 'read_restart_ncdf()'
      end if
      LEVEL3 'hotstart file has topo.nc size'
      ! Range in bathy index:
      il = max(imin+ioff,1); ih = min(imax+ioff,iextr)
      jl = max(jmin+joff,1); jh = min(jmax+joff,jextr)
      start(1) = il
      start(2) = jl
      ! Range in local sub-domain index:
      iloc = max(imin-ioff,1); jloc = max(jmin-joff,1)
      istart=max(imin-ioff,1); jstart=max(jmin-joff,1)
   end if
   ilen = ih-il+1 ; istop=istart+ilen-1
   jlen = jh-jl+1 ; jstop=jstart+jlen-1
   edges(1) = ih-il+1
   edges(2) = jh-jl+1
#ifndef NO_3D
   start(3) = 1  ; edges(3) = kmax+1
#endif

!  z is required
   status = &
   nf90_get_var(ncid,z_id,z(istart:istop,jstart:jstop),start(1:2),edges(1:2))
   LEVEL4 ' mpi-read returned with status=',status
   if (status .NE. NF90_NOERR) go to 10
   call update_2d_halo(z,z,az,imin,jmin,imax,jmax,H_TAG)
   call wait_halo(H_TAG)

   status = &
   nf90_get_var(ncid,zo_id,zo(istart:istop,jstart:jstop),start(1:2),edges(1:2))
   if (status .NE. NF90_NOERR) then
      LEVEL3 "read_restart_ncdf(): setting zo=z"
      zo=z
   else
      call update_2d_halo(zo,zo,az,imin,jmin,imax,jmax,H_TAG)
      call wait_halo(H_TAG)
   end if

   status = &
   nf90_get_var(ncid,U_id,U(istart:istop,jstart:jstop),start(1:2),edges(1:2))
   if (status .NE. NF90_NOERR) then
      LEVEL3 "read_restart_ncdf(): setting U=0"
      U=_ZERO_
   else
      call update_2d_halo(U,U,au,imin,jmin,imax,jmax,U_TAG)
      call wait_halo(U_TAG)
   end if

   status = &
   nf90_get_var(ncid,SlUx_id,SlUx(istart:istop,jstart:jstop),start(1:2),edges(1:2))
   if (status .NE. NF90_NOERR) then
      LEVEL3 "read_restart_ncdf(): setting SlUx=0"
      SlUx=_ZERO_
   else
      call update_2d_halo(SlUx,SlUx,au,imin,jmin,imax,jmax,U_TAG)
      call wait_halo(U_TAG)
   end if

   status = &
   nf90_get_var(ncid,Slru_id,Slru(istart:istop,jstart:jstop),start(1:2),edges(1:2))
   if (status .NE. NF90_NOERR) then
      LEVEL3 "read_restart_ncdf(): setting Slru=0"
   else
      call update_2d_halo(Slru,Slru,au,imin,jmin,imax,jmax,U_TAG)
      call wait_halo(U_TAG)
   end if

   status = &
   nf90_get_var(ncid,V_id,V(istart:istop,jstart:jstop),start(1:2),edges(1:2))
   if (status .NE. NF90_NOERR) then
      LEVEL3 "read_restart_ncdf(): setting V=0"
      V=_ZERO_
   else
      call update_2d_halo(V,V,av,imin,jmin,imax,jmax,V_TAG)
      call wait_halo(V_TAG)
   end if

   status = &
   nf90_get_var(ncid,SlVx_id,SlVx(istart:istop,jstart:jstop),start(1:2),edges(1:2))
   if (status .NE. NF90_NOERR) then
      LEVEL3 "read_restart_ncdf(): setting SlVx=0"
      SlVx=_ZERO_
   else
      call update_2d_halo(SlVx,SlVx,av,imin,jmin,imax,jmax,V_TAG)
      call wait_halo(V_TAG)
   end if

   status = &
   nf90_get_var(ncid,Slrv_id,Slrv(istart:istop,jstart:jstop),start(1:2),edges(1:2))
   if (status .NE. NF90_NOERR) then
      LEVEL3 "read_restart_ncdf(): setting Slrv=0"
      Slrv=_ZERO_
   else
      call update_2d_halo(Slrv,Slrv,av,imin,jmin,imax,jmax,V_TAG)
      call wait_halo(V_TAG)
   end if

#ifndef NO_3D
   if (runtype .ge. 2)  then
      status = &
      nf90_get_var(ncid,ssen_id,ssen(istart:istop,jstart:jstop),start(1:3),edges(1:3))
      if (status .NE. NF90_NOERR) then
         LEVEL3 "read_restart_ncdf(): setting ssen=0"
         ssen=_ZERO_
      else
         call update_2d_halo(ssen,ssen,az,imin,jmin,imax,jmax,H_TAG)
         call wait_halo(H_TAG)
      end if

      status = &
      nf90_get_var(ncid,ssun_id,ssun(istart:istop,jstart:jstop),start(1:3),edges(1:3))
      if (status .NE. NF90_NOERR) then
         LEVEL3 "read_restart_ncdf(): setting ssun=0"
         ssun=_ZERO_
         where (au .eq. 0)
            ssun=10.2
         end where
      else
         call update_2d_halo(ssun,ssun,au,imin,jmin,imax,jmax,U_TAG)
         call wait_halo(U_TAG)
      end if

      status = &
      nf90_get_var(ncid,ssvn_id,ssvn(istart:istop,jstart:jstop),start(1:3),edges(1:3))
      if (status .NE. NF90_NOERR) then
         LEVEL3 "read_restart_ncdf(): setting ssvn=0"
         ssvn=_ZERO_
         where (av .eq. 0)
            ssvn=10.2
         end where
      else
         call update_2d_halo(ssvn,ssvn,av,imin,jmin,imax,jmax,V_TAG)
         call wait_halo(V_TAG)
      end if

      status = &
      nf90_get_var(ncid,sseo_id,sseo(istart:istop,jstart:jstop),start(1:3),edges(1:3))
      if (status .NE. NF90_NOERR) then
         LEVEL3 "read_restart_ncdf(): setting sseo=0"
         sseo=_ZERO_
      else
         call update_2d_halo(sseo,sseo,az,imin,jmin,imax,jmax,H_TAG)
         call wait_halo(H_TAG)
      end if

      status = &
      nf90_get_var(ncid,ssuo_id,ssuo(istart:istop,jstart:jstop),start(1:3),edges(1:3))
      if (status .NE. NF90_NOERR) then
         LEVEL3 "read_restart_ncdf(): setting ssuo=0"
         ssuo=_ZERO_
      else
         call update_2d_halo(ssuo,ssuo,au,imin,jmin,imax,jmax,U_TAG)
         call wait_halo(U_TAG)
      end if

      status = &
      nf90_get_var(ncid,ssvo_id,ssvo(istart:istop,jstart:jstop),start(1:3),edges(1:3))
      if (status .NE. NF90_NOERR) then
         LEVEL3 "read_restart_ncdf(): setting ssvo=0"
         ssvo=_ZERO_
      else
         call update_2d_halo(ssvo,ssvo,av,imin,jmin,imax,jmax,V_TAG)
         call wait_halo(V_TAG)
      end if

      if (Uint_id .ne. -1) then
         status = &
         nf90_get_var(ncid,Uint_id,Uint(istart:istop,jstart:jstop),start(1:3),edges(1:3))
         if (status .NE. NF90_NOERR) go to 10
         call update_2d_halo(Uint,Uint,au,imin,jmin,imax,jmax,U_TAG)
         call wait_halo(U_TAG)
      else
         LEVEL3 "read_restart_ncdf(): setting Uint=0"
         LEVEL3 "                     (invalid within 2d cycle!)"
         Uint=_ZERO_
      end if

      if (Vint_id .ne. -1) then
         status = &
         nf90_get_var(ncid,Vint_id,Vint(istart:istop,jstart:jstop),start(1:3),edges(1:3))
         if (status .NE. NF90_NOERR) go to 10
         call update_2d_halo(Vint,Vint,av,imin,jmin,imax,jmax,V_TAG)
         call wait_halo(V_TAG)
      else
         LEVEL3 "read_restart_ncdf(): setting Vint=0"
         LEVEL3 "                     (invalid within 2d cycle!)"
         Vint=_ZERO_
      end if

      if (Uadv_id .ne. -1) then
         status = &
         nf90_get_var(ncid,Uadv_id,Uadv(istart:istop,jstart:jstop),start(1:3),edges(1:3))
         if (status .NE. NF90_NOERR) go to 10
         call update_2d_halo(Uadv,Uadv,au,imin,jmin,imax,jmax,U_TAG)
         call wait_halo(U_TAG)
      else if (Uinto_id .ne. -1) then
         status = &
         nf90_get_var(ncid,Uinto_id,Uadv(istart:istop,jstart:jstop),start(1:3),edges(1:3))
         if (status .NE. NF90_NOERR) go to 10
         call update_2d_halo(Uadv,Uadv,au,imin,jmin,imax,jmax,U_TAG)
         call wait_halo(U_TAG)
      else
         LEVEL3 "read_restart_ncdf(): setting Uadv=0"
         Uadv=_ZERO_
      end if

      if (Vadv_id .ne. -1) then
         status = &
         nf90_get_var(ncid,Vadv_id,Vadv(iloc:ilen,jloc:jlen),start,edges)
         if (status .NE. NF90_NOERR) go to 10
         call update_2d_halo(Vadv,Vadv,av,imin,jmin,imax,jmax,V_TAG)
         call wait_halo(V_TAG)
      else if (Vinto_id .ne. -1) then
         status = &
         nf90_get_var(ncid,Vinto_id,Vadv(iloc:ilen,jloc:jlen),start,edges)
         if (status .NE. NF90_NOERR) go to 10
         call update_2d_halo(Vadv,Vadv,av,imin,jmin,imax,jmax,V_TAG)
         call wait_halo(V_TAG)
      else
         LEVEL3 "read_restart_ncdf(): setting Vadv=0"
         Vadv=_ZERO_
      end if

      status = &
      nf90_get_var(ncid,uu_id,uu(istart:istop,jstart:jstop,0:kmax),start(1:3),edges(1:3))
      if (status .NE. NF90_NOERR) then
         LEVEL3 "read_restart_ncdf(): setting uu=0"
         uu=_ZERO_
      else
         call update_3d_halo(uu,uu,au,imin,jmin,imax,jmax,kmax,U_TAG)
         call wait_halo(U_TAG)
      end if

      status = &
      nf90_get_var(ncid,vv_id,vv(istart:istop,jstart:jstop,0:kmax),start(1:3),edges(1:3))
      if (status .NE. NF90_NOERR) then
         LEVEL3 "read_restart_ncdf(): setting vv=0"
         vv=_ZERO_
      else
         call update_3d_halo(vv,vv,av,imin,jmin,imax,jmax,kmax,V_TAG)
         call wait_halo(V_TAG)
      end if

      status = &
      nf90_get_var(ncid,ww_id,ww(istart:istop,jstart:jstop,0:kmax),start(1:3),edges(1:3))
      if (status .NE. NF90_NOERR) then
         LEVEL3 "read_restart_ncdf(): setting ww=0"
         ww=_ZERO_
      else
         call update_3d_halo(ww,ww,az,imin,jmin,imax,jmax,kmax,H_TAG)
         call wait_halo(H_TAG)
      end if

      status = &
      nf90_get_var(ncid,uuEx_id,uuEx(istart:istop,jstart:jstop,0:kmax),start(1:3),edges(1:3))
      if (status .NE. NF90_NOERR) then
         LEVEL3 "read_restart_ncdf(): setting uuEx=0"
         uuEx=_ZERO_
      else
         call update_3d_halo(uuEx,uuEx,au,imin,jmin,imax,jmax,kmax,U_TAG)
         call wait_halo(U_TAG)
      end if

      status = &
      nf90_get_var(ncid,vvEx_id,vvEx(istart:istop,jstart:jstop,0:kmax),start(1:3),edges(1:3))
      if (status .NE. NF90_NOERR) then
         LEVEL3 "read_restart_ncdf(): setting vvEx=0"
         vvEx=_ZERO_
      else
         call update_3d_halo(vvEx,vvEx,av,imin,jmin,imax,jmax,kmax,V_TAG)
         call wait_halo(V_TAG)
      end if

!     tke is required
      status = &
      nf90_get_var(ncid,tke_id,tke(istart:istop,jstart:jstop,0:kmax),start(1:3),edges(1:3))
      if (status .NE. NF90_NOERR) go to 10
      call update_3d_halo(tke,tke,az,imin,jmin,imax,jmax,kmax,H_TAG)
      call wait_halo(H_TAG)

!     eps is required
      status = &
      nf90_get_var(ncid,eps_id,eps(istart:istop,jstart:jstop,0:kmax),start(1:3),edges(1:3))
      if (status .NE. NF90_NOERR) go to 10
      call update_3d_halo(eps,eps,az,imin,jmin,imax,jmax,kmax,H_TAG)
      call wait_halo(H_TAG)

!     num is required
      status = &
      nf90_get_var(ncid,num_id,num(istart:istop,jstart:jstop,0:kmax),start(1:3),edges(1:3))
      if (status .NE. NF90_NOERR) go to 10
      call update_3d_halo(num,num,az,imin,jmin,imax,jmax,kmax,H_TAG)
      call wait_halo(H_TAG)

!     nuh is required
      status = &
      nf90_get_var(ncid,nuh_id,nuh(istart:istop,jstart:jstop,0:kmax),start(1:3),edges(1:3))
      if (status .NE. NF90_NOERR) go to 10
      call update_3d_halo(nuh,nuh,az,imin,jmin,imax,jmax,kmax,H_TAG)
      call wait_halo(H_TAG)

      if (ho_id .ne. -1) then
         status = &
         nf90_get_var(ncid,ho_id,ho(iloc:ilen,jloc:jlen,0:kmax),start,edges)
         if (status .NE. NF90_NOERR) go to 10
         call update_3d_halo(ho,ho,az,imin,jmin,imax,jmax,kmax,H_TAG)
         call wait_halo(H_TAG)
      endif

      if (hn_id .ne. -1) then
         status = &
         nf90_get_var(ncid,hn_id,hn(istart:istop,jstart:jstop,0:kmax),start(1:3),edges(1:3))
         if (status .NE. NF90_NOERR) go to 10
         call update_3d_halo(hn,hn,az,imin,jmin,imax,jmax,kmax,H_TAG)
         call wait_halo(H_TAG)
      endif

#ifndef NO_BAROCLINIC
      if (runtype .ge. 3)  then
!        T is required
         status = &
         nf90_get_var(ncid,T_id,T(istart:istop,jstart:jstop,0:kmax),start(1:3),edges(1:3))
         if (status .NE. NF90_NOERR) go to 10
         call update_3d_halo(T,T,az,imin,jmin,imax,jmax,kmax,H_TAG)
         call wait_halo(H_TAG)

!        S is required
         status = &
         nf90_get_var(ncid,S_id,S(istart:istop,jstart:jstop,0:kmax),start(1:3),edges(1:3))
         if (status .NE. NF90_NOERR) go to 10
         call update_3d_halo(S,S,az,imin,jmin,imax,jmax,kmax,H_TAG)
         call wait_halo(H_TAG)
      end if
#endif

      if (nonhyd_method .ne. 0) then
         if (nonhyd_method .eq. 1) then
            if (minus_bnh_id .eq. -1) then
               LEVEL3 "read_restart_ncdf(): setting minus_bnh=0"
               minus_bnh = _ZERO_
            else
               status = &
               nf90_get_var(ncid,minus_bnh_id,minus_bnh(istart:istop,jstart:jstop,0:kmax),start(1:3),edges(1:3))
               if (status .NE. NF90_NOERR) go to 10
               call update_3d_halo(minus_bnh,minus_bnh,az,imin,jmin,imax,jmax,kmax,H_TAG)
               call wait_halo(H_TAG)
            end if
         end if

         if (wco_id .eq. -1) then
            LEVEL3 "read_restart_ncdf(): setting wco=0"
            wco = _ZERO_
         else
            status = &
            nf90_get_var(ncid,wco_id,wco(istart:istop,jstart:jstop,0:kmax),start(1:3),edges(1:3))
            if (status .NE. NF90_NOERR) go to 10
            call update_3d_halo(wco,wco,az,imin,jmin,imax,jmax,kmax,H_TAG)
            call wait_halo(H_TAG)
         end if
      end if

#ifdef SPM
      if(spm_calc) then
        if (spm_hotstart) then
         status = &
         nf90_get_var(ncid,spm_id,spm(istart:istop,jstart:jstop,0:kmax),start(1:3),edges(1:3))
         if (status .NE. NF90_NOERR) then
            LEVEL3 "read_restart_ncdf(): setting spm=0"
            spm=_ZERO_
         else
            call update_3d_halo(spm,spm,az,imin,jmin,imax,jmax,kmax,H_TAG)
            call wait_halo(H_TAG)
         end if

         status = &
         nf90_get_var(ncid,spmpool_id,spm_pool(istart:istop,jstart:jstop),start(1:3),edges(1:3))
         if (status .NE. NF90_NOERR) then
            LEVEL3 "read_restart_ncdf(): setting spmpool=0"
            spm_pool=_ZERO_
         else
            call update_2d_halo(spm_pool,spm_pool,az,imin,jmin,imax,jmax,H_TAG)
            call wait_halo(H_TAG)
         end if
        else
         LEVEL3 'spm variables not read from hotstart file'
         LEVEL3 'set spm_init_method=0 to read them from hotstart file'
        end if
      end if
#endif
#ifdef _FABM_
      if(allocated(fabm_pel) .and. fabm_init_method .eq. 0) then
         start(4) = 1;  edges(4) = size(fabm_pel,4)
         status = &
         nf90_get_var(ncid,fabm_pel_id,fabm_pel(istart:istop,jstart:jstop,0:kmax,:),start(1:3),edges(1:3))
         if (status .NE. NF90_NOERR) go to 10

         start(3) = 1;  edges(3) = size(fabm_ben,3)

         if (fabm_ben_id .gt. 0) then
            status = &
            nf90_get_var(ncid,fabm_ben_id,fabm_ben(istart:istop,jstart:jstop,:),start(1:3),edges(1:3))
            if (status .NE. NF90_NOERR) go to 10
         end if
      end if
#endif
#ifdef GETM_BIO
      if(bio_calc .and. bio_init_method .eq. 0) then

         start(1) = 1;  edges(1) = numc
#ifdef _READ_HOT_HALOS_
STDERR 'needs a fix here - read_restart_ncdf()'
call getm_error('read_restart_ncdf()','needs a fix here HOT-HALO plus BIO problem')
stop
         start(2) = il; edges(2) = numc
         start(3) = jl; edges(3) = numc
#else
         start(2) = il; edges(2) = ih-il+1
         start(3) = jl; edges(3) = jh-jl+1
#endif
         start(4) = 1;  edges(4) = kmax+1

         status = &
         nf90_get_var(ncid,bio_id,cc3d(1:numc,istart:istop,jstart:jstop,0:kmax),start(1:3),edges(1:3))
         if (status .NE. NF90_NOERR) go to 10

      end if
#endif
   end if
#endif

   LEVEL3 'Done reading - closing hotstart file'
   status = nf90_close(ncid)
   if (status .NE. NF90_NOERR) go to 10

   return

10 LEVEL3 'Something bad happened'
   FATAL 'read_restart_ncdf: ',nf90_strerror(status)
   call getm_error('read_restart_ncdf()',nf90_strerror(status))
   stop
   return

   end subroutine read_restart_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Karsten Bolding (BBH)                           !
!-----------------------------------------------------------------------
