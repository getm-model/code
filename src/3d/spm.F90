#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  suspended_matter
!
! !INTERFACE:
   module suspended_matter
!
! !DESCRIPTION:
!  This model for Suspended Particulate Matter (SPM) considers a single
!  class of non-cohesive SPM particles that do not interact with the mean
!  flow (no density effect of SPM is taken into account by default).
!  The concentration $C$ of SPM is modelled with the tracer equation.
!  At the bottom, the net SPM flux is the residual of erosion and
!  sedimentation fluxes:
!  \begin{equation}\label{Bottom_SPM}
!  \begin{array}{l}
!  -w_sC-\partial_z(\nu_t' \partial_z C)=F_e-F_s,
!  \end{array}
!  \end{equation}
!  where erosion and sedimentation fluxes are modelled following
!  \cite{KRONE62} as functions of the bottom shear stress $\tau_b$.
!  In (\ref{Bottom_SPM}), $w_s$ is a positive settling velocity. So far,
!  GETM is only coded for constant settling velocities. 
!  The erosion flux is only non-zero when
!  the bottom shear stress exceeds a critical shear stress $\tau_{ce}$:
!  \begin{equation}\label{ero_flux}
!  F_e = \left\{
!  \begin{array}{ll}
!  \displaystyle
!  \max \left\{\frac {c_e} {\rho_0}(|\tau_b|-\tau_{ce}),0\right\},
!  & \mbox{ for } B>0 \mbox{ and } |\tau_b|>\tau_{ce} \\ \\
!  0, & \mbox{ else}
!  \end{array}
!  \right.
!  \end{equation}
!  with $c_e$ erosion constant with units kg\,s\,m$^{-4}$
!  and the fluff layer SPM content $B$ (see below).
!  The sedimentation flux is only non-zero for
!  bottom shear stresses smaller than a critical shear stress $\tau_{cs}$.
!  This flux is limited by the near bottom concentration $C_b$:
!  \begin{equation}\label{sed_flux}
!  F_s = \max\left\{\frac {w_s C_b} {\tau_{cs}} (\tau_{cs}-|\tau_b|),0\right\}.
!  \end{equation}
!  Critical shear stresses for erosion and sedimentation
!  ($\tau_{ce}$ and $\tau_{cs}$ have as units N\,m$^{-2}$).
!  However, the SPM flux between the water column and the bed
!  may be switched off by setting {\tt spm\_method} in {\tt spm.inp} to zero.
!  A pool $B$ of non-dynamic particulate matter (fluff layer)
!  is assumed in order to take into account the
!  effects of depletion of erodible material at the bottom.
!  A horizontally homogeneous distribution with $B=B_0$ kg\,m$^{-2}$
!  is initially assumed. Sedimentation and erosion fill and
!  empty this pool, respectively:
!  \begin{equation}\label{SPM_pool}
!  \begin{array}{l}
!  \partial_t(B) = F_s-F_e
!  \end{array}
!  \end{equation}
!  and the erosion flux is constricted by the availability of SPM from
!  the pool (see eq.\ (\ref{ero_flux})).
!  The erosion and sedimentation fluxes are discretised using the
!  quasi-implicit \cite{PATANKAR80} approach, which guarantees positivity
!  of SPM, but only in the diffusion step, negative values might appear
!  after the advection step, although these negative values should be small.
!  The settling of SPM is linearly reduced towards zero when the water
!  depth is between the critical and the minimum water depth. This is
!  done by means of multiplication of the settling velocity with $\alpha$,
!  (see the definition in equation (\ref{alpha})). 
!
!  It is possible to take into account the impact of sediments on density by
!  setting {\tt spm\_dens} to {\tt .true}. The modified density is computed as:
!
!  \begin{equation}\label{SPM_density}
!  {\rho}={\rho}_{T,S,p}+\left(1- \frac{{\rho}_{T,S,p}}{{\rho}_{spm}}\right) C.
!  \end{equation}
!
! !USES:
   use exceptions
   use domain, only: imin,jmin,imax,jmax,kmax,ioff,joff
#ifdef TRACER_POSITIVE
   use m2d, only : z,D
#endif
   use domain, only: H,az
   use parameters, only: rho_0,g
   use variables_3d, only: hn,taub,adv_schemes,spm,spm_ws,spm_pool
   use halo_zones, only: update_3d_halo,wait_halo,D_TAG
   IMPLICIT NONE
!
   private
!
! !PUBLIC DATA MEMBERS:
   public init_spm, do_spm
   logical, public           :: spm_calc=.false.
   logical, public           :: spm_save=.true.
   logical, public           :: spm_hotstart=.false.
!
! !PRIVATE DATA MEMBERS:
   integer                 :: spm_method=1
   integer                 :: spm_init_method=1, spm_format=2
   character(len=PATH_MAX) :: spm_file="spm.nc"
   character(len=32)       :: spm_name='spm'
   integer                 :: spm_hor_adv=1,spm_ver_adv=1,spm_adv_split=0
   REALTYPE                :: spm_AH = -_ONE_
   REALTYPE                :: spm_const= _ZERO_
   REALTYPE                :: spm_init= _ZERO_
   integer                 :: spm_ws_method = 0
   REALTYPE                :: spm_ws_const=0.001
   REALTYPE                :: spm_erosion_const, spm_tauc_sedimentation
   REALTYPE                :: spm_tauc_erosion, spm_pool_init
   REALTYPE                :: spm_porosity=_ZERO_
   REALTYPE                :: spm_rho= 2650.
   logical                 :: spm_dens=.false.
!  For erosion-sedimentation flux
   REALTYPE                :: Erosion_flux , Sedimentation_flux
   logical                 :: erosed_flux =.false.
! For flocculation (not yet in namelist)
   REALTYPE                :: spm_gellingC=0.08        !(g/l or kg/m3)
   REALTYPE                :: spm_part_density=2650.   !(g/l or kg/m3)
   integer                 :: spm_mfloc=4
!
! !REVISION HISTORY:
!  Original author(s): Manuel Ruiz Villarreal, Karsten Bolding 
!                      and Hans Burchard
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_spm
!
! !INTERFACE:
   subroutine init_spm(nml_file,runtype)
!
! !DESCRIPTION:
!
! Here, the suspended matter equation is initialised. First, the namelist
! {\tt spm} is read from {\tt getm.inp}. Then, depending on the
! {\tt spm\_init\_method}, the suspended matter field is read from a
! hotstart file ({\tt spm\_init\_method}=0), initialised with a constant value
! ({\tt spm\_init\_method}=1), initialised and interpolated
! with horizontally homogeneous
! suspended matter from a given suspended matter profile 
! ({\tt spm\_init\_method}=2),
! or read in and interpolated from a 3D netCDF field 
! ({\tt spm\_init\_method}=3).
! Then, some specifications for the SPM bottom pool are given, such as that
! there should be no initial SPM pool on tidal flats.  
!
! As the next step, a number of sanity checks is performed for the chosen
! suspended matter advection schemes.
!
! Finally, the settling velocity is directly prescibed or calculated 
! by means of the \cite{ZANKE77} formula.
!
! !USES:
!  For initialization of spm in intertidal flats
   use domain,only: min_depth
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)   :: nml_file
!   logical                   :: hotstart_spm
   integer, intent(in)            :: runtype
!
! !REVISION HISTORY:
!  See revision for the module
!
! !LOCAL VARIABLES:
   integer         :: i,j,k,n
   integer         :: rc
   integer, parameter     :: nmax=100
   REALTYPE        :: zlev(nmax),prof(nmax)
!  No initial pool of spm at intertidal flats
   logical         :: intertidal_spm0=.false.
   namelist /spm_nml/   spm_calc,spm_save,spm_method,spm_init_method,  &
                        spm_const,spm_format,spm_file,spm_name,        &
                        spm_hor_adv, spm_ver_adv,spm_adv_split,        &
                        spm_AH,spm_ws_method,spm_ws_const,             &
                        spm_erosion_const, spm_tauc_sedimentation,     &
                        spm_tauc_erosion, spm_porosity, spm_pool_init, &
                        spm_rho,spm_dens
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_spm() # ',Ncall
#endif

   LEVEL1 'init_spm()'
   open(NAMLST2,status='unknown',file=trim(nml_file))
   read(NAMLST2,spm_nml)
   close(NAMLST2)

   if (spm_calc) then
      if (spm_dens .and. runtype .ne. 4)    &
         stop 'spm_dens=.true. only makes sense for runtype=4'
      LEVEL2 'spm_init_method= ',spm_init_method
      select case (spm_init_method)
         case(0)
            LEVEL3 'getting initial SPM from hotstart'
            spm_hotstart=.true.
         case(1)
            LEVEL3 'setting to constant value'
            forall(i=imin:imax,j=jmin:jmax, az(i,j) .ne. 0) &
                   spm(i,j,:) = spm_const
         case(2)
            LEVEL3 'using profile'
            call read_profile(spm_file,nmax,zlev,prof,n)
            call ver_interpol(n,zlev,prof,imin,jmin,imax,jmax,kmax, &
                              az,H,hn,spm)
         case(3)
            LEVEL3 'interpolating from 3D field'
            call get_field(spm_file,spm_name,spm)
         case default
            FATAL 'Not valid spm_init_method specified'
            stop 'init_spm'
      end select

      LEVEL2 'spm_method= ',spm_method
      select case (spm_method)
          case(0)
             erosed_flux = .false.
          case(1)
!         erosion-sedimentation flux
             erosed_flux = .true.
!         A pool of nondynamic particulate matter is assumed that initially
!         has a concentration per m of spm_pool_init
             if(spm_init_method /= 0) then
                forall(i=imin:imax,j=jmin:jmax, az(i,j) .eq. 1) &
                   spm_pool(i,j) = spm_pool_init
                   if(.not. intertidal_spm0) then
                      LEVEL3 'No spm pool in intertidal flats'
                      forall(i=imin:imax,j=jmin:jmax, H(i,j) <= 1.35 ) &
                         spm_pool(i,j) = _ZERO_
                   end if
             end if
          case default
             FATAL 'Not valid spm_method specified'
             stop 'init_spm'
      end select

      LEVEL2 'spm_hor_adv=   ',spm_hor_adv
      LEVEL2 'spm_ver_adv=   ',spm_ver_adv
      LEVEL2 'spm_adv_split= ',spm_adv_split

      if(spm_hor_adv .eq. 1) then
         spm_adv_split=-1
         if(spm_ver_adv .ne. 1) then
            LEVEL3 "setting spm_ver_adv to 1 - since spm_hor_adv is 1"
            spm_ver_adv=1
         end if
      end if
      LEVEL3 "horizontal: ",trim(adv_schemes(spm_hor_adv))," of spm"
      LEVEL3 "vertical:   ",trim(adv_schemes(spm_ver_adv))," of spm"

      select case (spm_adv_split)
         case (-1)
         case (0)
            select case (spm_hor_adv)
               case (2,3,4,5,6)
               case default
                  call getm_error("init_3d()", &
                       "spm_adv_split=0: spm_hor_adv not valid (2-6)")
            end select
            select case (spm_ver_adv)
               case (2,3,4,5,6)
               case default
                  call getm_error("init_3d()", &
                       "spm_adv_split=0: spm_ver_adv not valid (2-6)")
            end select
            LEVEL3 "1D split --> full u, full v, full w"
         case (1)
            select case (spm_hor_adv)
               case (2,3,4,5,6)
               case default
                  call getm_error("init_3d()", &
                       "spm_adv_split=1: spm_hor_adv not valid (2-6)")
            end select
            select case (spm_ver_adv)
               case (2,3,4,5,6)
               case default
                  call getm_error("init_3d()", &
                       "spm_adv_split=1: spm_ver_adv not valid (2-6)")
            end select
            LEVEL3 "1D split --> half u, half v, full w, half v, half u"
         case (2)
            select case (spm_hor_adv)
            case (2,7)
            case default
               call getm_error("init_3d()", &
                       "spm_adv_split=2: spm_hor_adv not valid (2,7)")
            end select
            select case (spm_ver_adv)
               case (2,3,4,5,6)
               case default
                  call getm_error("init_3d()", &
                       "spm_adv_split=2: spm_ver_adv not valid (2-6)")
            end select
            LEVEL3 "2D-hor, 1D-vert split --> full uv, full w"
         case default
      end select
      spm_ws = _ZERO_

!  Compute settling velocity
      LEVEL2 'spm_ws_method= ',spm_ws_method
      select case(spm_ws_method)
         case(0)
!        This will have a function if spm_ws is changing and will go to hotstart
            LEVEL3 'Case 0 not valid for spm_ws_method'
            stop 'init_spm'
         case(1) !constant
            LEVEL3 'spm_ws_const= ',spm_ws_const
            do k=0,kmax
               spm_ws(:,:,k) = spm_ws_const
            end do
         case(2)
!     Zanke formula for fall velocity of sediment
            FATAL 'Zanke method for settling velocity not yet coded'
            stop  'init_spm'
!         case(3)
!            LEVEL3 'spm_ws_const= ',spm_ws_const
!            do k=0,kmax
!               spm_ws(:,:,k) = spm_ws_const
!            end do
      end select
      LEVEL2 'spm_erosion_const= ',spm_erosion_const
      LEVEL2 'spm_tauc_sedimentation= ',spm_tauc_sedimentation
      LEVEL2 'spm_tauc_erosion= ',spm_tauc_erosion
      LEVEL2 'spm_pool_init= ',spm_pool_init
#ifdef TRACER_POSITIVE
      LEVEL3 'Positivity of spm concentration will be tested'
      LEVEL3 'and negative values written to fort.101'
      write(101,*) 'Points where spm is negative'
      write(101,*) 'Negative values of less than 1e-12 probably do not imply problems'
      write(101,*) 'i,j,k,spm(i,j,k),D(i,j),z(i,j),ww(i,j,k),ww(i,j,kmax),1 or 2'
      write(101,*) 'In column 9 1 indicates negative value after advection step'
      write(101,*) 'In column 9 2 indicates negative value after vertical diffusion step'
#endif

#ifdef INTERLEAVING_TEST
   spm(2:6,2,1:20) = 1.
#endif

   else
      LEVEL2 'No suspended sediment calculations'
      spm_save=.false.
   end if
#ifdef DEBUG
   write(debug,*) 'Leaving init_spm()'
   write(debug,*)
#endif
   return
   end subroutine init_spm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_spm - suspended matter equation \label{sec-do-spm}
!
! !INTERFACE:
   subroutine do_spm()
!
! !DESCRIPTION:
!
! Here, one time step for the suspended matter equation is performed.
! First, preparations for the call to the advection schemes are
! made, i.e.\ calculating the necessary metric coefficients
! and the relevant vertical velocity, which is here composed of the
! grid-related vertical flow velocity and the settling velocity.
! Some lines of code allow here for consideration of
! flocculation processes.
! After the call to the advection schemes, which actually perform
! the advection (and horizontal diffusion) step as an operational
! split step, the fluxes between bottom SPM pool and the suspended matter
! in the water column are calculated. 
! Afterwards, the tri-diagonal matrix for calculating the
! new suspended matter by means of a semi-implicit central scheme for the
! vertical diffusion is set up.
! There are no source terms on the right hand sides.
! The subroutine is completed by solving the tri-diagonal linear
! equation by means of a tri-diagonal solver.
! 
! Optionally, the density of the sediment-laden water may be 
! corrected by the sediment density, see eq.\ (\ref{SPM_density}).
!
! Finally, some special settings for single test cases are made via
! compiler options.
!
! !USES:
   use advection_3d, only: do_advection_3d
   use variables_3d, only: dt,cnpar,hun,hvn,ho,nuh,uu,vv,ww
#ifndef NO_BAROCLINIC
   use variables_3d, only: rho
#endif
   use domain,       only: au,av
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxu,dxv,dyu,dyv,arcd1
#ifdef ELBE_TEST
   use domain, only :dxc
#endif
#else
   use domain, only: dx,dy,ard1
#endif
   use domain, only: dry_z
   IMPLICIT NONE
!
! !LOCAL VARIABLES:
   integer         :: i,j,k,rc
   REALTYPE        :: spmtot
   REALTYPE        :: Res(0:kmax)
   REALTYPE        :: auxn(1:kmax-1),auxo(1:kmax-1)
   REALTYPE        :: a1(0:kmax),a2(0:kmax)
   REALTYPE        :: a3(0:kmax),a4(0:kmax)
   REALTYPE        :: delxu(I2DFIELD),delxv(I2DFIELD)
   REALTYPE        :: delyu(I2DFIELD),delyv(I2DFIELD)
   REALTYPE        :: area_inv(I2DFIELD)
   REALTYPE        :: bed_flux
   REALTYPE        :: c
   REALTYPE        :: volCmud,volCpart
   integer         :: k2
   logical         :: patankar=.true.
#ifdef TRACER_POSITIVE
   logical         :: kk
#endif
   REALTYPE, allocatable, dimension (:,:,:) :: ww_aux
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_spm() # ',Ncall
#endif

#if defined(SPHERICAL) || defined(CURVILINEAR)
   delxu=dxu
   delxv=dxv
   delyu=dyu
   delyv=dyv
   area_inv=arcd1
#else
   delxu=dx
   delxv=dx
   delyu=dy
   delyv=dy
   area_inv=ard1
#endif
   allocate(ww_aux(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_spm: Error allocating memory (ww_aux)'

!  Update settling velocity if flocculation is considered
   select case(spm_ws_method)
      case(3) !Floculation according to winterwerp
         do i=imin,imax
         do j=jmin,jmax
            if(az(i,j) == 1) then
               do k=1,kmax-1
                  volCmud=spm(i,j,k)/spm_gellingC
                  volCpart=spm(i,j,k)/spm_part_density
                  volCpart=1.-volCpart
                  spm_ws(i,j,k)=spm_ws_const*volCpart
                  spm_ws(i,j,k)=spm_ws(i,j,k)*(_ONE_-min(_ONE_,volCmud))**spm_mfloc
                  spm_ws(i,j,k)=spm_ws(i,j,k)/(_ONE_+2.5*volCmud)
               end do
            end if
         end do
         end do
         spm_ws(:,:,0)=spm_ws(:,:,1)
         spm_ws(:,:,kmax)=spm_ws(:,:,kmax-1)
         LEVEL3 'stop. You cannot use a variable sinking velocity '
         LEVEL3 'without changing advection_3d. The velocity ww used in '
         LEVEL3 'the advection step for continuity is only the hydrodynamical '
         stop
   end select
!  The vertical velocity to be used in the advection routine for spm is ww-ws
!  In drying grid boxes, the settling velocity is reduced.
   do i=imin,imax
      do j=jmin,jmax
         do k=1,kmax-1
            ww_aux(i,j,k) = ww(i,j,k) - spm_ws(i,j,k)*dry_z(i,j)
         end do
      end do
   end do
   call do_advection_3d(dt,spm,uu,vv,ww_aux,hun,hvn,ho,hn,   &
                        delxu,delxv,delyu,delyv,area_inv,az,au,av,     &
                        spm_hor_adv,spm_ver_adv,spm_adv_split,spm_AH)
#ifdef TRACER_POSITIVE
   kk= .false.
   do i=imin,imax
      do j=jmin,jmax
         do k=1,kmax
            if(spm(i,j,k) < _ZERO_ .and. az(i,j) >= 1) then
               if(spm(i,j,k) < -1.e-12) then
                  write(101,*) i,j,k,spm(i,j,k),D(i,j),z(i,j),ww(i,j,k),ww(i,j,kmax),'1'
               end if
               spm(i,j,k) = _ZERO_
            end if
         end do
      end do
   end do
#endif

!  Vertical diffusion of spm
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1 ) then
            if (kmax .gt. 1) then
!              We impose a flux condition on bottom sediments
!              where flux is computed as the result of erosion and sedimentation
!              from a pool of nondynamic particulate matter
               if (erosed_flux) then
!                 1-spm_porosity is the fractional bed concentration
                  if(spm_pool(i,j) > _ZERO_) then
!                    if there are sediments in the pool
                     Erosion_Flux = spm_erosion_const / rho_0                  &
                            *(1-spm_porosity) * (taub(i,j)*rho_0-spm_tauc_erosion )
                     Erosion_Flux = max(Erosion_Flux,_ZERO_)
                  else
                     Erosion_Flux = _ZERO_
                  end if
                  Sedimentation_Flux = spm_ws(i,j,1) * spm(i,j,1) *  &
                                 (1.-taub(i,j)*rho_0 / spm_tauc_sedimentation)
                  Sedimentation_Flux = max(Sedimentation_Flux,_ZERO_)
                  bed_flux = Erosion_Flux-Sedimentation_Flux
                  if (bed_flux > _ZERO_) then
                     bed_flux = min(spm_pool(i,j)/dt,bed_flux)
                  else
                     Erosion_Flux = min(spm_pool(i,j)/dt,Erosion_Flux)
                     spm_pool(i,j) = spm_pool(i,j) -  dt * Erosion_Flux
                     spm_pool(i,j) = max(spm_pool(i,j),_ZERO_)
                  end if
                  if(.not. patankar) then
                     spm_pool(i,j) = spm_pool(i,j) -  dt * bed_flux
                     spm_pool(i,j) = max(spm_pool(i,j),_ZERO_)
                  end if
               else
                  bed_flux = _ZERO_
               end if
!     Auxiliary terms, old and new time level,
               do k=1,kmax-1
                  auxo(k)=2.*(1-cnpar)*dt*nuh(i,j,k)/(hn(i,j,k+1)+hn(i,j,k))
                  auxn(k)=2.*   cnpar *dt*nuh(i,j,k)/(hn(i,j,k+1)+hn(i,j,k))
               end do

!        Matrix elements for surface layer
               k=kmax
               a1(k)=-auxn(k-1)
               a2(k)=hn(i,j,k)+auxn(k-1)
               a4(k)=spm(i,j,k)*(hn(i,j,k)-auxo(k-1))+spm(i,j,k-1)*auxo(k-1)

!        Matrix elements for inner layers
               do k=2,kmax-1
                  a1(k)=-auxn(k-1)
                  a2(k)=+hn(i,j,k)+auxn(k)+auxn(k-1)
                  a3(k)=-auxn(k  )
                  a4(k)=spm(i,j,k+1)*auxo(k)                          &
                       +spm(i,j,k  )*(hn(i,j,k)-auxo(k)-auxo(k-1))    &
                       +spm(i,j,k-1)*auxo(k-1)
               end do

!        Matrix elements for bottom layer
               k=1
               a2(k)=hn(i,j,k)+auxn(k)
               a3(k)=-auxn(k  )
               a4(k)=spm(i,j,k+1)*auxo(k)+spm(i,j,k  )*(hn(i,j,k)-auxo(k))

!        Patankar trick
               if(patankar) then
                  if (bed_flux .ge. _ZERO_) then
                     a4(k)=a4(k)+bed_flux*dt
                  else
                     a2(k)=a2(k)-dt*bed_flux/spm(i,j,k)
                  end if
               else
                  a4(k)=a4(k)+bed_flux*dt
               end if

               call getm_tridiagonal(kmax,1,kmax,a1,a2,a3,a4,Res)
               do k=1,kmax
                  spm(i,j,k)=Res(k)
               end do

!        spm_pool is updated
               if (patankar) then
                  if (bed_flux .ge. _ZERO_) then
                     spm_pool(i,j) = spm_pool(i,j) -  dt * bed_flux
                  else
                     spm_pool(i,j) = spm_pool(i,j) +  dt * spm_ws(i,j,1) * spm(i,j,1)       &
                           *max(_ZERO_, 1.-taub(i,j)*rho_0 / spm_tauc_sedimentation)
                  end if
               end if
               spm_pool(i,j) = max(_ZERO_,spm_pool(i,j))
            end if
         end if
      end do
   end do

!  Boundary conditions. This should be done via do_bdy_3d but not ready yet
#if 1

#ifdef SALTWEDGE_TEST
!  Valid for saltwedge case, lateral zero-gradient BC for WE boundary
!   call do_bdy_3d(2,spm)
   do k=1,kmax
      if(uu(imin,2,k) > _ZERO_) then
         spm(imin,2,k)=spm(imin+1,2,k)
      end if
   end do
#endif

#ifndef SALTWEDGE_TEST
!  BC at open boundary points (no flux, can be other condition)
   do i=imin,imax
      if (az(i,jmin).eq.2) spm(i,jmin,:)=spm(i,jmin+1,:)
      if (az(i,jmax).eq.2) spm(i,jmax,:)=spm(i,jmax-1,:)
   end do

   do j=jmin,jmax
      if (az(imin,j).eq.2) spm(imin,j,:)=spm(imin+1,j,:)
      if (az(imax,j).eq.2) spm(imax,j,:)=spm(imax-1,j,:)
   end do
#endif

#ifdef ELBE_TEST
   do j=jmin,jmax
      if(az(1,j) /= 0) then
          do k=1 ,kmax
             if (uu(1,j,k) .gt. _ZERO_) then
                spm(1,j,k)=max(spm(1,j,k),_ZERO_)
             else
                spm(1,j,k)=spm(1,j,k)-dt*area_inv(1,j)*dyu(1,j)*uu(1,j,k)/hun(1,j,k)*(spm(2,j,k)-spm(1,j,k))
             end if
          end do
      end if
!     River boundary
      if(az(imax,j) /= 0) then
         do k=1,kmax
            spm(imax,j,k)=_ZERO_
         end do
      end if
   end do
#endif

#endif
!  End of BC

#ifndef NO_BAROCLINIC
!  contribution of spm to density
   if (spm_dens) then
      do j=jmin,jmax
         do i=imin,imax
            if (az(i,j) .eq. 1 ) then
               do k=1,kmax
                  rho(i,j,k)=rho(i,j,k)+(_ONE_-rho(i,j,k)/spm_rho)*spm(i,j,k)
               end do
            end if
         end do
      end do
   end if
#endif

#ifdef TRACER_POSITIVE
   kk=.false.
   do i=imin,imax
      do j=jmin,jmax
         do k=1,kmax
            if(spm(i,j,k) < _ZERO_ .and. az(i,j) >= 1) then
!               if(spm(i,j,k) < -1.e-12)
               write(101,*) i,j,k,spm(i,j,k),D(i,j),z(i,j),ww(i,j,k),ww(i,j,kmax),'2'
!               kk=.true.
!               spm(i,j,k) = _ZERO_
            end if
         end do
      end do
   end do
!   if(kk .eq. .true.) stop 'Negative spm concentration'
#endif

   call update_3d_halo(spm,spm,az,imin,jmin,imax,jmax,kmax,D_TAG)
   call wait_halo(D_TAG)

#ifdef FORTRAN90
   deallocate(ww_aux,stat=rc)
   if (rc /= 0) stop 'upstream_adv: Error de-allocating memory (ww_aux)'
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving do_spm()'
   write(debug,*)
#endif
   return
   end subroutine do_spm
!EOC

!-----------------------------------------------------------------------

   end module suspended_matter

!-----------------------------------------------------------------------
! Copyright (C) 2004 - Manuel Ruiz, Hans Burchard and Karsten Bolding  !
!-----------------------------------------------------------------------
