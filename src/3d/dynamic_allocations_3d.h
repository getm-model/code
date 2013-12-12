  allocate(uu(I3DFIELD),stat=rc)    ! 3D field for u-velocity
  if (rc /= 0) stop 'init_3d: Error allocating memory (uu)'

  allocate(vv(I3DFIELD),stat=rc)    ! 3D field for v-velocity
  if (rc /= 0) stop 'init_3d: Error allocating memory (vv)'

  allocate(ww(I3DFIELD),stat=rc)    ! 3D field for w-velocity
  if (rc /= 0) stop 'init_3d: Error allocating memory (ww)'

#ifdef _MOMENTUM_TERMS_
  allocate(tdv_u(I3DFIELD),stat=rc) ! 3D field for tdv_u
  if (rc /= 0) stop 'init_3d: Error allocating memory (tdv_u)'

  allocate(adv_u(I3DFIELD),stat=rc) ! 3D field for adv_u
  if (rc /= 0) stop 'init_3d: Error allocating memory (adv_u)'

  allocate(vsd_u(I3DFIELD),stat=rc) ! 3D field for vsd_u
  if (rc /= 0) stop 'init_3d: Error allocating memory (vsd_u)'

  allocate(hsd_u(I3DFIELD),stat=rc) ! 3D field for hsd_u
  if (rc /= 0) stop 'init_3d: Error allocating memory (hsd_u)'

  allocate(cor_u(I3DFIELD),stat=rc) ! 3D field for cor_u
  if (rc /= 0) stop 'init_3d: Error allocating memory (cor_u)'

  allocate(epg_u(I3DFIELD),stat=rc) ! 3D field for epg_u
  if (rc /= 0) stop 'init_3d: Error allocating memory (epg_u)'

  allocate(ipg_u(I3DFIELD),stat=rc) ! 3D field for ipg_u
  if (rc /= 0) stop 'init_3d: Error allocating memory (ipg_u)'

  allocate(tdv_v(I3DFIELD),stat=rc) ! 3D field for tdv_v
  if (rc /= 0) stop 'init_3d: Error allocating memory (tdv_v)'

  allocate(adv_v(I3DFIELD),stat=rc) ! 3D field for adv_v
  if (rc /= 0) stop 'init_3d: Error allocating memory (adv_v)'

  allocate(vsd_v(I3DFIELD),stat=rc) ! 3D field for vsd_v
  if (rc /= 0) stop 'init_3d: Error allocating memory (vsd_v)'

  allocate(hsd_v(I3DFIELD),stat=rc) ! 3D field for hsd_v
  if (rc /= 0) stop 'init_3d: Error allocating memory (hsd_v)'

  allocate(cor_v(I3DFIELD),stat=rc) ! 3D field for cor_v
  if (rc /= 0) stop 'init_3d: Error allocating memory (cor_v)'

  allocate(epg_v(I3DFIELD),stat=rc) ! 3D field for epg_v
  if (rc /= 0) stop 'init_3d: Error allocating memory (epg_v)'

  allocate(ipg_v(I3DFIELD),stat=rc) ! 3D field for ipg_v
  if (rc /= 0) stop 'init_3d: Error allocating memory (ipg_v)'


#endif

#ifdef STRUCTURE_FRICTION
  allocate(sf(I3DFIELD),stat=rc)    ! 3D field for velocity in T-points
  if (rc /= 0) stop 'init_3d: Error allocating memory (sf)'
#endif

  allocate(ho(I3DFIELD),stat=rc)    ! 3D field for old box height (z-column)
  if (rc /= 0) stop 'init_3d: Error allocating memory (ho)'

  allocate(hn(I3DFIELD),stat=rc)    ! 3D field for new box height (z-column)
  if (rc /= 0) stop 'init_3d: Error allocating memory (hn)'

  allocate(hvel(I3DFIELD),stat=rc)    ! 3D field for intermediate box height (z-column)
  if (rc /= 0) stop 'init_3d: Error allocating memory (hvel)'

  allocate(huo(I3DFIELD),stat=rc)   ! 3D field for old box height (u-column)
  if (rc /= 0) stop 'init_3d: Error allocating memory (huo)'

  allocate(hun(I3DFIELD),stat=rc)   ! 3D field for new box height (u-column)
  if (rc /= 0) stop 'init_3d: Error allocating memory (hun)'

  allocate(hvo(I3DFIELD),stat=rc)   ! 3D field for old box height (v-column)
  if (rc /= 0) stop 'init_3d: Error allocating memory (hvo)'

  allocate(hvn(I3DFIELD),stat=rc)   ! 3D field for new box height (v-column)
  if (rc /= 0) stop 'init_3d: Error allocating memory (hvn)'

  allocate(hcc(I3DFIELD),stat=rc)    ! Hydrostatic consistency check
  if (rc /= 0) stop 'init_3d: Error allocating memory (hcc)'

  allocate(uuEx(I3DFIELD),stat=rc)  ! 3D field for explicit terms in u-equation
  if (rc /= 0) stop 'init_3d: Error allocating memory (uuEx)'

  allocate(vvEx(I3DFIELD),stat=rc)  ! 3D field for explicit terms in v-equation
  if (rc /= 0) stop 'init_3d: Error allocating memory (vvEx)'

  allocate(num(I3DFIELD),stat=rc)   ! 3D field for vertical eddy viscosity
  if (rc /= 0) stop 'init_3d: Error allocating memory (num)'

  allocate(nuh(I3DFIELD),stat=rc)   ! 3D field for vertical eddy diffusivity
  if (rc /= 0) stop 'init_3d: Error allocating memory (nuh)'

  allocate(tke(I3DFIELD),stat=rc)   ! 3D field for turbulent kinetic energy
  if (rc /= 0) stop 'init_3d: Error allocating memory (tke)'

  allocate(eps(I3DFIELD),stat=rc)   ! 3D field for dissipation
  if (rc /= 0) stop 'init_3d: Error allocating memory (eps)'

  allocate(SS(I3DFIELD),stat=rc)   ! 3D field for shear frequency
  if (rc /= 0) stop 'init_3d: Error allocating memory (SS)'

#ifndef NO_BAROCLINIC
  allocate(NN(I3DFIELD),stat=rc)   ! 3D field for buoyancy frequency
  if (rc /= 0) stop 'init_3d: Error allocating memory (NN)'

  allocate(S(I3DFIELD),stat=rc)     ! 3D field for salinity
  if (rc /= 0) stop 'init_3d: Error allocating memory (S)'

  allocate(T(I3DFIELD),stat=rc)     ! 3D field for temperature
  if (rc /= 0) stop 'init_3d: Error allocating memory (T)'

  allocate(rho(I3DFIELD),stat=rc)  ! 3D field for density
  if (rc /= 0) stop 'init_3d: Error allocating memory (rho)'

  allocate(alpha(I3DFIELD),stat=rc)  ! 3D field for alpha
  if (rc /= 0) stop 'init_3d: Error allocating memory (alpha)'

  allocate(beta(I3DFIELD),stat=rc)  ! 3D field for beta
  if (rc /= 0) stop 'init_3d: Error allocating memory (beta)'

  allocate(buoy(I3DFIELD),stat=rc)  ! 3D field for buoyancy
  if (rc /= 0) stop 'init_3d: Error allocating memory (buoy)'

  allocate(idpdx(I3DFIELD),stat=rc) ! Internal pressure gradient - x
  if (rc /= 0) stop 'init_3d: Error allocating memory (idpdx)'

  allocate(idpdy(I3DFIELD),stat=rc) ! Internal pressure gradient - y
  if (rc /= 0) stop 'init_3d: Error allocating memory (idpdy)'

  allocate(rad(I3DFIELD),stat=rc) ! Solar radiation
  if (rc /= 0) stop 'init_3d: Error allocating memory (rad)'

  allocate(light(I3DFIELD),stat=rc) ! light advection velocity
  if (rc /= 0) stop 'init_3d: Error allocating memory (light)'
#endif

#ifdef SPM
  allocate(spm(I3DFIELD),stat=rc) ! Suspended particulate matter
  if (rc /= 0) stop 'init_3d: Error allocating memory (spm)'

  allocate(spm_ws(I3DFIELD),stat=rc) ! Sinking velocity
  if (rc /= 0) stop 'init_3d: Error allocating memory (spm_ws)'

  allocate(spm_pool(I2DFIELD),stat=rc) ! Pool of spm
  if (rc /= 0) stop 'init_3d: Error allocating memory (spm_pool)'
#endif

! 2D fields in the 3D domain
  allocate(sseo(I2DFIELD),stat=rc)  ! Elevation before macro time step (z-column)
  if (rc /= 0) stop 'init_3d: Error allocating memory (sseo)'

  allocate(ssen(I2DFIELD),stat=rc)  ! Elevation after  macro time step (z-column)
  if (rc /= 0) stop 'init_3d: Error allocating memory (ssen)'

  allocate(ssuo(I2DFIELD),stat=rc)  ! Elevation before macro time step (u-column)
  if (rc /= 0) stop 'init_3d: Error allocating memory (ssuo)'

  allocate(ssun(I2DFIELD),stat=rc)  ! Elevation after  macro time step (u-column)
  if (rc /= 0) stop 'init_3d: Error allocating memory (ssun)'

  allocate(ssvo(I2DFIELD),stat=rc)  ! Elevation before macro time step (v-column)
  if (rc /= 0) stop 'init_3d: Error allocating memory (ssvo)'

  allocate(ssvn(I2DFIELD),stat=rc)  ! Elevation after  macro time step (v-column)
  if (rc /= 0) stop 'init_3d: Error allocating memory (ssvn)'

  allocate(Dn(I2DFIELD),stat=rc)  ! depth after  macro time step (z-column)
  if (rc /= 0) stop 'init_3d: Error allocating memory (Dn)'

  allocate(Dveln(I2DFIELD),stat=rc)  ! depth during  macro time step (z-column)
  if (rc /= 0) stop 'init_3d: Error allocating memory (Dveln)'

  allocate(Dun(I2DFIELD),stat=rc)  ! depth after  macro time step (u-column)
  if (rc /= 0) stop 'init_3d: Error allocating memory (Dun)'

  allocate(Dvn(I2DFIELD),stat=rc)  ! depth after  macro time step (v-column)
  if (rc /= 0) stop 'init_3d: Error allocating memory (Dvn)'

   allocate(Uadv(I2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_3d: Error allocating memory (Uadv)'

   allocate(Vadv(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_3d: Error allocating memory (Vadv)'

  allocate(rru(I2DFIELD),stat=rc)   ! Bottom drag term in u-vel. points (3D)
  if (rc /= 0) stop 'init_3d: Error allocating memory (rru)'

  allocate(rrv(I2DFIELD),stat=rc)   ! Bottom drag term in v-vel. points (3D)
  if (rc /= 0) stop 'init_3d: Error allocating memory (rrv)'

  allocate(zub(I2DFIELD),stat=rc)   ! bottom roughness length in U-point (3D)
  if (rc /= 0) stop 'init_3d: Error allocating memory (zub)'

  allocate(zvb(I2DFIELD),stat=rc)   ! bottom roughness length in V-point (3D)
  if (rc /= 0) stop 'init_3d: Error allocating memory (zvb)'

  allocate(taus(I2DFIELD),stat=rc)  ! Absolute Value of surface stress
  if (rc /= 0) stop 'init_3d: Error allocating memory (taus)'

  allocate(taubx(I2DFIELD),stat=rc)  ! x-component of bottom stress
  if (rc /= 0) stop 'init_3d: Error allocating memory (taubx)'

  allocate(tauby(I2DFIELD),stat=rc)  ! y-component of bottom stress
  if (rc /= 0) stop 'init_3d: Error allocating memory (tauby)'

  allocate(taub(I2DFIELD),stat=rc)  ! Absolute Value of bottom stress
  if (rc /= 0) stop 'init_3d: Error allocating memory (taub)'

  allocate(kmin(I2DFIELD),stat=rc) ! Bottom index for vertical z-columns
  if (rc /= 0) stop 'init_3d: Error allocating memory (kmin)'

  allocate(kumin(I2DFIELD),stat=rc) ! Bottom index for vertical u-columns
  if (rc /= 0) stop 'init_3d: Error allocating memory (kumin)'

  allocate(kvmin(I2DFIELD),stat=rc) ! Bottom index for vertical v-columns
  if (rc /= 0) stop 'init_3d: Error allocating memory (kvmin)'

! Below are allocations for Poor Mans Z-coordinates - PMZ
  allocate(kmin_pmz(I2DFIELD),stat=rc)
  if (rc /= 0) stop 'init_3d: Error allocating memory (kmin_pmz)'

  allocate(kumin_pmz(I2DFIELD),stat=rc)
  if (rc /= 0) stop 'init_3d: Error allocating memory (kumin_pmz)'

  allocate(kvmin_pmz(I2DFIELD),stat=rc)
  if (rc /= 0) stop 'init_3d: Error allocating memory (kvmin_pmz)'

! for light attenuation
  allocate(A(I2DFIELD),stat=rc)
  if (rc /= 0) stop 'init_3d: Error allocating memory (A)'

  allocate(g1(I2DFIELD),stat=rc)
  if (rc /= 0) stop 'init_3d: Error allocating memory (g1)'

  allocate(g2(I2DFIELD),stat=rc)
  if (rc /= 0) stop 'init_3d: Error allocating memory (g2)'
