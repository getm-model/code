  allocate(uu(I3DFIELD),stat=rc)    ! 3D field for u-velocity
  if (rc /= 0) stop 'init_3d: Error allocating memory (uu)'

  allocate(vv(I3DFIELD),stat=rc)    ! 3D field for v-velocity
  if (rc /= 0) stop 'init_3d: Error allocating memory (vv)'

  allocate(ww(I3DFIELD),stat=rc)    ! 3D field for w-velocity
  if (rc /= 0) stop 'init_3d: Error allocating memory (ww)'

#ifdef STRUCTURE_FRICTION
  allocate(sf(I3DFIELD),stat=rc)    ! 3D field for velocity in T-points
  if (rc /= 0) stop 'init_3d: Error allocating memory (sf)'
#endif

  allocate(ho(I3DFIELD),stat=rc)    ! 3D field for old box height (z-column)
  if (rc /= 0) stop 'init_3d: Error allocating memory (ho)'

  allocate(hn(I3DFIELD),stat=rc)    ! 3D field for new box height (z-column)
  if (rc /= 0) stop 'init_3d: Error allocating memory (hn)'

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

  allocate(buoy(I3DFIELD),stat=rc)  ! 3D field for buoyancy
  if (rc /= 0) stop 'init_3d: Error allocating memory (buoy)'

  allocate(idpdx(I3DFIELD),stat=rc) ! Internal pressure gradient - x
  if (rc /= 0) stop 'init_3d: Error allocating memory (idpdx)'

  allocate(idpdy(I3DFIELD),stat=rc) ! Internal pressure gradient - y
  if (rc /= 0) stop 'init_3d: Error allocating memory (idpdy)'

  allocate(rad(I3DFIELD),stat=rc) ! Solar radiation
  if (rc /= 0) stop 'init_3d: Error allocating memory (rad)'

#endif

#ifdef SPM
  allocate(spm(I3DFIELD),stat=rc) ! Suspended particulate matter
  if (rc /= 0) stop 'init_3d: Error allocating memory (spm)'

  allocate(spm_ws(I3DFIELD),stat=rc) ! Sinking velocity
  if (rc /= 0) stop 'init_3d: Error allocating memory (spm_ws)'

  allocate(spm_pool(I2DFIELD),stat=rc) ! Pool of spm
  if (rc /= 0) stop 'init_3d: Error allocating memory (spm_pool)'
#endif

  allocate(light(I3DFIELD),stat=rc) ! light advection velocity
  if (rc /= 0) stop 'init_3d: Error allocating memory (light)'

#ifdef UV_TVD
  allocate(uadv(I3DFIELD),stat=rc) ! u advection velocity
  if (rc /= 0) stop 'init_3d: Error allocating memory (uadv)'

  allocate(vadv(I3DFIELD),stat=rc) ! v advection velocity
  if (rc /= 0) stop 'init_3d: Error allocating memory (vadv)'

  allocate(wadv(I3DFIELD),stat=rc) ! w advection velocity
  if (rc /= 0) stop 'init_3d: Error allocating memory (wadv)'

  allocate(huadv(I3DFIELD),stat=rc) ! layer height for u advection
  if (rc /= 0) stop 'init_3d: Error allocating memory (huadv)'

  allocate(hvadv(I3DFIELD),stat=rc) ! layer height for v advection
  if (rc /= 0) stop 'init_3d: Error allocating memory (hvadv)'

  allocate(hoadv(I3DFIELD),stat=rc) ! old layer height for advection
  if (rc /= 0) stop 'init_3d: Error allocating memory (hoadv)'

  allocate(hnadv(I3DFIELD),stat=rc) ! new layer height for advection
  if (rc /= 0) stop 'init_3d: Error allocating memory (hnadv)'
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

  allocate(rru(I2DFIELD),stat=rc)   ! Bottom drag term in u-vel. points (3D)   
  if (rc /= 0) stop 'init_3d: Error allocating memory (rru)'

  allocate(rrv(I2DFIELD),stat=rc)   ! Bottom drag term in v-vel. points (3D)   
  if (rc /= 0) stop 'init_3d: Error allocating memory (rrv)'

  allocate(taus(I2DFIELD),stat=rc)  ! Absolute Value of surface stress
  if (rc /= 0) stop 'init_3d: Error allocating memory (taus)'

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
