!  mask
   allocate(az(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (az)'

   allocate(au(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (au)'

   allocate(av(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (av)'

   allocate(ax(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (ax)'


!  bathymetry
   allocate(H(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (H)'

   allocate(HU(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (HU)'
   HU = _ZERO_

   allocate(HV(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (HV)'
   HV = _ZERO_

   allocate(dry_z(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (dry_z)'
   dry_z = _ONE_

   allocate(dry_u(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (dry_u)'
   dry_u = _ONE_

   allocate(dry_v(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (dry_v)'
   dry_v = _ONE_


!  coriolis terms
   allocate(cor(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (cor)'

   allocate(coru(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (coru)'

   allocate(corv(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (corv)'


!  lat/lon
   allocate(lonc(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (lonc)'

   allocate(latc(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (latc)'

   allocate(lonx(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (lonx)'

   allocate(latx(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (latx)'

   allocate(lonu(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (lonu)'

   allocate(latu(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (latu)'

   allocate(lonv(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (lonv)'

   allocate(latv(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (latv)'



!  grid convergence
   allocate(convc(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (convc)'

   allocate(convx(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (convx)'

   allocate(angle(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (angle)'


! grid points
   allocate(xx(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (xx)'

   allocate(yx(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (yx)'

   allocate(xc(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (xc)'

   allocate(yc(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (yc)'

   allocate(xu(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (xu)'

   allocate(yu(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (yu)'

   allocate(xv(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (xv)'

   allocate(yv(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (yv)'



!  metric parameters
   allocate(dxdyc(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (dxdyc)'

   allocate(dydxc(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (dydxc)'

   allocate(dxc(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (dxc)'

   allocate(dxu(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (dxu)'

   allocate(dxv(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (dxv)'

   allocate(dxx(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (dxx)'

   allocate(dyc(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (dyc)'

   allocate(dyu(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (dyu)'

   allocate(dyv(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (dyv)'

   allocate(dyx(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (dyx)'

   allocate(arcd1(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (arcd1)'

   allocate(arud1(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (arud1)'

   allocate(arvd1(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (arvd1)'

!  bottom roughness
   allocate(z0(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (z0)'
