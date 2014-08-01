!  coordinate axes - in case of grid-type = 1 or 2
   allocate(xcord(imin-HALO:imax+HALO),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (xcord) - egoon'

   allocate(ycord(jmin-HALO:jmax+HALO),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (ycord)'

!  pseudo coordinate axes - in case of grid-type = 3 or 4
   allocate(xxcord(imin-HALO-1:imax+HALO),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (xxcord)'

   allocate(yxcord(jmin-HALO-1:jmax+HALO),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (yxcord)'

!  bathymetry
   allocate(H(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (H)'
   H = -10.

! grid points
   allocate(xx(E2DXFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (xx)'

   allocate(yx(E2DXFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (yx)'

   allocate(xc(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (xc)'

   allocate(yc(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (yc)'


!  lat/lon
   allocate(lonc(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (lonc)'
   lonc = -999.

   allocate(latc(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (latc)'
   latc = -999.

   allocate(lonx(E2DXFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (lonx)'
   lonx = -999.

   allocate(latx(E2DXFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (latx)'
   latx = -999.

!  grid convergence
   allocate(convc(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (convc)'
   convc = _ZERO_

   allocate(convx(E2DXFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (convx)'
   convx = _ZERO_

#if 0
   allocate(angle(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (angle)'
   angle = -999.
#endif

!  masks
   allocate(mask(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (mask)'

   allocate(az(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (az)'

   allocate(au(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (au)'

   allocate(av(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (av)'

   allocate(ax(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (ax)'

!  depth at U and V points
   allocate(HU(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (HU)'
   HU = -10.

   allocate(HV(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_domain: Error allocating memory (HV)'
   HV = -10.

!  dry fields
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
   cor = _ZERO_

   allocate(coru(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (coru)'
   coru = _ZERO_

   allocate(corv(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (corv)'
   corv = _ZERO_

!  lon and lat for U and V points
   allocate(lonu(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (lonu)'
   lonu = -999.

   allocate(latu(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (latu)'
   latu = -999.

   allocate(lonv(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (lonv)'
   lonv = -999.

   allocate(latv(E2DFIELD),stat=rc)
   if (rc /=0) stop 'init_domain: Error allocating memory (latv)'
   latv = -999.

!  positions of X-points
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
