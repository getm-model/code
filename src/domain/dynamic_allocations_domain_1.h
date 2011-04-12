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
