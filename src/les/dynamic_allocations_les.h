   allocate(Am2d(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (Am2d)'

   allocate(AmX2d(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (AmX2d)'

#ifndef NO_3D
   if (runtype .ne. 1) then
      allocate(Am3d(I3DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (Am3d)'

      allocate(AmX3d(I3DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (AmX3d)'

      allocate(AmU3d(I3DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (AmU3d)'

      allocate(AmV3d(I3DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (AmV3d)'
   end if
#endif
