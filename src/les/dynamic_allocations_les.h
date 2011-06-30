   allocate(Am_2d(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (Am_2d)'

   allocate(AmX_2d(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (AmX_2d)'

#ifndef NO_3D
   if (runtype .ne. 1) then
      allocate(Am_3d(I3DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (Am_3d)'

      allocate(AmX_3d(I3DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (AmX_3d)'

      allocate(AmU_3d(I3DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (AmU_3d)'

      allocate(AmV_3d(I3DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (AmV_3d)'
   end if
#endif
