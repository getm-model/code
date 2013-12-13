#ifdef USE_BREAKS
   allocate(break_mask(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (iteration_mask)'

   allocate(break_stat(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (iteration_stat)'
#endif

   allocate(t_zo(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (t_zo)'

   allocate(t_z(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (t_z)'

   allocate(D(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (D)'

   allocate(Dvel(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (Dvel)'

   allocate(DU(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (DU)'

   allocate(DV(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (DV)'

   allocate(U(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (U)'

   allocate(V(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (V)'

   allocate(UEx(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (UEx)'

   allocate(VEx(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (VEx)'

   allocate(ru(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (ru)'

   allocate(rv(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (rv)'

!   if (runtype .gt. 1) then
   if (runtype .gt. 0) then
      allocate(Uint(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (Uint)'

      allocate(Vint(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (Vint)'
   end if

   allocate(res_du(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (res_du)'

   allocate(res_u(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (res_u)'

   allocate(res_dv(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (res_dv)'

   allocate(res_v(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (res_v)'

   allocate(SlUx(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (SlUx)'

   allocate(SlVx(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (SlVx)'

   allocate(Slru(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (Slru)'

   allocate(Slrv(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (Slrv)'

   allocate(An(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (An)'

   allocate(AnX(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (AnX)'

   allocate(fwf(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (fwf)'

   allocate(fwf_int(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (fwf_int)'

   allocate(EWbdy(jmax),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (EWbdy)'

   allocate(ENbdy(imax),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (ENbdy)'

   allocate(EEbdy(jmax),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (EEbdy)'

   allocate(ESbdy(imax),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (ESbdy)'

