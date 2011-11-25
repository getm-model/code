#ifdef USE_BREAKS
   allocate(break_mask(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (iteration_mask)'

   allocate(break_stat(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (iteration_stat)'
#endif

   allocate(D(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (D)'

   allocate(DU(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (DU)'

   allocate(DV(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (DV)'

   allocate(z(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (z)'

   allocate(zo(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (zo)'

   allocate(U(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (U)'

   allocate(V(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (V)'

   allocate(UEx(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (UEx)'

   allocate(VEx(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (VEx)'

   allocate(fU(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (fU)'

   allocate(fV(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (fV)'

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

      allocate(Uinto(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (Uinto)'

      allocate(Vinto(E2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_2d: Error allocating memory (Vinto)'
   end if

   allocate(res_du(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (res_du)'

   allocate(res_u(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (res_u)'

   allocate(res_dv(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (res_dv)'

   allocate(res_v(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (res_v)'

!kbk
   allocate(ruu(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (ruu)'

   allocate(rvv(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (rvv)'

   allocate(PP(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (PP)'
!kbk

   allocate(SlUx(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (SlUx)'

   allocate(SlVx(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (SlVx)'

   allocate(Slru(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (Slru)'

   allocate(Slrv(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (Slrv)'

   allocate(zub(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (zub)'

   allocate(zvb(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (zvb)'

   allocate(zub0(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (zub0)'

   allocate(zvb0(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (zvb0)'

   allocate(surfdiv(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (surfdiv)'

   allocate(fwf(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (surfdiv)'

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

!  input arrays for do_advection
   allocate(fadv(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (fadv)'

   allocate(Uadv(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (Uadv)'

   allocate(Vadv(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (Vadv)'

   allocate(DUadv(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (DUadv)'

   allocate(DVadv(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_2d: Error allocating memory (DVadv)'

