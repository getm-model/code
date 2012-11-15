! Remember to update this value if you add more 2D arrays.
   integer, parameter :: n2d_fields=27
!
#ifdef USE_BREAKS
   integer  break_mask(E2DFIELD)
   integer  break_stat(E2DFIELD)
#endif
   REALTYPE,dimension(E2DFIELD),target :: t_D,t_Dlast,DU,DV,t_z,t_zo
   REALTYPE U(E2DFIELD)
   REALTYPE V(E2DFIELD)
   REALTYPE UEx(E2DFIELD)
   REALTYPE VEx(E2DFIELD)
   REALTYPE fU(E2DFIELD)
   REALTYPE fV(E2DFIELD)
   REALTYPE ru(E2DFIELD)
   REALTYPE rv(E2DFIELD)
   REALTYPE Uint(E2DFIELD)
   REALTYPE Vint(E2DFIELD)
   REALTYPE res_du(E2DFIELD)
   REALTYPE res_u(E2DFIELD)
   REALTYPE res_dv(E2DFIELD)
   REALTYPE res_v(E2DFIELD)
   REALTYPE SlUx(E2DFIELD)
   REALTYPE SlVx(E2DFIELD)
   REALTYPE Slru(E2DFIELD)
   REALTYPE Slrv(E2DFIELD)
   REALTYPE fwf(E2DFIELD)
   REALTYPE fwf_int(E2DFIELD)
   REALTYPE EWbdy(jmax),ENbdy(imax),EEbdy(jmax),ESbdy(imax)

