!$Id: gen_bdy_ncdf.F90,v 1.1 2003-04-03 06:11:26 gotm Exp $

subroutine create_nc_file(fn,nbdyp,ts,ncid,time_id,jul_id,secs_id,elev_id)

   implicit none

   character(len=*), intent(in)        :: fn
   integer, intent(in)                 :: nbdyp
   character(len=*), intent(in)        :: ts
   integer, intent(out)                :: ncid,time_id,jul_id,secs_id,elev_id

   include 'netcdf.inc'

   integer  nbdyp_len
   integer, parameter :: time_len = NF_UNLIMITED

   integer  nbdyp_dim, time_dim

   integer  time_dims(1)
   integer  elev_dims(2)

   real  realval(1)

   integer iret

   character(len=80) :: string

   nbdyp_len = nbdyp


   iret = nf_create(trim(fn), NF_CLOBBER, ncid)
!kbk   call check_err(iret)
   write(0,*) ncid

   iret = nf_def_dim(ncid, 'nbdyp', nbdyp_len, nbdyp_dim)
   !kbkcall check_err(iret)
   iret = nf_def_dim(ncid, 'time', NF_UNLIMITED, time_dim)
   !kbkcall check_err(iret)


   time_dims(1) = time_dim
   iret = nf_def_var(ncid, 'time', NF_REAL, 1, time_dims, time_id)

   iret = nf_def_var(ncid, 'julday', NF_INT, 1, time_dims, jul_id)
   iret = nf_def_var(ncid, 'secs', NF_INT, 1, time_dims, secs_id)

   elev_dims(2) = time_dim
   elev_dims(1) = nbdyp_dim
   iret = nf_def_var(ncid, 'elev', NF_REAL, 2, elev_dims, elev_id)
!kbk    call check_err(iret)

   string = 'seconds since '//trim(ts)
   iret = nf_put_att_text(ncid,time_id, 'units', len_trim(string), trim(string))
!kbk   call check_err(iret)

   iret = nf_put_att_text(ncid, jul_id, 'long_name', 15, 'true julian day')
   iret = nf_put_att_text(ncid, jul_id, 'units', 4, 'days')

   iret = nf_put_att_text(ncid, secs_id, 'long_name', 22, 'seconds since midnight')
   iret = nf_put_att_text(ncid, secs_id, 'units', 4, 'seconds')

   iret = nf_put_att_text(ncid, elev_id, 'long_name', 16, 'tidal elevations')
!kbk     call check_err(iret)
   iret = nf_put_att_text(ncid, elev_id, 'units', 6, 'meters')
!kbk    call check_err(iret)
   realval(1) = -15
   iret = nf_put_att_real(ncid, elev_id, 'valid_min', NF_REAL, 1, realval)
!kbk    call check_err(iret)
   realval(1) = 15
   iret = nf_put_att_real(ncid, elev_id, 'valid_max', NF_REAL, 1, realval)
!kbk    call check_err(iret)
   realval(1) = -99
   iret = nf_put_att_real(ncid, elev_id, 'missing_value', NF_REAL, 1,realval)
!kbk    call check_err(iret)

   iret = nf_enddef(ncid)
!kbk    call check_err(iret)

   return

end
