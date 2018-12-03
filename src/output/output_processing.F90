#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: output_processing
!
! !INTERFACE:
   module output_processing
!
! !DESCRIPTION:
!  This modules serves as a container for processing output variables. 

! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax
   use domain, only: az, au, av
   use processed_variables
   IMPLICIT NONE
!
! !PUBLIC DATA FUNCTIONS:
   public init_output_processing, do_output_processing
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_output_processing - read required variables
!
! !INTERFACE:
   subroutine init_output_processing
!
! !USES:
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
!    character(len=*), intent(in)        :: filename
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                             :: rc
!EOP
!-------------------------------------------------------------------------
!BOC
   allocate(u_2d(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_output_processing: Error allocating memory (u_2d)'
   allocate(v_2d(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_output_processing: Error allocating memory (v_2d)'
   allocate(u_2d_destag(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_output_processing: Error allocating memory (u_2d_destag)'
   allocate(v_2d_destag(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_output_processing: Error allocating memory (v_2d_destag)'

#if 0
   allocate(u_3d(I3DFIELD),stat=rc)
   if (rc /= 0) stop 'init_output_processing: Error allocating memory (u_3d)'
   allocate(v_3d(I3DFIELD),stat=rc)
   if (rc /= 0) stop 'init_output_processing: Error allocating memory (v_3d)'
   allocate(u_3d_destag(I3DFIELD),stat=rc)
   if (rc /= 0) stop 'init_output_processing: Error allocating memory (u_3d_destag)'
   allocate(v_3d_destag(I3DFIELD),stat=rc)
   if (rc /= 0) stop 'init_output_processing: Error allocating memory (v_3d_destag)'
#endif
   return
   end subroutine init_output_processing
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: register_processed_variables()
!
! !INTERFACE:
   subroutine register_processed_variables(fm)
!
! !DESCRIPTION:
!
! !USES:
!   use variables_2d
   use field_manager
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_field_manager) :: fm
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
! !LOCAL VARIABLES:
   logical :: used
   integer,parameter :: rk = kind(_ONE_)
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'register_processed_variables()'

   call fm%register('u2d', 'm/s', 'velocity in local x-direction', standard_name='', data2d=u_2d(_2D_W_), category='velocities')
   call fm%register('v2d', 'm/s', 'velocity in local y-direction', standard_name='', data2d=v_2d(_2D_W_), category='velocities')
   call fm%register('u2d_destag', 'm/s', 'velocity in local x-direction(destag)', standard_name='', data2d=u_2d_destag(_2D_W_), category='velocities',output_level=output_level_debug)
   call fm%register('v2d_destag', 'm/s', 'velocity in local y-direction(destag)', standard_name='', data2d=v_2d_destag(_2D_W_), category='velocities',output_level=output_level_debug)

   return
   end subroutine register_processed_variables
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: do_output_processing - read required variables
!
! !INTERFACE:
  subroutine do_output_processing
!
! !USES:
   use variables_2d, only: z,D
   use variables_2d, only: U,V,DU,DV
   use variables_3d, only: kmin,hn,uu,hun,vv,hvn
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
!    character(len=*), intent(in)        :: filename
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL VARIABLES:
!   integer                             :: ncid
!KB - maybe this should go
   REALTYPE, parameter                 :: vel_missing        =-9999.
!EOP
!-------------------------------------------------------------------------
!BOC
   if (allocated(u_2d) .and. allocated(v_2d)) then
      call to_2d_vel(imin,jmin,imax,jmax,au,U,DU,vel_missing,       &
                     imin,jmin,imax,jmax,u_2d)
      call to_2d_vel(imin,jmin,imax,jmax,av,V,DV,vel_missing,       &
                     imin,jmin,imax,jmax,v_2d)
   end if
   if (allocated(u_2d_destag) .and. allocated(v_2d_destag)) then
      call to_2d_u(imin,jmin,imax,jmax,az,U,DU,vel_missing,      &
                   imin,jmin,imax,jmax,u_2d_destag)
      call to_2d_v(imin,jmin,imax,jmax,az,V,DV,vel_missing,      &
                   imin,jmin,imax,jmax,v_2d_destag)
   end if

#ifndef NO_3D
   if (allocated(u_3d) .and. allocated(v_3d)) then
      call to_3d_uu(imin,jmin,imax,jmax,kmin,kmax,az, &
                    hun,uu,vel_missing,u_3d)
      call to_3d_vv (imin,jmin,imax,jmax,kmin,kmax,az, &
                     hvn,vv,vel_missing,v_3d)
   end if

   if (allocated(u_3d_destag) .and. allocated(v_3d_destag)) then
      call to_3d_vel(imin,jmin,imax,jmax,kmin,kmax,au, &
                     hun,uu,vel_missing,u_3d_destag)
      call to_3d_vel(imin,jmin,imax,jmax,kmin,kmax,av, &
                     hvn,vv,vel_missing,v_3d_destag)
   end if
#endif

   return
   end subroutine do_output_processing
!EOC

   end module output_processing

!-----------------------------------------------------------------------
! Copyright (C) 2019 - Karsten Bolding & Jorn Bruggeman (BB)           !
!-----------------------------------------------------------------------
