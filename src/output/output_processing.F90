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
   use field_manager
   IMPLICIT NONE
!
! !PUBLIC DATA FUNCTIONS:
   public init_output_processing, do_output_processing
!
! !PUBLIC DATA MEMBERS:
   REALTYPE, dimension(:,:),   allocatable, target :: u2d, v2d
   REALTYPE, dimension(:,:),   allocatable, target :: u2d_destag, v2d_destag
   REALTYPE, dimension(:,:,:), allocatable, target :: u3d, v3d
   REALTYPE, dimension(:,:,:), allocatable, target :: u3d_destag, v3d_destag
!
! !PRIVATE DATA MEMBERS:
   logical, target:: u2d_use, v2d_use
   logical, target:: u2d_destag_use, v2d_destag_use
   logical, target:: u3d_use, v3d_use
   logical, target:: u3d_destag_use, v3d_destag_use
   integer, parameter :: rk = kind(_ONE_)
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
   allocate(u2d(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_output_processing: Error allocating memory (u2d)'
   u2d = 0._rk
   allocate(v2d(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_output_processing: Error allocating memory (v2d)'
   v2d = 0._rk
   allocate(u2d_destag(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_output_processing: Error allocating memory (u2d_destag)'
   u2d_destag = 0._rk
   allocate(v2d_destag(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_output_processing: Error allocating memory (v2d_destag)'
   v2d_destag = 0._rk

#if 0
   allocate(u3d(I3DFIELD),stat=rc)
   if (rc /= 0) stop 'init_output_processing: Error allocating memory (u3d)'
   allocate(v3d(I3DFIELD),stat=rc)
   if (rc /= 0) stop 'init_output_processing: Error allocating memory (v3d)'
   allocate(u3d_destag(I3DFIELD),stat=rc)
   if (rc /= 0) stop 'init_output_processing: Error allocating memory (u3d_destag)'
   allocate(v3d_destag(I3DFIELD),stat=rc)
   if (rc /= 0) stop 'init_output_processing: Error allocating memory (v3d_destag)'
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
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_field_manager) :: fm
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'register_processed_variables()'

   call fm%register('u2d', 'm/s', 'velocity in local x-direction', standard_name='', data2d=u2d(_2D_W_), fill_value=-9999._rk, category='velocities', used_now=u2d_use)
   call fm%register('v2d', 'm/s', 'velocity in local y-direction', standard_name='', data2d=v2d(_2D_W_), fill_value=-9999._rk, category='velocities', used_now=v2d_use)
   call fm%register('u2d-destag', 'm/s', 'velocity in local x-direction(destag)', standard_name='', data2d=u2d_destag(_2D_W_), fill_value=-9999._rk, category='velocities',output_level=output_level_debug, used_now=u2d_destag_use)
   call fm%register('v2d-destag', 'm/s', 'velocity in local y-direction(destag)', standard_name='', data2d=v2d_destag(_2D_W_), fill_value=-9999._rk, category='velocities',output_level=output_level_debug, used_now=v2d_destag_use)

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
   use domain, only: az, au, av
   use variables_2d, only: z,D
   use variables_2d, only: U,V,DU,DV
   use variables_3d, only: kmin,hn,uu,hun,vv,hvn
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE, parameter                 :: vel_missing        =-9999.
!EOP
!-------------------------------------------------------------------------
!BOC

!  2D - velocities
   if (u2d_use .and. v2d_use) then
      call to_2d_vel(imin,jmin,imax,jmax,au,U,DU,vel_missing,       &
                     imin,jmin,imax,jmax,u2d)
      call to_2d_vel(imin,jmin,imax,jmax,av,V,DV,vel_missing,       &
                     imin,jmin,imax,jmax,v2d)
   end if
   if (u2d_destag_use .and. v2d_destag_use) then
      call to_2d_u(imin,jmin,imax,jmax,az,U,DU,vel_missing,      &
                   imin,jmin,imax,jmax,u2d_destag)
      call to_2d_v(imin,jmin,imax,jmax,az,V,DV,vel_missing,      &
                   imin,jmin,imax,jmax,v2d_destag)
   end if

#if 0
!  3D - velocities
#ifndef NO_3D
   if (allocated(u3d) .and. allocated(v3d)) then
      call to_3d_uu(imin,jmin,imax,jmax,kmin,kmax,az, &
                    hun,uu,vel_missing,u3d)
      call to_3d_vv (imin,jmin,imax,jmax,kmin,kmax,az, &
                     hvn,vv,vel_missing,v3d)
   end if

   if (allocated(u3d_destag) .and. allocated(v3d_destag)) then
      call to_3d_vel(imin,jmin,imax,jmax,kmin,kmax,au, &
                     hun,uu,vel_missing,u3d_destag)
      call to_3d_vel(imin,jmin,imax,jmax,kmin,kmax,av, &
                     hvn,vv,vel_missing,v3d_destag)
   end if
#endif
#endif

   return
   end subroutine do_output_processing
!EOC

   end module output_processing

!-----------------------------------------------------------------------
! Copyright (C) 2019 - Karsten Bolding & Jorn Bruggeman (BB)           !
!-----------------------------------------------------------------------
