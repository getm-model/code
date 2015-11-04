#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: register_all_variables() - register GETM variables.
!
! !INTERFACE:
   subroutine register_all_variables(runtype,fm)
!
! !DESCRIPTION:
!
! !USES:
   use field_manager
   use domain
   use meteo
   use variables_2d
#ifndef NO_3D
   use variables_3d
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)               :: runtype
!
! !INPUT/OUTPUT PARAMETERS:
!KB   class (type_field_manager), target intent(inout) :: fm
   type (type_field_manager), target :: fm
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer,parameter :: rk = kind(_ONE_)
   character(len=16)         :: xname=''
   character(len=16)         :: xlongname=''
   character(len=16)         :: xunits=''
   character(len=16)         :: yname=''
   character(len=16)         :: ylongname=''
   character(len=16)         :: yunits=''
   character(len=16)         :: zname='sigma'
   character(len=64)         :: zlongname='sigma'
   character(len=16)         :: zunits='sigma'
!EOP
!-----------------------------------------------------------------------
!BOC

   LEVEL1 'register_all_variables()'

   select case (grid_type)
      case (1)
         xname     = 'xc'
         xlongname = 'x'
         xunits    = 'm'
         yname     = 'yc'
         ylongname = 'y'
         yunits    = 'm'
      case (2)
         xname     = 'lonc'
         xlongname = 'longitude'
         xunits    = 'degrees_east'
         yname     = 'latc'
         ylongname = 'latitude'
         yunits    = 'degrees_north'
      case (3)
         xname     = 'xic'
         xlongname = 'xic'
         yname     = 'etac'
         ylongname = 'etac'
   end select

#ifndef NO_3D
   if (runtype .ge. 2) then
      select case (vert_cord)
         case (1)
            zname     = 'sigma'
            zlongname = 'sigma layers'
            zunits    = 'sigma_level'
         case (2)
            zname     = 'z'
            zlongname = 'geopotential'
            zunits    = 'm'
         case (3,4,5)
            zname  = 'level'
            zlongname  = 'general vertical coordinates'
            zunits = 'level'
         case default
      end select
   end if
#endif


!  register - dimensions
   call fm%register_dimension(trim(xname),imax-imin+1,global_length=iextr,offset=ioff,id=id_dim_lon)
   call fm%register_dimension(trim(yname),jmax-jmin+1,global_length=jextr,offset=joff,id=id_dim_lat)
#ifndef NO_3D
   if (runtype .ge. 2) then
      call fm%register_dimension(trim(zname),kmax+1,id=id_dim_z)
   end if
#endif

   call fm%register_dimension('time',id=id_dim_time)
   call fm%initialize(prepend_by_default=(/id_dim_lon,id_dim_lat/),append_by_default=(/id_dim_time/))


!  register - domain
   call fm%register(trim(xname),trim(xunits),trim(xlongname),dimensions=(/id_dim_lon/),no_default_dimensions=.true.,data1d=xcord(_IRANGE_NO_HALO_),coordinate_dimension=id_dim_lon)
   call fm%register(trim(yname),trim(yunits),trim(ylongname),dimensions=(/id_dim_lat/),no_default_dimensions=.true.,data1d=ycord(_JRANGE_NO_HALO_),coordinate_dimension=id_dim_lat)
#ifndef NO_3D
   if (runtype .ge. 2) then
      call fm%register(trim(zname),trim(zunits),trim(zlongname),dimensions=(/id_dim_z/),no_default_dimensions=.true.,data1d=ga,coordinate_dimension=id_dim_z)
   end if
#endif

   call fm%register('bathymetry', 'm', 'bathymetry', standard_name='bathymetry', dimensions=(/id_dim_lon,id_dim_lat/), no_default_dimensions=.true., fill_value=-10._rk, data2d=H(_2D_W_), category='domain',output_level=output_level_required)

!  register -  metric
   call fm%register('dxc', 'm', 'dx at T-points', dimensions=(/id_dim_lon,id_dim_lat/), no_default_dimensions=.true., data2d=dxc(_2D_W_), category="metrics", output_level=output_level_debug)
   call fm%register('dyc', 'm', 'dy at T-points', dimensions=(/id_dim_lon,id_dim_lat/), no_default_dimensions=.true., data2d=dyc(_2D_W_), category="metrics", output_level=output_level_debug)


!  category - 2d
   call fm%register('z', 'm', 'sea surface elevation', standard_name='sea surface elevation', fill_value=10.05_rk, data2d=z(_2D_W_), category="2d")

   call fm%register('zo', 'm', 'sea surface elevation', standard_name='sea surface elevation', fill_value=10.05_rk, data2d=zo(_2D_W_), category="2d", output_level=output_level_debug)
   call fm%register('D', 'm', 'water depth', standard_name='water depth', fill_value=10.05_rk, data2d=D(_2D_W_), category="2d")

#ifndef NO_3D
!  category - 3d
   if (runtype .ge. 2) then
      call fm%register('hn', 'm', 'layer thickness', standard_name='cell_thickness', dimensions=(/id_dim_z/),data3d=hn(_3D_W_), category='grid')
      call fm%register('hun', 'm', 'layer thickness - U-points', standard_name='cell_thickness', dimensions=(/id_dim_z/),data3d=hun(_3D_W_), category='grid', output_level=output_level_debug)
      call fm%register('hvn', 'm', 'layer thickness - V-points', standard_name='cell_thickness', dimensions=(/id_dim_z/),data3d=hvn(_3D_W_), category='grid', output_level=output_level_debug)
      call fm%register('ssen', 'm', 'elevtion at T-points (3D)', standard_name='', data2d=ssen(_2D_W_), category='3d', output_level=output_level_debug)
      call fm%register('ssun', 'm', 'elevtion at U-points (3D)', standard_name='', data2d=ssun(_2D_W_), category='3d', output_level=output_level_debug)
      call fm%register('ssvn', 'm', 'elevtion at V-points (3D)', standard_name='', data2d=ssvn(_2D_W_), category='3d', output_level=output_level_debug)
   end if
#ifndef NO_BAROCLINIC
   if (runtype .ge. 3) then
      call fm%register('temp', 'Celsius', 'temperature', standard_name='', dimensions=(/id_dim_z/),data3d=T(_3D_W_), category='baroclinic')
      call fm%register('salt', 'PSU', 'salinity', standard_name='', dimensions=(/id_dim_z/),data3d=S(_3D_W_), category='baroclinic')
      call fm%register('idpdx', 'm', 'baroclinic pressure gradient - x', standard_name='', dimensions=(/id_dim_z/),data3d=idpdx(_3D_W_), category='baroclinic', output_level=output_level_debug)
      call fm%register('idpdy', 'm', 'baroclinic pressure gradient - y', standard_name='', dimensions=(/id_dim_z/),data3d=idpdy(_3D_W_), category='baroclinic', output_level=output_level_debug)
   end if
#endif
#endif

   return
   end subroutine register_all_variables
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2015 - Karsten Bolding and Jorn Bruggeman(BB)          !
!-----------------------------------------------------------------------
