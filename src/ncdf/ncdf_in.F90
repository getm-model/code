!$Id: ncdf_in.F90,v 1.2 2003-04-07 15:36:08 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  ncdfin - input in NetCDF format
!
! !INTERFACE:
   module ncdfin
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
   private
!
   public :: get_field_ncdf
!
! !PUBLIC DATA MEMBERS:
   integer, public	:: ncbathy
   integer, public	:: h_id,dx_id,dy_id
   integer, public	:: dims,dimids(2),start(2),edges(2)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ncdf_in.F90,v $
!  Revision 1.2  2003-04-07 15:36:08  kbk
!  not including interface for get_bathymetry
!
!  Revision 1.1.1.1  2002/05/02 14:01:47  gotm
!  recovering after CVS crash
!
!  Revision 1.6  2001/05/25 19:26:22  bbh
!  ncdf_meteo.F90
!
!  Revision 1.5  2001/05/14 12:49:13  bbh
!  Removed init_2d_bdy_ncdf and get_2d_bdy_ncdf
!
!  Revision 1.4  2001/05/10 11:38:29  bbh
!  Added get_field_ncdf() + various small bug fixes
!
!  Revision 1.3  2001/05/07 14:43:01  bbh
!  Removed arrays julday and secs - not used anyway
!
!  Revision 1.2  2001/05/06 18:51:55  bbh
!  Towards proper implementation of specified 2D bdy.
!
!  Revision 1.1.1.1  2001/04/17 08:43:07  bbh
!  initial import into CVS
!
!EOP
!-----------------------------------------------------------------------

   interface
#ifdef STATIC
      subroutine get_dimensions(fname,rc)
         character(len=*), intent(in)	:: fname
         integer,intent(out)		:: rc
#else
      subroutine get_dimensions(fname,iextr,jextr,rc)
         character(len=*), intent(in)	:: fname
         integer,intent(out)		:: iextr,jextr,rc
#endif
      end subroutine get_dimensions
   end interface
#if 0
   interface
      subroutine get_bathymetry(H,Hland,dx,dy,imin,imax,jmin,jmax,rc)
         integer, intent(inout)		:: imin,imax,jmin,jmax
         REALTYPE, intent(out)		:: H(E2DFIELD)
         REALTYPE, intent(out)		:: Hland,dx,dy
         integer, intent(out)		:: rc
      end subroutine get_bathymetry
   end interface
#endif

   interface
      subroutine get_field_ncdf(fname,var,f)
         character(len=*), intent(in)	:: fname,var
         REALTYPE, intent(out)		:: f
      end subroutine get_field_ncdf
   end interface

!-----------------------------------------------------------------------

   end module ncdfin

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
