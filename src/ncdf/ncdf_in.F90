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
   integer, public                     :: ncbathy
   integer, public                     :: h_id,dx_id,dy_id
   integer, public                     :: dims,dimids(2),start(2),edges(2)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-----------------------------------------------------------------------

   interface
#ifdef STATIC
      subroutine get_dimensions(fname,rc)
         character(len=*), intent(in)  :: fname
         integer,intent(out)           :: rc
#else
      subroutine get_dimensions(fname,iextr,jextr,rc)
         character(len=*), intent(in)  :: fname
         integer,intent(out)           :: iextr,jextr,rc
#endif
      end subroutine get_dimensions
   end interface
#if 0
   interface
      subroutine get_bathymetry(H,Hland,dx,dy,imin,imax,jmin,jmax,rc)
         integer, intent(inout)        :: imin,imax,jmin,jmax
         REALTYPE, intent(out)         :: H(E2DFIELD)
         REALTYPE, intent(out)         :: Hland,dx,dy
         integer, intent(out)          :: rc
      end subroutine get_bathymetry
   end interface
#endif

   interface
      subroutine get_field_ncdf(fname,var,f)
         character(len=*), intent(in)  :: fname,var
         REALTYPE, intent(out)         :: f
      end subroutine get_field_ncdf
   end interface

!-----------------------------------------------------------------------

   end module ncdfin

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
