!$Id: ncdf_out.F90,v 1.1 2002-05-02 14:01:52 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdf_out - interfaces for NetCDF IO subroutines
!
! !INTERFACE:
   module ncdf_out
!
! !DESCRIPTION:
!  This file contains the interface for the files in the ../ncdf/ directory
!  that does the actual saving in the NetCDF format. To avoid circular
!  dependencies it has to be here instead of - more naturally - in the
!  ncdf-directory.
!
! !USE:
   IMPLICIT NONE
!
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_2d_ncdf,save_2d_ncdf
   public init_3d_ncdf,save_3d_ncdf
!
! !PUBLIC DATA MEMBERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ncdf_out.F90,v $
!  Revision 1.1  2002-05-02 14:01:52  gotm
!  Initial revision
!
!  Revision 1.2  2001/10/17 14:43:52  bbh
!  save_meteo not passed as argument
!
!  Revision 1.1  2001/09/13 14:50:02  bbh
!  Cleaner and smaller NetCDF implementation + better axis support
!
!
!EOP
!-----------------------------------------------------------------------

   interface
      subroutine init_2d_ncdf(fn,title,starttime)
         character(len=*), intent(in)	:: fn,title,starttime
      end subroutine init_2d_ncdf
   end interface

   interface
      subroutine save_2d_ncdf(secs)
         REALTYPE, intent(in)		:: secs
      end subroutine save_2d_ncdf
   end interface

   interface
      subroutine init_3d_ncdf(fn,title,starttime)
         character(len=*), intent(in)	:: fn,title,starttime
      end subroutine init_3d_ncdf
   end interface

   interface
      subroutine save_3d_ncdf(secs)
         REALTYPE, intent(in)		:: secs
      end subroutine save_3d_ncdf
   end interface

!-----------------------------------------------------------------------

   end module ncdf_out

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Karsten Bolding and Hans Burchard (BBH)         !
!-----------------------------------------------------------------------
