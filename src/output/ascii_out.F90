!$Id: ascii_out.F90,v 1.1 2002-05-02 14:01:52 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ascii_out - interfaces for NetCDF IO subroutines
!
! !INTERFACE:
   module ascii_out
!
! !DESCRIPTION:
!  This file contains the interface for the files in the ../ascii/ directory
!  that does the actual saving in the NetCDF format. To avoid circular
!  dependencies it has to be here instead of - more naturally - in the
!  ascii-directory.
!
! !USE:
   IMPLICIT NONE
!
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_2d_ascii,save_2d_ascii
   public init_3d_ascii,save_3d_ascii
!
! !PUBLIC DATA MEMBERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ascii_out.F90,v $
!  Revision 1.1  2002-05-02 14:01:52  gotm
!  Initial revision
!
!  Revision 1.1  2001/09/13 14:50:34  bbh
!  Stubs for ascii output
!
!
!EOP
!-----------------------------------------------------------------------

   interface
      subroutine init_2d_ascii(fn,title,starttime)
         character(len=*), intent(in)	:: fn,title,starttime
      end subroutine init_2d_ascii
   end interface

   interface
      subroutine save_2d_ascii(secs,save_meteo)
         REALTYPE, intent(in)		:: secs
         logical, intent(in)		:: save_meteo
      end subroutine save_2d_ascii
   end interface

   interface
      subroutine init_3d_ascii(fn,title,starttime)
         character(len=*), intent(in)	:: fn,title,starttime
      end subroutine init_3d_ascii
   end interface

   interface
      subroutine save_3d_ascii(secs)
         REALTYPE, intent(in)		:: secs
      end subroutine save_3d_ascii
   end interface

!-----------------------------------------------------------------------

   end module ascii_out

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Karsten Bolding and Hans Burchard (BBH)         !
!-----------------------------------------------------------------------
