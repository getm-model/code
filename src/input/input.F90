!$Id: input.F90,v 1.1.1.1 2002-05-02 14:01:33 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  input - input specifications
!
! !INTERFACE:
   module input
!
! !DESCRIPTION:
!
! !USES:
   use meteo, only: metforcing,method,meteo_file
   use rivers, only: river_method,nriver,river_data
   use m2d, only: bdy2d,bdyfile_2d,bdyfmt_2d
   use m3d, only: bdy3d,bdyfile_3d,bdyfmt_3d
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: input.F90,v $
!  Revision 1.1.1.1  2002-05-02 14:01:33  gotm
!  recovering after CVS crash
!
!  Revision 1.5  2001/10/07 14:50:22  bbh
!  Reading river data implemented - NetCFD
!
!  Revision 1.4  2001/07/26 13:57:14  bbh
!  Meteo working - needs some polishing
!
!  Revision 1.3  2001/05/25 19:03:02  bbh
!  New method for - input - all is done via do_input()
!
!  Revision 1.2  2001/05/10 11:33:48  bbh
!  Added wrapper - get_field
!
!  Revision 1.1.1.1  2001/04/17 08:43:09  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------

   interface
      subroutine init_2d_bdy(fn,fmt)
         character(len=*), intent(in)	:: fn
         integer, intent(in)	:: fmt
      end subroutine init_2d_bdy
   end interface

   interface
      subroutine get_2d_bdy(fmt,n)
         integer, intent(in)	:: fmt,n
      end subroutine get_2d_bdy
   end interface

   interface
      subroutine init_3d_bdy(fn,fmt)
         character(len=*), intent(in) :: fn
         integer, intent(in)	:: fmt
      end subroutine init_3d_bdy
   end interface

   interface
      subroutine get_3d_bdy(fmt,n)
         integer, intent(in)	:: fmt,n
      end subroutine get_3d_bdy
   end interface

   interface
      subroutine get_field(fname,var,f)
	 character(len=*), intent(in) :: fname,var
	 REALTYPE, intent(out)	:: f
      end subroutine get_field
   end interface

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_input - initialise all external files and units
!
! !INTERFACE:
   subroutine init_input(n)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)	:: n
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_input() # ',Ncall
#endif

   LEVEL1 'init_input'

   if (metforcing .and. method .eq. 2) then
      call init_meteo_input(meteo_file,n)
   end if

   if (river_method .gt. 0 .and. nriver .gt. 0) then
      call init_river_input(river_data,n)
   end if

   if(bdy2d) then
      call init_2d_bdy(bdyfile_2d,bdyfmt_2d)
   end if

   if(bdy3d) then
      call init_3d_bdy(bdyfile_3d,bdyfmt_3d)
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving init_input()'
   write(debug,*)
#endif
   return
   end subroutine init_input
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_input - cleans up after run
!
! !INTERFACE:
   subroutine do_input(n)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  To be written
!
! !INPUT PARAMETERS:
   integer, intent(in) :: n
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_input() # ',Ncall
#endif

   if(metforcing .and. method .eq. 2) then
      call get_meteo_data(n)
   end if

   if(river_method .eq. 2) then
      call get_river_data(n)
   end if

   if(bdy2d) then
      call get_2d_bdy(bdyfmt_2d,n)
   end if

   if(bdy3d) then
      call get_3d_bdy(bdyfmt_3d,n)
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving do_input()'
   write(debug,*)
#endif
   return
   end subroutine do_input
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clean_input - cleans up after run
!
! !INTERFACE:
   subroutine clean_input()
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Writes calculated fields to files.
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
! !LOCAL VARIABLES:
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'clean_input() # ',Ncall
#endif

   LEVEL1 'clean_input'

#ifdef DEBUG
   write(debug,*) 'Leaving clean_io()'
   write(debug,*)
#endif
   return
   end subroutine clean_input
!EOC

!-----------------------------------------------------------------------

   end module input

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
