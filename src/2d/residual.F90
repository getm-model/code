!$Id: residual.F90,v 1.4 2004-01-03 16:40:28 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: do_residual() - barotropic residual currents.
!
! !INTERFACE:
   subroutine do_residual(finish)
!
! !DESCRIPTION:
!
! !USES:
   use variables_2d, only: u,v,res_du,res_u,res_dv,res_v,du,dv
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: finish
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: residual.F90,v $
!  Revision 1.4  2004-01-03 16:40:28  kbk
!  bug fix
!
!  Revision 1.3  2003/04/23 12:09:44  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 15:30:58  kbk
!  needs to be fixed
!
!  Revision 1.1.1.1  2002/05/02 14:00:45  gotm
!  recovering after CVS crash
!
!  Revision 1.2  2001/05/03 19:35:01  bbh
!  Use of variables_2d
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   if(finish .eq. 1) then
      LEVEL2 'Calculating residual currents'
      where ( res_du .ne. _ZERO_ )
         res_u=res_u/res_du
      elsewhere
         res_u= _ZERO_
      end where
      where ( res_dv .ne. _ZERO_ )
         res_v=res_v/res_dv
      elsewhere
         res_v= _ZERO_
      end where
   else
      res_du=res_du+du
      res_u=res_u+u
      res_dv=res_dv+dv
      res_v=res_v+v
   end if

   return
   end subroutine do_residual
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
