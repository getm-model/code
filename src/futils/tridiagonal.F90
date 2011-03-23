#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: getm_tridiagonal() - solves a linear system of equations.
!
! !INTERFACE:
   subroutine getm_tridiagonal(kmax,fi,lt,au,bu,cu,du,value)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT none
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: kmax
   integer, intent(in)                 :: fi,lt
   REALTYPE, intent(in)                :: au(0:kmax),bu(0:kmax)
   REALTYPE, intent(in)                :: cu(0:kmax),du(0:kmax)
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: value(0:kmax)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   REALTYPE                  :: ru(0:kmax),qu(0:kmax)
   integer                   :: i
!EOP
!-----------------------------------------------------------------------
!BOC
! OMP-NOTE: As all local storage (i, ru, and qu) is on stack, this routine
!  should be thread-safe. BJB 2009-09-21.
   ru(lt)=au(lt)/bu(lt)
   qu(lt)=du(lt)/bu(lt)

   do i=lt-1,fi+1,-1
      ru(i)=au(i)/(bu(i)-cu(i)*ru(i+1))
      qu(i)=(du(i)-cu(i)*qu(i+1))/(bu(i)-cu(i)*ru(i+1))
   end do

   qu(fi)=(du(fi)-cu(fi)*qu(fi+1))/(bu(fi)-cu(fi)*ru(fi+1))

   value(fi)=qu(fi)
   do i=fi+1,lt
      value(i)=qu(i)-ru(i)*value(i-1)
   end do

   return
   end subroutine getm_tridiagonal
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
