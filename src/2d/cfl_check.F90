!$Id: cfl_check.F90,v 1.4 2003-04-23 12:09:43 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: cfl_check() - barotropic residual currents.
!
! !INTERFACE:
   subroutine cfl_check()
!
! !DESCRIPTION:
!
! !USES:
   use parameters, only: g
   use domain, only: imin,imax,jmin,jmax,H,az
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dyc,dxc
#else
   use domain, only: dy,dx
#endif
   use m2d, only: dtm
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: cfl_check.F90,v $
!  Revision 1.4  2003-04-23 12:09:43  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.3  2003/04/03 07:01:49  gotm
!  fixed CFL calc. for non cartesian grid
!
!  Revision 1.2  2002/10/04 13:56:58  gotm
!  Uses Becker and Deleersnijder (1993) CFL criterion with Coriolis
!
!  Revision 1.1.1.1  2002/05/02 14:00:41  gotm
!  recovering after CVS crash
!
!  Revision 1.7  2001/09/19 11:25:14  gotm
!  Removed E2DFIELD when de-allocating lmask
!
!  Revision 1.6  2001/09/19 11:20:32  bbh
!  Explicit de-allocates memory when -DFORTRAN90
!
!  Revision 1.5  2001/09/01 17:07:58  bbh
!  Removed some print statements
!
!  Revision 1.4  2001/08/01 08:25:52  bbh
!  CURVILINEAR now implemented
!
!  Revision 1.3  2001/04/20 14:03:39  bbh
!  should not have deleted imin,imax,jmin,jmax :-)
!
!  Revision 1.2  2001/04/20 13:47:33  bbh
!  Fixed bug concerning optional argument
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
   integer                   :: pos(2),max_pos(2),rc,i,j
   REALTYPE                  :: h_max=-99.,c,max_dt,dtt
   logical, dimension(:,:), allocatable :: lmask
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'CFL check'

   allocate(lmask(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'cfl_check: Error allocating memory (lmask)'

   lmask = .false.
   lmask = (az .gt. 0)
   h_max = maxval(H,mask = lmask)
   pos = maxloc(H,mask = lmask)
#if defined(SPHERICAL) || defined(CURVILINEAR)
   max_dt=1000000000.
   do i=imin,imax
      do j=jmin,jmax
         if (az(i,j) .ge. 1 .and. H(i,j) .gt. _ZERO_) then
#if 0
            dtt=min(dxc(i,j),dyc(i,j))/sqrt(2.*g*H(i,j))
#else
            c = sqrt(g*H(i,j))
            dtt = (dxc(i,j)*dyc(i,j))/ &
                   (sqrt(2.0)*c*sqrt(dxc(i,j)*dxc(i,j)+dyc(i,j)*dyc(i,j)))
 
#endif
            if (dtt .lt. max_dt) then
               max_dt=dtt
               max_pos(1)=i
               max_pos(2)=j
            end if
         end if
      end do
   end do
   LEVEL3 'max CFL number at depth=',H(max_pos(1),max_pos(2)),' at ',max_pos
   LEVEL3 'at this position, dx = ',dxc(max_pos(1),max_pos(2)),' and dy =  ',dyc(max_pos(1),max_pos(2))
#else
   c = sqrt(g*h_max)
!  Becker and Deleersnijder
   max_dt = (dx*dy)/(sqrt(2.0)*c*sqrt(dx*dx+dy*dy)) 

#if 0
#ifdef POM_CFL
   max_dt=0.5/(c*sqrt( _ONE_/dx**2 + _ONE_/dy**2 ))
#endif
#ifdef OLD_GETM
   max_dt = min(dx,dy)/(c*sqrt(2.))
#endif
#endif

   LEVEL3 'max depth=',h_max,' at ',pos
#endif

   if (dtm .gt. max_dt) then
      FATAL 'reduce time-step (',dtm,') below ',max_dt
      stop 'cfl_check()'
   else
      LEVEL3 'used dt (',dtm,') is less than ',max_dt
   end if

#ifdef FORTRAN90
   deallocate(lmask,stat=rc)
   if (rc /= 0) stop 'cfl_check: Error allocating memory (lmask)'
#endif
   return
   end subroutine cfl_check
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
