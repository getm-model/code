#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: to_u() - calculates cell centered velocity in global x-direction
!
! !INTERFACE:
   subroutine to_u(imin,jmin,imax,jmax,az,                            &
                   dt,grid_type,                                      &
#if defined(CURVILINEAR) || defined(SPHERICAL)
                   dxv,dyu,arcd1,                                     &
#else
                   dx,dy,ard1,                                        &
#endif
                   xc,xu,xv,D,Dlast,U,DU,V,DV,wwm,wwp,missing,velx)
!
! !DESCRIPTION:
!
! Here the lateral velocity in the x-y-system at cell centres
! (but velocity time stages) is calculated.
!
! !USES:
!$ use omp_lib
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)                       :: imin,jmin,imax,jmax,grid_type
   integer,dimension(E2DFIELD),intent(in)   :: az
   REALTYPE,intent(in)                      :: dt
#if defined(CURVILINEAR) || defined(SPHERICAL)
   REALTYPE,dimension(E2DFIELD),intent(in)  :: dxv,dyu,arcd1
#else
   REALTYPE,intent(in)                      :: dx,dy,ard1
#endif
   REALTYPE,dimension(E2DFIELD),intent(in)  :: xc,xu,xv,D,Dlast,U,DU,V,DV,wwm,wwp
   REALTYPE,intent(in)                      :: missing
!
! !OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(out) :: velx
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
! !LOCAL VARIABLES:
   integer                                  :: i,j
   REALTYPE                                 :: dtm1
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'to_u() # ',Ncall
#endif
#ifdef SLICE_MODEL
!  Note (KK): this value MUST NOT be changed !!!
   j = jmax/2
#endif

   dtm1 = _ONE_/dt

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j)                                         &
!$OMP          PRIVATE(i)


   select case(grid_type)

      case (1,2)

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
         do j=jmin-HALO,jmax+HALO
#endif
            do i=imin-HALO+1,imax+HALO
               if (az(i,j) .eq. 1) then
                  velx(i,j) = ( U(i-1,j) + U(i,j) ) / ( Dlast(i,j) + D(i,j) )
               else
                  velx(i,j) = missing
               end if
            end do
#ifndef SLICE_MODEL
         end do
#endif
!$OMP END DO

      case (3)

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
         do j=jmin-HALO+1,jmax+HALO
#endif
            do i=imin-HALO+1,imax+HALO
               if (az(i,j) .eq. 1) then
                  velx(i,j) = (                                   &
                                 (                                &
                                     ( D(i,j) - Dlast(i,j) )*dtm1  &
                                   + ( wwp(i,j) - wwm(i,j) )      &
                                 )                                &
                                 *xc(i,j)                         &
                               + (                                &
                                    U(i  ,j  )*xu(i  ,j  )*DYU    &
                                  - U(i-1,j  )*xu(i-1,j  )*DYUIM1 &
#ifndef SLICE_MODEL
                                  + V(i  ,j  )*xv(i  ,j  )*DXV    &
                                  - V(i  ,j-1)*xv(i  ,j-1)*DXVJM1 &
#endif
                                 )                                &
                                 *ARCD1                           &
                              )                                   &
                              /(_HALF_*(Dlast(i,j)+D(i,j)))
               else
                  velx(i,j) = missing
               end if
            end do
#ifndef SLICE_MODEL
         end do
#endif
!$OMP END DO

      case (4)
         stop 'to_u: grid_type=4 not implemented yet'

      case default
         stop 'to_u: invalid grid_type'

   end select

!$OMP END PARALLEL

#ifdef SLICE_MODEL
   velx(:,j+1) = velx(:,j)
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving to_u()'
   write(debug,*)
#endif
   return
   end subroutine to_u
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2012 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
