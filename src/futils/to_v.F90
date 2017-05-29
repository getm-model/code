#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: to_v() - calculates cell centered velocity in global y-direction
!
! !INTERFACE:
   subroutine to_v(imin,jmin,imax,jmax,az,                               &
                   dt,grid_type,                                         &
                   dxv,dyu,arcd1,                                        &
                   yc,yu,yv,z,zo,Dvel,U,DU,V,DV,wwm,wwp,missing,vely)
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
   REALTYPE,dimension(E2DFIELD),intent(in)  :: dxv,dyu,arcd1
   REALTYPE,dimension(E2DFIELD),intent(in)  :: yc,yu,yv,z,zo,Dvel,U,DU,V,DV,wwm,wwp
   REALTYPE,intent(in)                      :: missing
!
! !OUTPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(out) :: vely
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
   write(debug,*) 'to_v() # ',Ncall
#endif
#ifdef SLICE_MODEL
!  Note (KK): this value MUST NOT be changed !!!
   j = jmax/2
#endif

!  KK-TODO: completely remove mask checks

   dtm1 = _ONE_/dt

!$OMP PARALLEL DEFAULT(SHARED)                                         &
!$OMP          FIRSTPRIVATE(j)                                         &
!$OMP          PRIVATE(i)


   select case(grid_type)

      case (1,2)

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
         do j=jmin-HALO+1,jmax+HALO
#endif
            do i=imin-HALO,imax+HALO
               if (az(i,j) .eq. 1) then
                  vely(i,j) = _HALF_*( V(i,j-1) + V(i,j) )/Dvel(i,j)
               !else
               !   vely(i,j) = missing
               end if
            end do
#ifndef SLICE_MODEL
         end do
#endif
!$OMP END DO
         !vely(:,jmin-HALO) = missing

      case (3)

!$OMP DO SCHEDULE(RUNTIME)
#ifndef SLICE_MODEL
         do j=jmin-HALO+1,jmax+HALO
#endif
            do i=imin-HALO+1,imax+HALO
               if (az(i,j) .eq. 1) then
                  vely(i,j) = (                                   &
                                 (                                &
                                     ( z(i,j) - zo(i,j) )*dtm1    &
                                   + ( wwp(i,j) - wwm(i,j) )      &
                                 )                                &
                                 *yc(i,j)                         &
                               + (                                &
                                    U(i  ,j  )*yu(i  ,j  )*DYU    &
                                  - U(i-1,j  )*yu(i-1,j  )*DYUIM1 &
                                  + V(i  ,j  )*yv(i  ,j  )*DXV    &
                                  - V(i  ,j-1)*yv(i  ,j-1)*DXVJM1 &
                                 )                                &
                                 *ARCD1                           &
                              )                                   &
                              /Dvel(i,j)
               !else
               !   vely(i,j) = missing
               end if
            end do
#ifndef SLICE_MODEL
         end do
#endif
!$OMP END DO
         !vely(imin-HALO,:) = missing
         !vely(:,jmin-HALO) = missing

      case (4)
         stop 'to_v: grid_type=4 not implemented yet'

      case default
         stop 'to_v: invalid grid_type'

   end select

!$OMP END PARALLEL

#ifdef SLICE_MODEL
   vely(:,j+1) = vely(:,j)
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving to_v()'
   write(debug,*)
#endif
   return
   end subroutine to_v
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2012 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
