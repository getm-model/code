!$Id: eqstate.F90,v 1.9 2007-02-23 12:20:36 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  eqstate \label{sec-eqstate}
!
! !INTERFACE:
   module eqstate
!
! !DESCRIPTION:
!
! Documentation will follow when the equation of state calculations
! are updated. The idea is to use the respective routines from GOTM.
!
! !USES:
   use  parameters, only: g,rho_0
!  IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   public init_eqstate, do_eqstate
!
! !PRIVATE DATA MEMBERS:
   integer                   :: eqstate_method=1
   REALTYPE                  :: T0 = 10., S0 = 33.75, p0 = 0.
   REALTYPE                  :: dtr0 = -0.17, dsr0 = 0.78
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_eqstate
!
! !INTERFACE:
   subroutine init_eqstate()
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Reads the namelist and makes calls to the init functions of the
!  various model components.
!
! !LOCAL VARIABLES:
   namelist /eqstate/ eqstate_method,T0,S0,p0,dtr0,dsr0
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_eqstate() # ',Ncall
#endif

   LEVEL2 'init_eqstate()'
   read(NAMLST,nml=eqstate)

   select case (eqstate_method)
      case (1)
         LEVEL3 'linear equation of state'
         LEVEL4 'T0   = ',T0
         LEVEL4 'S0   = ',S0
         LEVEL4 'dtr0 = ',dtr0
         LEVEL4 'dsr0 = ',dsr0
      case (2)
         LEVEL3 'UNESCO - no pressure adjustment'
      case default
         FATAL 'init_eqstate(): not a valid eqstate_method'
         stop 'init_eqstate()'
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving init_eqstate()'
   write(debug,*)
#endif
   return
   end subroutine init_eqstate
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_eqstate - equation of state \label{sec-do-eqstate} 
!
! !INTERFACE:
   subroutine do_eqstate()
!
! !DESCRIPTION:
!
! Here, the equation of state is calculated for every 3D grid point.  
!
! !USES:
   use domain, only: iimin,iimax,jjmin,jjmax,kmax,az
   use variables_3d, only: T,S,rho,buoy
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   REALTYPE                  :: x
   REALTYPE                  :: KK
   REALTYPE                  :: T1,T2,T3,T4,T5,S1,S15,S2,S3,p2
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_eqstate() # ',Ncall
#endif

#define BUOYANCY
   select case (eqstate_method)
      case (1)
         forall(i=iimin-1:iimax+1,j=jjmin-1:jjmax+1,az(i,j) .gt. 0)  &
            rho(i,j,1:kmax) = rho_0 +                                &
                    dtr0*(T(i,j,1:kmax)-T0) + dsr0*(S(i,j,1:kmax)-S0)

#ifdef HAIDVOGEL_TEST
         forall(i=iimin-1:iimax+1,j=jjmin-1:jjmax+1,az(i,j) .gt. 0)  &
            rho(i,j,1:kmax) = 1000. + S(i,j,1:kmax)
#endif
#ifdef CONSTANCE_TEST
         forall(i=iimin-1:iimax+1,j=jjmin-1:jjmax+1,az(i,j) .gt. 0)  &
          rho(i,j,1:kmax) = 1000. + *dtr0*(T(i,j,1:kmax)-T0)
#endif
      case (2)
         do k = 1,kmax
            do j = jjmin-1,jjmax+1
               do i = iimin-1,iimax+1
                  if (az(i,j) .gt. 0) then
                     T1 = T(i,j,k)
                     T2 = T1*T1
                     T3 = T1*T2
                     T4 = T2*T2
                     T5 = T1*T4
                     S1 = S(i,j,k)
#ifdef DEBUG
                     if (S1 .lt. _ZERO_) then
                        STDERR 'Salinity at point ',i,',',j,',',k,' < 0.'
                        STDERR 'Value is S = ',S(i,j,k)
                        STDERR 'Programm continued, value set to zero ...'
                        S(i,j,k)= _ZERO_
                     end if
#endif
                     S15= S1**1.5
                     S2 = S1*S1
                     S3 = S1*S2

                     x=999.842594+6.793952e-02*T1-9.09529e-03*T2+1.001685e-04*T3
                     x=x-1.120083e-06*T4+6.536332e-09*T5
                     x=x+S1*(0.824493-4.0899e-03*T1+7.6438e-05*T2-8.2467e-07*T3)
                     x=x+S1*5.3875e-09*T4
                     x=x+sqrt(S3)*(-5.72466e-03+1.0227e-04*T1-1.6546e-06*T2)
                     x=x+4.8314e-04*S2

#ifdef UNPRESS
                     if ((p.gt.0)) then
                        p2=p*p
                        KK= 19652.21                                                &
                          +148.4206     *T1       -2.327105    *T2        &
                          +  1.360477E-2*T3       -5.155288E-5 *T4        &
                          +  3.239908      *p     +1.43713E-3  *T *p      &
                          +  1.16092E-4 *T2*p     -5.77905E-7  *T3*p      &
                          +  8.50935E-5    *p2    -6.12293E-6  *T *p2     &
                          +  5.2787E-8  *T2*p2                            &
                          + 54.6746           *S1 -0.603459    *T    *S1  &
                          +  1.09987E-2 *T2   *S1 -6.1670E-5   *T3   *S1  &
                          +  7.944E-2         *S15+1.6483E-2   *T    *S15 &
                          -  5.3009E-4  *T2   *S15+2.2838E-3      *p *S1  &
                          -  1.0981E-5  *T1*p *S1 -1.6078E-6   *T2*p *S1  &
                          +  1.91075E-4    *p *S15-9.9348E-7      *p2*S1  &
                          +  2.0816E-8  *T1*p2*S1   +9.1697E-10  *T2*p2*S1
                        x=x/(1.-p/KK)
                     end if
#endif
                     rho(i,j,k)=x
                  end if
               end do
            end do
         end do

      case default
   end select
#undef BUOYANCY

   x=-g/rho_0
   buoy=x*(rho-rho_0)

#ifdef DEBUG
   write(debug,*) 'leaving do_eqstate()'
   write(debug,*)
#endif
   return
   end subroutine do_eqstate
!EOC

!-----------------------------------------------------------------------

   end module eqstate

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------

!#define TEST_EQSTATE
#ifdef TEST_EQSTATE
   program test_eqstate

   use eqstate

   REALTYPE T(1,1,1),S(1,1,1),rho(1,1,1)

   eqstate_method=2

   call do_eqstate()

   end
#endif
