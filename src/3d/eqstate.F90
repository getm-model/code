!$Id: eqstate.F90,v 1.4 2004-01-02 13:54:24 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  eqstate
!
! !INTERFACE:
   module eqstate
!
! !DESCRIPTION:
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
!  $Log: eqstate.F90,v $
!  Revision 1.4  2004-01-02 13:54:24  kbk
!  read equation of state info from namelist - Ruiz
!
!  Revision 1.3  2003/04/23 12:16:34  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 13:19:48  kbk
!  parallel support
!
!  Revision 1.1.1.1  2002/05/02 14:00:59  gotm
!  recovering after CVS crash
!
!  Revision 1.9  2001/10/17 13:56:02  bbh
!  Density for Constance case
!
!  Revision 1.8  2001/09/03 13:04:29  bbh
!  Small change for HAIDVOGEL
!
!  Revision 1.7  2001/08/27 11:50:17  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.6  2001/07/26 13:23:50  bbh
!  Partial new eq. of state
!
!  Revision 1.5  2001/07/26 12:54:11  bbh
!  Testing advection schems - using ifdef HAIDVOGEL_TEST
!
!  Revision 1.4  2001/05/23 11:23:51  bbh
!  Removed BAROCLINIC
!
!  Revision 1.3  2001/05/11 18:41:23  bbh
!  define BAROCLINIC for testing
!
!  Revision 1.2  2001/05/11 13:46:00  bbh
!  Added linear equation of state
!
!  Revision 1.1  2001/05/10 11:30:16  bbh
!  Added further support for baroclinicity
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
! !REVISION HISTORY:
!  See the log for the module
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
! !IROUTINE:  do_eqstate()
!
! !INTERFACE:
   subroutine do_eqstate()
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: iimin,iimax,jjmin,jjmax,kmax,az
   use variables_3d, only: T,S,rho
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   REALTYPE                  :: x
#if 0
!  Brydon et. al. - Table 2 - narrow
   REALTYPE, parameter       :: a1=-1.36471e-1,b1= 5.06423e-1,g1=-5.52640e-4
   REALTYPE, parameter       :: a2=-4.68181e-2,b2=-3.57109e-3,g2= 4.88584e-6
   REALTYPE, parameter       :: a3= 8.07004e-1,b3=-8.76148e-4,g3= 9.96027e-7
   REALTYPE, parameter       :: a4=-7.45353e-3,b4= 5.25243e-5,g4=-7.25139e-8
   REALTYPE, parameter       :: a5=-2.94418e-3,b5= 1.57976e-5,g5=-3.98736e-9
   REALTYPE, parameter       :: a6= 3.43570e-5,b6=-3.46686e-7,g6= 4.00631e-10
   REALTYPE, parameter       :: a7= 3.48658e-5,b7=-1.68764e-7,g7= 8.26368e-11
!
   REALTYPE                  :: p=50.,p2
   REALTYPE                  :: c1,c2,c3,c4,c5,c6,c7
!
   REALTYPE                  :: s1,t1,t2,t3
#endif
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
#ifdef DENSITY
         forall(i=iimin-1:iimax+1,j=jjmin-1:jjmax+1,az(i,j) .gt. 0)  &
            rho(i,j,1:kmax) = rho_0 +                                &
                    dtr0*(T(i,j,1:kmax)-T0) + dsr0*(S(i,j,1:kmax)-S0)
#endif
#undef DENSITY
#ifdef BUOYANCY
         x = -g/rho_0
         forall(i=iimin-1:iimax+1,j=jjmin-1:jjmax+1,az(i,j) .gt. 0)  &
            rho(i,j,1:kmax) =                                        &
                   x*(dtr0*(T(i,j,1:kmax)-T0) + dsr0*(S(i,j,1:kmax)-S0))
#endif

#ifdef HAIDVOGEL_TEST
         forall(i=iimin-1:iimax+1,j=jjmin-1:jjmax+1,az(i,j) .gt. 0)  &
            rho(i,j,1:kmax) = -g/rho_0*(S(i,j,1:kmax)+1000.-rho_0)
#endif
#ifdef CONSTANCE_TEST
         x = -g/rho_0
         forall(i=iimin-1:iimax+1,j=jjmin-1:jjmax+1,az(i,j) .gt. 0)  &
          rho(i,j,1:kmax) = x*dtr0*(T(i,j,1:kmax)-T0)
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
                     if (S1 .lt. _ZERO_) then
                        STDERR 'Salinity at point ',i,',',j,',',k,' < 0.'
                        STDERR 'Value is S = ',S(i,j,k)
                        STDERR 'Programm continued, value set to zero ...'
                        S(i,j,k)= _ZERO_
                     end if
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
                     rho(i,j,k)=-g*(x-rho_0)/rho_0
                  end if
               end do
            end do
         end do

#if 0
         do k = 1,kmax
            do j = jjmin-1:jjmax+1
               do i = iimin-1:iimax+1
                  if (az(i,j) .gt. 0) then
                     p = 50.
                     p2 = p*p
                     c1 = a1 + b1*p + g1*p2
                     c2 = a2 + b2*p + g2*p2
                     c3 = a3 + b3*p + g3*p2
                     c4 = a4 + b4*p + g4*p2
                     c5 = a5 + b5*p + g5*p2
                     c6 = a6 + b6*p + g6*p2
                     c7 = a7 + b7*p + g7*p2
                     s1 = S(i,j,k)
                     t1 = T(i,j,k)
                     t2 = t1*t1
                     t3 = t1*t2
                     rho(i,j,k)=c1+c2*t1+c3*s1+c4*t2+c5*s1*t1+c6*t3+c7*s1*t2
                     rho(i,j,k)=-g/rho_0*(rho(i,j,k)-rho_0)
                  end if
               end do
            end do
         end do
#endif
      case default
   end select
#undef BUOYANCY

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
