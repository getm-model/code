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
   use domain, only: imin,imax,jmin,jmax,kmax,az
   use  parameters, only:   g,rho_0
   use  variables_3d, only: T,S,rho
   use salinity, only: nonnegsalt_method
!  IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   public init_eqstate, do_eqstate
!
! !PRIVATE DATA MEMBERS:
   integer                   :: eqstate_method=3
   REALTYPE                  :: T0 = 10., S0 = 33.75, p0 = 0.
   REALTYPE                  :: dtr0 = -0.17, dsr0 = 0.78
   REALTYPE,dimension(:,:,:),allocatable,target :: S_eos
   REALTYPE,dimension(:,:,:),pointer            :: pS
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
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
! !DESCRIPTION:
!  Reads the namelist and makes calls to the init functions of the
!  various model components.
!
! !LOCAL VARIABLES:
   integer :: rc
   namelist /eqstate/ eqstate_method,T0,S0,p0,dtr0,dsr0
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
      case (3)
         LEVEL3 'Jackett ea 2006 - without pressure in rho'
      case default
         FATAL 'init_eqstate(): not a valid eqstate_method'
         stop 'init_eqstate()'
   end select

   if (nonnegsalt_method .eq. 1) then
      allocate(S_eos(I3DFIELD),stat=rc)
      if (rc /= 0) stop 'init_eqstate: Error allocating memory (S_eos)'
      pS => S_eos
   else
      pS => S
   end if

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
   use parameters, only: rho_0
   use domain, only: imin,imax,jmin,jmax,kmax,az
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: arcd1
#else
   use domain, only: ard1
#endif
   use variables_3d, only: kmin,T,S,rho,buoy,hn,alpha,beta
   use getm_timers, only: tic, toc, TIM_EQSTATE
!$ use omp_lib
   IMPLICIT NONE
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   integer                   :: negpoints
   REALTYPE                  :: negvol,negsalt,negsalt_min
   REALTYPE                  :: x
   REALTYPE                  :: p1,s1,t1
   REALTYPE                  :: th,densp
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_eqstate() # ',Ncall
#endif
   call tic(TIM_EQSTATE)

!$OMP PARALLEL DEFAULT(SHARED)          &
!$OMP PRIVATE(i,j,k, x,p1,s1,t1,th,densp)

   if (nonnegsalt_method .gt. 0) then
!$OMP SINGLE
      negpoints = 0
      negvol = _ZERO_
      negsalt = _ZERO_
      negsalt_min = _ZERO_
!$OMP END SINGLE
! OMP-NOTE (KK): SINGLE implies BARRIER
      do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
         do j=jmin-HALO,jmax+HALO
            do i=imin-HALO,imax+HALO
               if (az(i,j) .ne. 0) then
                  s1 = S(i,j,k)
                  if (s1 .lt. _ZERO_) then
                     pS(i,j,k) = _ZERO_
                     if (      imin.le.i .and. i.le.imax &
                         .and. jmin.le.j .and. j.le.jmax ) then
!$OMP ATOMIC
                        negpoints = negpoints+1
!$OMP END ATOMIC
                        x = hn(i,j,k)/ARCD1
!$OMP ATOMIC
                        negvol = negvol + x
!$OMP END ATOMIC
!$OMP ATOMIC
                        negsalt = negsalt + s1*x
!$OMP END ATOMIC
!$OMP CRITICAL
                        negsalt_min = min(s1,negsalt_min)
#ifdef DEBUG
                        STDERR 'Salinity at point ',i,',',j,',',k,' < 0.'
                        STDERR 'Value is S = ',s1
                        STDERR 'Programm continued, value set to zero ...'
#endif
!$OMP END CRITICAL
                     end if
                  else if(nonnegsalt_method .eq. 1) then
                     pS(i,j,k) = s1
                  end if
               end if
            end do
         end do
!$OMP END DO NOWAIT
      end do
!$OMP BARRIER
!$OMP MASTER
      if (negpoints .gt. 0) then
         if (nonnegsalt_method .eq. 1) then
            STDERR "clipped ",negpoints," negative salinities (only inside EOS)"
         else
            STDERR "clipped ",negpoints," negative salinities"
         end if
         STDERR "(min:",real(negsalt_min)," psu", &
                ",mean:",real(negsalt/x)," psu",  &
                ",sum:",real(negsalt*rho_0/1000.0d0)," kg)"
      end if
!$OMP END MASTER
   end if

!  KK-TODO: do we still need BUOYANCY ?!
#define BUOYANCY
   select case (eqstate_method)
      case (1)
         do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
            do j=jmin-HALO,jmax+HALO
               do i=imin-HALO,imax+HALO
                  rho(i,j,k) = rho_0 +                                 &
                       dtr0*(T(i,j,k)-T0) + dsr0*(pS(i,j,k)-S0)
               end do
            end do
!$OMP END DO NOWAIT
         end do

#ifndef _OLD_BVF_
! First finished threas should just start intiializing arrays...
!$OMP SINGLE
         alpha = dtr0
!$OMP END SINGLE NOWAIT
!$OMP SINGLE
         beta  = dsr0
!$OMP END SINGLE NOWAIT
#endif
! This OMP barrier is necessary for the NOWAIT speedup of the rho-loop above
!$OMP BARRIER

      case (2)
         do k = 1,kmax
!$OMP DO SCHEDULE(RUNTIME)
            do j = jmin-HALO,jmax+HALO
               do i = imin-HALO,imax+HALO
                  if (az(i,j) .gt. 0) then
                     call rho_from_theta_unesco80(T(i,j,k),pS(i,j,k), &
                                                  rho(i,j,k))
                  end if
               end do
            end do
!$OMP END DO NOWAIT
         end do
! This OMP barrier is necessary for the NOWAIT speedup of the previous loop:
! OMP-TODO (KK): why not move this barrier behind #ifdef _OLD_BVF_ ?!
!$OMP BARRIER
#undef BUOYANCY

#ifndef _OLD_BVF_
!$OMP DO SCHEDULE(RUNTIME)
         do j = jmin-HALO,jmax+HALO
            do i = imin-HALO,imax+HALO
               if (az(i,j) .gt. 0) then

!                 not used anyway - but looks better
                  alpha(i,j,kmax) = _ZERO_
                  beta(i,j,kmax)  = _ZERO_

                  p1 = _ZERO_
                  do k = kmax-1,kmin(i,j),-1
                     th = _HALF_ * (T(i,j,k+1)+T(i,j,k))
                     s1 = _HALF_ * (pS(i,j,k+1)+pS(i,j,k))
                     p1 = p1+hn(i,j,k+1)
                     call eosall_from_theta(s1,th,p1,  &
                                            beta(i,j,k),alpha(i,j,k))
                  end do
               end if
            end do
         end do
!$OMP END DO
#endif
      case (3)
!        first calculate potential density
         do k = 1,kmax
!$OMP DO SCHEDULE(RUNTIME)
            do j = jmin-HALO,jmax+HALO
               do i = imin-HALO,imax+HALO
                  if (az(i,j) .gt. 0) then
                     call rho_from_theta(pS(i,j,k),T(i,j,k),_ZERO_,rho(i,j,k),densp)
                  end if
               end do
            end do
!$OMP END DO NOWAIT
         end do

#undef BUOYANCY

#ifndef _OLD_BVF_
!        calculate at SS interface, rho, alpha, beta
!$OMP DO SCHEDULE(RUNTIME)
         do j = jmin-HALO,jmax+HALO
            do i = imin-HALO,imax+HALO
               if (az(i,j) .gt. 0) then

!                 not used anyway - but looks better
                  alpha(i,j,kmax) = _ZERO_
                  beta(i,j,kmax)  = _ZERO_

                  p1 = _ZERO_
                  do k = kmax-1,kmin(i,j),-1
                     th = _HALF_ * (T(i,j,k+1)+T(i,j,k))
                     s1 = _HALF_ * (pS(i,j,k+1)+pS(i,j,k))
                     p1 = p1+hn(i,j,k+1)
                     call eosall_from_theta(s1,th,p1,  &
                                            beta(i,j,k),alpha(i,j,k))
                  end do
               end if
            end do
         end do
!$OMP END DO NOWAIT
#endif

! This OMP barrier is necessary for the NOWAIT speedup of the rho-loop above
!$OMP BARRIER

      case default
   end select

! GETM still uses surface potential density for buoy
! to be used in pressure gradient routines
   x= -g/rho_0
   do k=1,kmax
!$OMP DO SCHEDULE(RUNTIME)
      do j=jmin-HALO,jmax+HALO
         do i=imin-HALO,imax+HALO
            buoy(i,j,k)=x*(rho(i,j,k)-rho_0)
         end do
      end do
!$OMP END DO NOWAIT
   end do

!$OMP END PARALLEL

   call toc(TIM_EQSTATE)
#ifdef DEBUG
   write(debug,*) 'leaving do_eqstate()'
   write(debug,*)
#endif
   return
   end subroutine do_eqstate
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: rho_from_theta_unesco80

! !INTERFACE:
   subroutine rho_from_theta_unesco80(T,S,rho)
!
! !DESCRIPTION:
! Here, the equation of state is calculated using the UNESCO 1980 code.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)      :: T   ! potential temperature degC
   REALTYPE, intent(in)      :: S   ! salinity PSU
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)     :: rho ! density [kg/m3]
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
   REALTYPE                  :: x,T1,T2,T3,T4,T5,S1,S2,S3
!EOP
!-----------------------------------------------------------------------
!BOC
   T1 = T
   T2 = T1*T1
   T3 = T1*T2
   T4 = T2*T2
   T5 = T1*T4
   S1 = S
   S2 = S1*S1
   S3 = S1*S2

   x=999.842594+6.793952e-02*T1-9.09529e-03*T2+1.001685e-04*T3
   x=x-1.120083e-06*T4+6.536332e-09*T5
   x=x+S1*(0.824493-4.0899e-03*T1+7.6438e-05*T2-8.2467e-07*T3)
   x=x+S1*5.3875e-09*T4
   x=x+sqrt(S3)*(-5.72466e-03+1.0227e-04*T1-1.6546e-06*T2)
   x=x+4.8314e-04*S2
   rho = x
   return
   end subroutine rho_from_theta_unesco80
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: rho_from_theta
!
! !INTERFACE:
   subroutine rho_from_theta(s,th,p,dens0,densp)
!
! !DESCRIPTION:
! Here, the equation of state is calculated
! Uses Jackett ea 2006 algorithm specific for potential temperature
! Checkvalue
!   S=35
!   T=25
!   p=10000
!   rho_from_theta = 1062.53817
!
!   s                : salinity                           (psu)
!   th               : potential temperature              (deg C, ITS-90)
!   p                : gauge pressure                     (dbar)
!                      (absolute pressure - 10.1325 dbar)
!
!   rho_from_theta   : in-situ density                    (kg m^-3)
!
!   check value      : rho_from_theta(20,20,1000) = 1017.728868019642
!
!   based on DRJ on 10/12/03
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)      :: th  ! potential temperature degC
   REALTYPE, intent(in)      :: s  ! in situ salinity PSU
   REALTYPE, intent(in)      :: p  ! pressure in dbars
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)     :: dens0  ! density at 0.0 dbars
   REALTYPE, intent(out)     :: densp  ! density at p dbars
!
! !REVISION HISTORY:
!  AS 2009 based on code provided by Jackett 2005
!  See the log for the module
!
! !LOCAL VARIABLES
   REALTYPE th2,sqrts,anum,aden,pth
!EOP
!-----------------------------------------------------------------------
!BOC
   th2 = th*th; sqrts = sqrt(s)

   anum =          9.9984085444849347d+02 +    &
              th*( 7.3471625860981584d+00 +    &
              th*(-5.3211231792841769d-02 +    &
              th*  3.6492439109814549d-04)) +  &
               s*( 2.5880571023991390d+00 -    &
              th*  6.7168282786692355d-03 +    &
               s*  1.9203202055760151d-03)

   aden =          1.0000000000000000d+00 +    &
              th*( 7.2815210113327091d-03 +    &
              th*(-4.4787265461983921d-05 +    &
              th*( 3.3851002965802430d-07 +    &
              th*  1.3651202389758572d-10))) + &
               s*( 1.7632126669040377d-03 -    &
              th*( 8.8066583251206474d-06 +    &
             th2*  1.8832689434804897d-10) +   &
           sqrts*( 5.7463776745432097d-06 +    &
             th2*  1.4716275472242334d-09))

   dens0 = anum/aden
   densp = dens0

   if(p .ne. _ZERO_) then

      pth = p*th

      anum = anum +        p*( 1.1798263740430364d-02 +   &
                         th2*  9.8920219266399117d-08 +   &
                           s*  4.6996642771754730d-06 -   &
                           p*( 2.5862187075154352d-08 +   &
                         th2*  3.2921414007960662d-12))

      aden = aden +        p*( 6.7103246285651894d-06 -   &
                    pth*(th2*  2.4461698007024582d-17 +   &
                           p*  9.1534417604289062d-18))

      densp = anum/aden
   end if

   return
   end subroutine rho_from_theta
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: eosall_from_theta
!
! !INTERFACE:
   subroutine eosall_from_theta(s,th,p,rho_s,rho_th)
!
! !DESCRIPTION:
!   in-situ density and its derivatives (only 2) as functions of
!   salinity, potential temperature and pressure as in
!   Jackett, McDougall, Feistel, Wright and Griffies (2006), JAOT
!
!   s                : salinity                           (psu)
!   th               : potential temperature              (deg C, ITS-90)
!   p                : gauge pressure                     (dbar)
!                      (absolute pressure - 10.1325 dbar)
!
!   rho              : in-situ density                    (kg m^-3)
!   rho_s            : partial derivative wrt s           (kg m^-3 psu^-1)
!   rho_th           : partial derivative wrt th          (kg m^-3 deg C^-1)
!
!   check values     : eosall_from_theta(20,20,1000,...) gives
!
!                               rho =  1017.728868019642
!                               rho_s =   0.7510471164699279
!                               rho_th = -0.2570255211349140
!
!   based on DRJ on 10/12/03

! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)      :: th ! potential temperature degC
   REALTYPE, intent(in)      :: s  ! in situ salinity PSU
   REALTYPE, intent(in)      :: p  ! pressure in dbars
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)     :: rho_s   ! partial derivative wrt s
   REALTYPE, intent(out)     :: rho_th  ! partial derivative wrt th

! !REVISION HISTORY:
!  AS 2009 based on code provided by Jackett 2005
!  See the log for the module
!
! !LOCAL VARIABLES
   REALTYPE        :: th2,sqrts,anum,aden,pth
   REALTYPE        :: rho,anum_s,aden_s,anum_th,aden_th,rec_aden
!EOP
!-----------------------------------------------------------------------
!BOC
   th2 = th*th; sqrts = sqrt(s)

   anum =       9.9984085444849347d+02 +         &
           th*( 7.3471625860981584d+00 +         &
           th*(-5.3211231792841769d-02 +         &

           th*  3.6492439109814549d-04)) +       &
            s*( 2.5880571023991390d+00 -         &
           th*  6.7168282786692355d-03 +         &
            s*  1.9203202055760151d-03)

   aden =       1.0000000000000000d+00 +         &
           th*( 7.2815210113327091d-03 +         &
           th*(-4.4787265461983921d-05 +         &
           th*( 3.3851002965802430d-07 +         &
           th*  1.3651202389758572d-10))) +      &
            s*( 1.7632126669040377d-03 -         &
           th*( 8.8066583251206474d-06 +         &
          th2*  1.8832689434804897d-10) +        &
        sqrts*( 5.7463776745432097d-06 +         &
          th2*  1.4716275472242334d-09))


   anum_s =     2.5880571023991390d+00 -         &
           th*  6.7168282786692355d-03 +         &
            s*  3.8406404111520300d-03

   aden_s =     1.7632126669040377d-03 +        &
           th*(-8.8066583251206470d-06 -        &
          th2*  1.8832689434804897d-10) +       &
        sqrts*( 8.6195665118148150d-06 +        &
          th2*  2.2074413208363504d-09)

   anum_th =    7.3471625860981580d+00 +        &
           th*(-1.0642246358568354d-01 +        &
           th*  1.0947731732944364d-03)-        &
           s*  6.7168282786692355d-03

   aden_th =    7.2815210113327090d-03 +        &
           th*(-8.9574530923967840d-05 +        &
           th*( 1.0155300889740728d-06 +        &
           th*  5.4604809559034290d-10)) +      &
            s*(-8.8066583251206470d-06 -        &
          th2*  5.6498068304414700d-10 +        &
           th*sqrts*  2.9432550944484670d-09)

   if(p .ne. _ZERO_) then

     pth = p*th

     anum = anum +   p*( 1.1798263740430364d-02 +     &
                   th2*  9.8920219266399117d-08 +     &
                     s*  4.6996642771754730d-06 -     &
                     p*( 2.5862187075154352d-08 +     &
                   th2*  3.2921414007960662d-12))

     aden = aden + p*( 6.7103246285651894d-06 -       &
                 pth*(th2*  2.4461698007024582d-17 +  &
                   p*  9.1534417604289062d-18))


     anum_s = anum_s +  p*  4.6996642771754730d-06

     anum_th = anum_th + pth*( 1.9784043853279823d-07 - &
                           p*  6.5842828015921320d-12)

     aden_th = aden_th -                                &
                    p*p*(th2*  7.3385094021073750d-17 + &
                      p*  9.1534417604289060d-18)

   end if

   rec_aden = _ONE_ / aden

   rho = anum*rec_aden

   rho_s = (anum_s-aden_s*rho)*rec_aden

   rho_th = (anum_th-aden_th*rho)*rec_aden

!  saline contraction coefficient is rho_s/rho
!  thermal expansion coefficient is -rho_th/rho

   return
   end subroutine eosall_from_theta
!EOC

!-----------------------------------------------------------------------

   end module eqstate

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------

