#include "cppdefs.h"
!----------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_mean_fields() - produces averaged output.
!
! !INTERFACE:
   subroutine calc_mean_fields(n,write_mean,runtype)
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imax,imin,jmax,jmin,kmax
   use domain, only: az,au,av
   use meteo, only: swr
   use m3d, only: M,update_temp,update_salt
   use m3d, only: nonhyd_method
   use variables_3d, only: do_numerical_analyses_3d
   use variables_3d, only: hn,uu,hun,vv,hvn,ww,taub
#ifndef NO_BAROCLINIC
   use variables_3d, only: S,T,rho
   use variables_3d, only: buoy
#endif
   use variables_3d, only: minus_bnh
   use variables_3d, only: nummix_S,nummix_T
   use variables_3d, only: nummix_S_old,nummix_S_int,nummix_T_old,nummix_T_int
   use variables_3d, only: phymix_S,phymix_S_int,phymix_T,phymix_T_int
   use variables_3d, only: numdis_3d,numdis_3d_old,numdis_int
   use variables_3d, only: phydis_3d,phydis_int
#ifdef GETM_BIO
   use bio, only: bio_calc
   use bio_var, only: numc
   use variables_3d,  only: cc3d
#endif
#ifdef _FABM_
   use getm_fabm, only: fabm_pel,fabm_ben,fabm_diag,fabm_diag_hz
#endif
   use output, only: save_s,save_t,save_rho
   use diagnostic_variables
   use getm_timers, only: tic, toc, TIM_CALCMEANF
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: n,runtype
   logical, intent(in)  :: write_mean
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Adolf Stips
!
! !LOCAL VARIABLES:
   integer         :: i,j,k,rc
   REALTYPE        :: tmpf(I3DFIELD)
   integer,save    :: step=0
   logical,save    :: first=.true.
   logical,save    :: fabm_mean=.false.
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'calc_mean_fields() # ',Ncall
#endif
   call tic(TIM_CALCMEANF)

   if (first ) then
      LEVEL3 'calc_mean_fields(): initialising variables'
      allocate(swrmean(E2DFIELD),stat=rc)
      if (rc /= 0) &
          stop 'calc_mean_fields.F90: Error allocating memory (swrmean)'
      allocate(ustarmean(E2DFIELD),stat=rc)
      if (rc /= 0) &
          stop 'calc_mean_fields.F90: Error allocating memory (ustarmean)'
      allocate(ustar2mean(E2DFIELD),stat=rc)
      if (rc /= 0) &
          stop 'calc_mean_fields.F90: Error allocating memory (ustar2mean)'
      allocate(uumean(I3DFIELD),stat=rc)
      if (rc /= 0) &
          stop 'calc_mean_fields.F90: Error allocating memory (uumean)'
      allocate(vvmean(I3DFIELD),stat=rc)
      if (rc /= 0) &
          stop 'calc_mean_fields.F90: Error allocating memory (vvmean)'
      allocate(wmean(I3DFIELD),stat=rc)
      if (rc /= 0) &
          stop 'calc_mean_fields.F90: Error allocating memory (wmean)'
      allocate(humean(I3DFIELD),stat=rc)
      if (rc /= 0) &
          stop 'calc_mean_fields.F90: Error allocating memory (humean)'
      allocate(hvmean(I3DFIELD),stat=rc)
      if (rc /= 0) &
          stop 'calc_mean_fields.F90: Error allocating memory (hvmean)'
      allocate(hmean(I3DFIELD),stat=rc)
      if (rc /= 0) &
          stop 'calc_mean_fields.F90: Error allocating memory (hmean)'
#ifndef NO_BAROCLINIC
      if (save_t) then
         allocate(Tmean(I3DFIELD),stat=rc)
         if (rc /= 0) &
             stop 'calc_mean_fields.F90: Error allocating memory (Tmean)'
      end if
      if (save_s) then
         allocate(Smean(I3DFIELD),stat=rc)
         if (rc /= 0) &
             stop 'calc_mean_fields.F90: Error allocating memory (Smean)'
      end if
      if (save_rho) then
         allocate(rhomean(I3DFIELD),stat=rc)
         if (rc /= 0) &
             stop 'calc_mean_fields.F90: Error allocating memory (rhomean)'
      end if
#endif
      if (nonhyd_method .ne. 0) then
         allocate(bnhmean(I3DFIELD),stat=rc)
           if (rc /= 0) &
              stop 'calc_mean_fields.F90: Error allocating memory (bnhmean)'
      end if
      if (do_numerical_analyses_3d) then
         allocate(numdis_3d_mean(I3DFIELD),stat=rc)
           if (rc /= 0) &
              stop 'calc_mean_fields.F90: Error allocating memory (numdis_3d_mean)'
#ifdef _NUMERICAL_ANALYSES_OLD_
         allocate(numdis_3d_old_mean(I3DFIELD),stat=rc)
           if (rc /= 0) &
              stop 'calc_mean_fields.F90: Error allocating memory (numdis_3d_old_mean)'
         allocate(numdis_int_mean(I2DFIELD),stat=rc)
           if (rc /= 0) &
              stop 'calc_mean_fields.F90: Error allocating memory (numdis_int_mean)'
#endif
         allocate(phydis_3d_mean(I3DFIELD),stat=rc)
           if (rc /= 0) &
              stop 'calc_mean_fields.F90: Error allocating memory (phydis_3d_mean)'
         allocate(phydis_int_mean(I2DFIELD),stat=rc)
           if (rc /= 0) &
              stop 'calc_mean_fields.F90: Error allocating memory (phydis_int_mean)'
         if (update_temp) then
            allocate(nummix_T_mean(I3DFIELD),stat=rc)
            if (rc /= 0) &
               stop 'calc_mean_fields.F90: Error allocating memory (nummix_T_mean)'
#ifdef _NUMERICAL_ANALYSES_OLD_
            allocate(nummix_T_old_mean(I3DFIELD),stat=rc)
            if (rc /= 0) &
               stop 'calc_mean_fields.F90: Error allocating memory (nummix_T_old_mean)'
            allocate(nummix_T_int_mean(I2DFIELD),stat=rc)
            if (rc /= 0) &
               stop 'calc_mean_fields.F90: Error allocating memory (nummix_T_int_mean)'
#endif
            allocate(phymix_T_mean(I3DFIELD),stat=rc)
            if (rc /= 0) &
               stop 'calc_mean_fields.F90: Error allocating memory (nummix_T_mean)'
            allocate(phymix_T_int_mean(I2DFIELD),stat=rc)
            if (rc /= 0) &
               stop 'calc_mean_fields.F90: Error allocating memory (nummix_T_int_mean)'
         end if
         if (update_salt) then
            allocate(nummix_S_mean(I3DFIELD),stat=rc)
            if (rc /= 0) &
               stop 'calc_mean_fields.F90: Error allocating memory (nummix_S_mean)'
#ifdef _NUMERICAL_ANALYSES_OLD_
            allocate(nummix_S_old_mean(I3DFIELD),stat=rc)
            if (rc /= 0) &
               stop 'calc_mean_fields.F90: Error allocating memory (nummix_S_old_mean)'
            allocate(nummix_S_int_mean(I2DFIELD),stat=rc)
            if (rc /= 0) &
               stop 'calc_mean_fields.F90: Error allocating memory (nummix_S_int_mean)'
#endif
            allocate(phymix_S_mean(I3DFIELD),stat=rc)
            if (rc /= 0) &
               stop 'calc_mean_fields.F90: Error allocating memory (nummix_S_mean)'
            allocate(phymix_S_int_mean(I2DFIELD),stat=rc)
            if (rc /= 0) &
               stop 'calc_mean_fields.F90: Error allocating memory (nummix_S_int_mean)'
         end if
      end if
#ifdef GETM_BIO
      allocate(cc3dmean(numc,I3DFIELD),stat=rc)
      if (rc /= 0) &
          stop 'calc_mean_fields.F90: Error allocating memory (cc3dmean)'
#endif
#ifdef _FABM_
      if (allocated(fabm_pel)) then
         fabm_mean=.true.
         allocate(fabmmean_pel(I3DFIELD,ubound(fabm_pel,4)),stat=rc)
         if (rc /= 0) &
            stop 'calc_mean_fields.F90: Error allocating memory (fabmmean_pel)'
         allocate(fabmmean_ben(I2DFIELD,ubound(fabm_ben,3)),stat=rc)
         if (rc /= 0) &
            stop 'calc_mean_fields.F90: Error allocating memory (fabmmean_ben)'
         allocate(fabmmean_diag(I3DFIELD,ubound(fabm_diag,4)),stat=rc)
         if (rc /= 0) &
            stop 'calc_mean_fields.F90: Error allocating memory (fabmmean_diag)'
         allocate(fabmmean_diag_hz(I2DFIELD,ubound(fabm_diag_hz,3)),stat=rc)
         if (rc /= 0) &
            stop 'calc_mean_fields.F90: Error allocating memory (fabmmean_diag_hz)'
      else
         fabm_mean=.false.
      end if
#endif
      first = .false.
   end if


!  Sum every macro time step, even less would be okay
   if(mod(n,M) .eq. 0) then

!     reset to start new meanout period
      if (step .eq. 0) then
         uumean=_ZERO_; vvmean=_ZERO_; wmean=_ZERO_
         humean=_ZERO_; hvmean=_ZERO_; hmean=_ZERO_
#ifndef NO_BAROCLINIC
         if (save_t) Tmean=_ZERO_
         if (save_s) Smean=_ZERO_
         if (save_rho) rhomean=_ZERO_
#endif
         if (nonhyd_method .ne. 0) bnhmean=_ZERO_
         if (do_numerical_analyses_3d) then
            numdis_3d_mean=_ZERO_
#ifdef _NUMERICAL_ANALYSES_OLD_
            numdis_3d_old_mean=_ZERO_; numdis_int_mean=_ZERO_
#endif
            phydis_3d_mean=_ZERO_; phydis_int_mean=_ZERO_
            if (update_temp) then
               nummix_T_mean=_ZERO_
#ifdef _NUMERICAL_ANALYSES_OLD_
               nummix_T_old_mean=_ZERO_; nummix_T_int_mean=_ZERO_
#endif
               phymix_T_mean=_ZERO_; phymix_T_int_mean=_ZERO_
            end if
            if (update_salt) then
               nummix_S_mean=_ZERO_
#ifdef _NUMERICAL_ANALYSES_OLD_
               nummix_S_old_mean=_ZERO_; nummix_S_int_mean=_ZERO_
#endif
               phymix_S_mean=_ZERO_; phymix_S_int_mean=_ZERO_
            end if
         end if
#ifdef GETM_BIO
         cc3dmean=_ZERO_
#endif
#ifdef _FABM_
         if (fabm_mean) then
            fabmmean_pel=_ZERO_
            fabmmean_ben=_ZERO_
            fabmmean_diag=_ZERO_
            fabmmean_diag_hz=_ZERO_
         end if
#endif
         ustarmean=_ZERO_; ustar2mean=_ZERO_; swrmean=_ZERO_
      end if

      swrmean = swrmean + swr
!     AS this has to be checked, if it is the correct ustar,
!     so we must not divide by rho_0 !!
      ustarmean = ustarmean + sqrt(taub)
      ustar2mean = ustar2mean + (taub)

      uumean = uumean + uu
      vvmean = vvmean + vv

!  calculate the real vertical velocities
!KBK - the towas done by Adolf Stips has some errors. For now the mean
!vertical velocity is the grid-ralated velocity.
#if 0
      tmpf=_ZERO_
      call towas(tmpf)
      wmean = wmean + tmpf
#else
      wmean = wmean + ww
#endif

      humean = humean + hun
      hvmean = hvmean + hvn
      hmean = hmean + hn

#ifndef NO_BAROCLINIC
      if (save_t) Tmean = Tmean + T*hn
      if (save_s) Smean = Smean + S*hn
      if (save_rho) rhomean = rhomean + rho*hn
#endif
      if (nonhyd_method .ne. 0) then
         if (runtype.eq.2 .or. nonhyd_method.eq.1) then
            bnhmean = bnhmean - minus_bnh*hn
         else
#ifndef NO_BAROCLINIC
            bnhmean = bnhmean + abs(minus_bnh/max(buoy,SMALL))
#endif
         end if
      end if

      if (do_numerical_analyses_3d) then
         numdis_3d_mean = numdis_3d_mean + numdis_3d*hn
#ifdef _NUMERICAL_ANALYSES_OLD_
         numdis_3d_old_mean = numdis_3d_old_mean + numdis_3d_old*hn
         numdis_int_mean = numdis_int_mean + numdis_int
#endif
         phydis_3d_mean = phydis_3d_mean + phydis_3d*hn
         phydis_int_mean = phydis_int_mean + phydis_int
         if (update_temp) then
            nummix_T_mean = nummix_T_mean + nummix_T*hn
#ifdef _NUMERICAL_ANALYSES_OLD_
            nummix_T_old_mean = nummix_T_old_mean + nummix_T_old*hn
            nummix_T_int_mean = nummix_T_int_mean + nummix_T_int
#endif
            phymix_T_mean = phymix_T_mean + phymix_T*hn
            phymix_T_int_mean = phymix_T_int_mean + phymix_T_int
         end if
         if (update_salt) then
            nummix_S_mean = nummix_S_mean + nummix_S*hn
#ifdef _NUMERICAL_ANALYSES_OLD_
            nummix_S_old_mean = nummix_S_old_mean + nummix_S_old*hn
            nummix_S_int_mean = nummix_S_int_mean + nummix_S_int
#endif
            phymix_S_mean = phymix_S_mean + phymix_S*hn
            phymix_S_int_mean = phymix_S_int_mean + phymix_S_int
         end if
      end if
#ifdef GETM_BIO
      if (bio_calc) cc3dmean=cc3dmean + cc3d
#endif
#ifdef _FABM_
      if (fabm_mean) then
          fabmmean_pel = fabmmean_pel + fabm_pel
          fabmmean_ben = fabmmean_ben + fabm_ben
          fabmmean_diag = fabmmean_diag + fabm_diag
          fabmmean_diag_hz = fabmmean_diag_hz + fabm_diag_hz
      end if
#endif
!  count them
      step = step + 1
   end if   ! here we summed them up

!  prepare for output
   if (write_mean) then

      if (step .gt. 1) then
         uumean = uumean / step
         vvmean = vvmean / step
         wmean = wmean / step
         humean = humean / step
         hvmean = hvmean / step
         hmean = hmean / step

#ifndef NO_BAROCLINIC
         if (save_t) Tmean = Tmean / step
         if (save_s) Smean = Smean / step
         if (save_rho) rhomean = rhomean / step
#endif

         if (nonhyd_method .ne. 0) then
           if (runtype.eq.2 .or. nonhyd_method.eq.1) then
             bnhmean = bnhmean / step / hmean
           else
             bnhmean = bnhmean / step
           end if
         end if

         if (do_numerical_analyses_3d) then
            numdis_3d_mean = numdis_3d_mean / step / hmean
#ifdef _NUMERICAL_ANALYSES_OLD_
            numdis_3d_old_mean = numdis_3d_old_mean / step / hmean
            numdis_int_mean = numdis_int_mean / step
#endif
            phydis_3d_mean = phydis_3d_mean / step / hmean
            phydis_int_mean = phydis_int_mean / step
            if (update_temp) then
               nummix_T_mean = nummix_T_mean / step / hmean
#ifdef _NUMERICAL_ANALYSES_OLD_
               nummix_T_old_mean = nummix_T_old_mean / step / hmean
               nummix_T_int_mean = nummix_T_int_mean / step
#endif
               phymix_T_mean = phymix_T_mean / step / hmean
               phymix_T_int_mean = phymix_T_int_mean / step
            end if
            if (update_salt) then
               nummix_S_mean = nummix_S_mean / step / hmean
#ifdef _NUMERICAL_ANALYSES_OLD_
               nummix_S_old_mean = nummix_S_old_mean / step / hmean
               nummix_S_int_mean = nummix_S_int_mean / step
#endif
               phymix_S_mean = phymix_S_mean / step / hmean
               phymix_S_int_mean = phymix_S_int_mean / step
            end if
         end if
#ifdef GETM_BIO
         if (bio_calc) cc3dmean = cc3dmean / step
#endif
#ifdef _FABM_
         if (fabm_mean) then
            fabmmean_pel = fabmmean_pel / step
            fabmmean_ben = fabmmean_ben / step
            fabmmean_diag = fabmmean_diag / step
            fabmmean_diag_hz = fabmmean_diag_hz / step
         end if
#endif
         ustarmean = ustarmean / step
         swrmean = swrmean / step

      end if

      if (step .ge. 1) then

#ifndef NO_BAROCLINIC
         if (save_t) then
            forall (i=imin:imax,j=jmin:jmax,az(i,j).ne.0)
               Tmean(i,j,1:) = Tmean(i,j,1:) / hmean(i,j,1:)
            end forall
         end if
         if (save_s) then
            forall (i=imin:imax,j=jmin:jmax,az(i,j).ne.0)
               Smean(i,j,1:) = Smean(i,j,1:) / hmean(i,j,1:)
            end forall
         end if
         if (save_rho) then
            forall (i=imin:imax,j=jmin:jmax,az(i,j).ne.0)
               rhomean(i,j,1:) = rhomean(i,j,1:) / hmean(i,j,1:)
            end forall
         end if
#endif


!  now calculate the velocities
         where ( humean .ne. _ZERO_ )
            uumean = uumean/humean
         elsewhere
            uumean =  _ZERO_
         end where

         where ( hvmean .ne. _ZERO_ )
            vvmean = vvmean/hvmean
         elsewhere
            vvmean = _ZERO_
         end where

!  we must destagger,  yes

         tmpf = _ZERO_
         do j=jmin,jmax
            do i=imin,imax
!  check if we are in the water
               if(au(i,j) .gt. 0 .and. au(i-1,j) .gt. 0) then
                  do k = 1, kmax
                     tmpf(i,j,k)=(uumean(i,j,k)+uumean(i-1,j,k))/2.0
                  end do !k
               end if
            end do
         end do
         uumean = tmpf

         tmpf = _ZERO_
         do j=jmin,jmax
            do i=imin,imax
!  check if we are in the water
               if(av(i,j) .gt. 0 .and. av(i,j-1) .gt. 0) then
                  do k = 1, kmax
                     tmpf(i,j,k)=(vvmean(i,j,k)+vvmean(i,j-1,k))/2.0
                  end do !k
               end if
            end do
         end do
         vvmean = tmpf

         tmpf = 0.0
         do j=jmin,jmax
            do i=imin,imax
!  check if we are in the water
               if(az(i,j) .gt. 0) then
                  tmpf(i,j,1)=wmean(i,j,1)/2.0
                  do k = 2, kmax
                     tmpf(i,j,k) = (wmean(i,j,k)+wmean(i,j,k-1))/2.0
                  end do
               end if
            end do
         end do
         wmean = tmpf
      end if
      step = 0
   end if

   call toc(TIM_CALCMEANF)
   return
   end subroutine calc_mean_fields
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2004 -  Adolf Stips  & Karsten Bolding                 !
!-----------------------------------------------------------------------
