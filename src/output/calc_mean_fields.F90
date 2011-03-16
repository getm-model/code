!$Id: calc_mean_fields.F90,v 1.5 2009-08-18 10:24:46 bjb Exp $
#include "cppdefs.h"
!----------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_mean_fields() - produces averaged output.
!
! !INTERFACE:
   subroutine calc_mean_fields(n,meanout)
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: imax,imin,jmax,jmin,kmax
   use domain, only: az,au,av
   use meteo, only: swr
   use m3d, only: M,calc_temp,calc_salt
   use variables_3d, only: do_mixing_analysis
   use variables_3d, only: hn,uu,hun,vv,hvn,ww,taub
#ifndef NO_BAROCLINIC
   use variables_3d, only: S,T
   use variables_3d, only: nummix3d_S,nummix2d_S,nummix3d_T,nummix2d_T
   use variables_3d, only: phymix3d_S,phymix2d_S,phymix3d_T,phymix2d_T
#endif
#ifdef GETM_BIO
   use bio, only: bio_calc
   use bio_var, only: numc
   use variables_3d,  only: cc3d
#endif
   use diagnostic_variables
   use getm_timers, only: tic, toc, TIM_CALCMEANF
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: n,meanout
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Adolf Stips
!
!  $Log: calc_mean_fields.F90,v $
!  Revision 1.5  2009-08-18 10:24:46  bjb
!  New getm_timers module
!
!  Revision 1.4  2007-06-07 10:25:20  kbk
!  iimin,iimax,jjmin,jjmax -> imin,imax,jmin,jmax
!
!  Revision 1.3  2006-01-27 16:07:57  kbk
!  now works with -DNO_BAROCLINIC
!
!  Revision 1.2  2004/04/21 09:11:44  lars
!  removed tab
!
!  Revision 1.1  2004/03/29 15:35:52  kbk
!  possible to store calculated mean fields
!
!
! !LOCAL VARIABLES:
   integer         :: i,j,k,rc
   REALTYPE        :: tmpf(I3DFIELD)
   REALTYPE,save   :: step=_ZERO_
   logical,save    :: first=.true.
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
      allocate(Tmean(I3DFIELD),stat=rc)
      if (rc /= 0) &
          stop 'calc_mean_fields.F90: Error allocating memory (Tmean)'
      allocate(Smean(I3DFIELD),stat=rc)
      if (rc /= 0) &
          stop 'calc_mean_fields.F90: Error allocating memory (Smean)'

      if (do_mixing_analysis) then
         if (calc_temp) then
            allocate(nummix3d_T_mean(I3DFIELD),stat=rc)
            if (rc /= 0) &
               stop 'calc_mean_fields.F90: Error allocating memory (nummix3d_T_mean)'
            allocate(nummix2d_T_mean(I2DFIELD),stat=rc)
            if (rc /= 0) &
               stop 'calc_mean_fields.F90: Error allocating memory (nummix2d_T_mean)'
            allocate(phymix3d_T_mean(I3DFIELD),stat=rc)
            if (rc /= 0) &
               stop 'calc_mean_fields.F90: Error allocating memory (nummix3d_T_mean)'
            allocate(phymix2d_T_mean(I2DFIELD),stat=rc)
            if (rc /= 0) &
               stop 'calc_mean_fields.F90: Error allocating memory (nummix2d_T_mean)'
         end if
         if (calc_salt) then
            allocate(nummix3d_S_mean(I3DFIELD),stat=rc)
            if (rc /= 0) &
               stop 'calc_mean_fields.F90: Error allocating memory (nummix3d_S_mean)'
            allocate(nummix2d_S_mean(I2DFIELD),stat=rc)
            if (rc /= 0) &
               stop 'calc_mean_fields.F90: Error allocating memory (nummix2d_S_mean)'
            allocate(phymix3d_S_mean(I3DFIELD),stat=rc)
            if (rc /= 0) &
               stop 'calc_mean_fields.F90: Error allocating memory (nummix3d_S_mean)'
            allocate(phymix2d_S_mean(I2DFIELD),stat=rc)
            if (rc /= 0) &
               stop 'calc_mean_fields.F90: Error allocating memory (nummix2d_S_mean)'
         end if
      end if
#endif
#ifdef GETM_BIO
      allocate(cc3dmean(numc,I3DFIELD),stat=rc)
      if (rc /= 0) &
          stop 'calc_mean_fields.F90: Error allocating memory (cc3dmean)'
#endif
      first = .false.
   end if

   if (step .eq. _ZERO_ ) then
      uumean=_ZERO_; vvmean=_ZERO_; wmean=_ZERO_
      humean=_ZERO_; hvmean=_ZERO_; hmean=_ZERO_
#ifndef NO_BAROCLINIC
      Tmean=_ZERO_; Smean=_ZERO_
      if (do_mixing_analysis) then
         if (calc_temp) then
            nummix3d_T_mean=_ZERO_; nummix2d_T_mean=_ZERO_
            phymix3d_T_mean=_ZERO_; phymix2d_T_mean=_ZERO_
         end if
         if (calc_salt) then
            nummix3d_S_mean=_ZERO_; nummix2d_S_mean=_ZERO_
            phymix3d_S_mean=_ZERO_; phymix2d_S_mean=_ZERO_
         end if
      end if
#endif
#ifdef GETM_BIO
      cc3dmean=_ZERO_
#endif
      ustarmean=_ZERO_; ustar2mean=_ZERO_; swrmean=_ZERO_
   end if

!  Sum every macro time step, even less would be okay
   if(mod(n,M) .eq. 0) then 

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
      Tmean = Tmean + T
      Smean = Smean + S
      if (do_mixing_analysis) then
         if (calc_temp) then
            nummix3d_T_mean = nummix3d_T_mean + nummix3d_T
            nummix2d_T_mean = nummix2d_T_mean + nummix2d_T
            phymix3d_T_mean = phymix3d_T_mean + phymix3d_T
            phymix2d_T_mean = phymix2d_T_mean + phymix2d_T
         end if
         if (calc_salt) then
            nummix3d_S_mean = nummix3d_S_mean + nummix3d_S
            nummix2d_S_mean = nummix2d_S_mean + nummix2d_S
            phymix3d_S_mean = phymix3d_S_mean + phymix3d_S
            phymix2d_S_mean = phymix2d_S_mean + phymix2d_S
         end if
      end if
#endif
#ifdef GETM_BIO
      if (bio_calc) cc3dmean=cc3dmean + cc3d
#endif
!  count them
      step = step + 1.0
   end if   ! here we summed them up

!  prepare for output
   if(meanout .gt. 0 .and. mod(n,meanout) .eq. 0) then

      if ( step .ge. 1.0) then
         uumean = uumean / step
         vvmean = vvmean / step
         wmean = wmean / step
         humean = humean / step
         hvmean = hvmean / step
         hmean = hmean / step

#ifndef NO_BAROCLINIC
         Tmean = Tmean / step
         Smean = Smean / step
         if (do_mixing_analysis) then
            if (calc_temp) then
               nummix3d_T_mean = nummix3d_T_mean / step
               nummix2d_T_mean = nummix2d_T_mean / step
               phymix3d_T_mean = phymix3d_T_mean / step
               phymix2d_T_mean = phymix2d_T_mean / step
            end if
            if (calc_salt) then
               nummix3d_S_mean = nummix3d_S_mean / step
               nummix2d_S_mean = nummix2d_S_mean / step
               phymix3d_S_mean = phymix3d_S_mean / step
               phymix2d_S_mean = phymix2d_S_mean / step
            end if
         end if
#endif
#ifdef GETM_BIO
         if (bio_calc) cc3dmean = cc3dmean / step
#endif
         ustarmean = ustarmean / step
         swrmean = swrmean / step

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
      step = _ZERO_
   end if

   call toc(TIM_CALCMEANF)
   return
   end subroutine calc_mean_fields
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2004 -  Adolf Stips  & Karsten Bolding                 !
!-----------------------------------------------------------------------
