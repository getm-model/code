#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  bdy_3d - 3D boundary conditions \label{bdy-3d}
!
! !INTERFACE:
   module bdy_3d
!
! !DESCRIPTION:
!
! Here, the three-dimensional boundary
! conditions for temperature and salinity are handled.
!
! !USES:
   use halo_zones, only : H_TAG,U_TAG,V_TAG
   use domain, only: imin,jmin,imax,jmax,kmax,H,az,au,av
   use domain, only: nsbvl,nbdy,NWB,NNB,NEB,NSB,bdy_index,bdy_index_l
   use domain, only: bdy_3d_desc,bdy_3d_type
   use domain, only: need_3d_bdy
   use domain, only: wi,wfj,wlj,nj,nfi,nli,ei,efj,elj,sj,sfi,sli
   use variables_3d
#ifdef _FABM_
   use getm_fabm, only: fabm_calc,model,fabm_pel,fabm_ben
#endif
   use exceptions
   IMPLICIT NONE
!
   private
!
! !PUBLIC DATA MEMBERS:
   public init_bdy_3d, do_bdy_3d
   character(len=PATH_MAX),public      :: bdyfile_3d
   integer,public                      :: bdyfmt_3d
   integer,public                      :: bdy3d_ramp
   integer,public                      :: bdy3d_sponge_size=3
   logical,public                      :: bdy3d_tmrlx=.false.
   REALTYPE,public                     :: bdy3d_tmrlx_ucut=_ONE_/50
   REALTYPE                            :: bdy3d_tmrlx_umin
   REALTYPE,public                     :: bdy3d_tmrlx_max=_ONE_/4
   REALTYPE,public                     :: bdy3d_tmrlx_min=_ZERO_

   REALTYPE,dimension(:,:),allocatable,public :: bdy_data_S,bdy_data_T
#ifdef _FABM_
   REALTYPE, public, allocatable       :: bio_bdy(:,:,:)
   integer, public, allocatable        :: have_bio_bdy_values(:)
#endif
!
! !PRIVATE DATA MEMBERS:
   private bdy3d_active
   REALTYPE,         allocatable       :: bdyvertS(:), bdyvertT(:)
   REALTYPE,         allocatable       :: rlxcoef(:),sp(:)
#ifdef _FABM_
   integer                             :: npel=-1,nben=-1
#endif
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
! !IROUTINE: init_bdy_3d - initialising 3D boundary conditions
! \label{sec-init-bdy-3d}
!
! !INTERFACE:
   subroutine init_bdy_3d(bdy3d,runtype,hotstart)
!
! !DESCRIPTION:
!
! Here, the necessary fields {\tt S\_bdy} and {\tt T\_bdy} for
! salinity and temperature, respectively, are allocated.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: runtype
   logical, intent(in)                 :: hotstart
!
! !INPUT/OUTPUT PARAMETERS:
   logical, intent(inout)              :: bdy3d
!
! !LOCAL VARIABLES:
   integer                   :: rc,i,j,k,l,n
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_bdy_3d() # ',Ncall
#endif

   LEVEL2 'init_bdy_3d()'

   if (runtype.eq.2) bdy3d=.false.
   if (bdy3d .and. runtype.eq.3) then
      LEVEL3 'reset bdy3d=.false. in runtype=3'
      bdy3d = .false.
   end if

   if (bdy3d) then
      do l=1,nbdy
         if (bdy3d_active(bdy_3d_type(l))) then
            need_3d_bdy = .true.
            exit
         end if
      end do
      if (.not. need_3d_bdy) then
         bdy3d = .false.
      end if
   else
      do l=1,nbdy
         if (bdy3d_active(bdy_3d_type(l))) then
            LEVEL3 'bdy3d=F resets local 3D bdy #',l
            LEVEL4 'old: ',trim(bdy_3d_desc(bdy_3d_type(l)))
            bdy_3d_type(l) = CONSTANT
            LEVEL4 'new: ',trim(bdy_3d_desc(bdy_3d_type(l)))
         end if
      end do
   end if


   if (bdy3d) then

      LEVEL3 'bdyfile_3d=',TRIM(bdyfile_3d)
      LEVEL3 'bdyfmt_3d=',bdyfmt_3d
      if (bdy3d_ramp .gt. 1) then
         LEVEL3 'bdy3d_ramp=',bdy3d_ramp
         if (hotstart) then
            LEVEL4 'WARNING: hotstart is .true. AND bdy3d_ramp .gt. 1'
            LEVEL4 'WARNING: .. be sure you know what you are doing ..'
         end if
      end if

      if (bdy3d_tmrlx) then
         LEVEL3 'bdy3d_tmrlx=.true.'
         LEVEL3 'bdy3d_tmrlx_max=   ',bdy3d_tmrlx_max
         LEVEL3 'bdy3d_tmrlx_min=   ',bdy3d_tmrlx_min
         LEVEL3 'bdy3d_tmrlx_ucut=  ',bdy3d_tmrlx_ucut
         if (bdy3d_tmrlx_min<_ZERO_ .or. bdy3d_tmrlx_min>_ONE_)          &
              call getm_error("init_3d()",                               &
              "bdy3d_tmrlx_min is out of valid range [0:1]")
         if (bdy3d_tmrlx_max<bdy3d_tmrlx_min .or. bdy3d_tmrlx_min>_ONE_) &
              call getm_error("init_3d()",                               &
              "bdy3d_tmrlx_max is out of valid range [bdy3d_tmrlx_max:1]")
         if (bdy3d_tmrlx_ucut<_ZERO_)                                    &
              call getm_error("init_3d()",                               &
              "bdy3d_tmrlx_max is out of valid range [0:inf[")

!        Hardcoding of lower limit of velocity cut-off for temporal relaxation.
!        Linear variation between bdy3d_tmrlx_umin and bdy3d_tmrlx_ucut.
         bdy3d_tmrlx_umin = -_QUART_*bdy3d_tmrlx_ucut
         LEVEL3 'bdy3d_tmrlx_umin=  ',bdy3d_tmrlx_umin
      end if

      allocate(bdy_data_S(0:kmax,nsbvl),stat=rc)
      if (rc /= 0) stop 'init_init_bdy_3d: Error allocating memory (bdy_data_S)'

      allocate(bdy_data_T(0:kmax,nsbvl),stat=rc)
      if (rc /= 0) stop 'init_init_bdy_3d: Error allocating memory (bdy_data_T)'

      allocate(bdyvertS(0:kmax),stat=rc)
      if (rc /= 0) stop 'init_init_bdy_3d: Error allocating memory (bdyvertS)'

      allocate(bdyvertT(0:kmax),stat=rc)
      if (rc /= 0) stop 'init_init_bdy_3d: Error allocating memory (bdyvertT)'

      allocate(rlxcoef(0:kmax),stat=rc)
      if (rc /= 0) stop 'init_init_bdy_3d: Error allocating memory (rlxcoef)'

      if (bdy3d_sponge_size .gt. 0) then
         allocate(sp(bdy3d_sponge_size),stat=rc)
         if (rc /= 0) stop 'init_init_bdy_3d: Error allocating memory (sp)'

!        Sponge layer factors according to Martinsen and Engedahl, 1987.
!        Note (KK): factor=1 (bdy cell) does not count for sponge size
!                   (in contrast to earlier GETM)
         LEVEL3 "sponge layer factors:"
         do i=1,bdy3d_sponge_size
            sp(i) = ((_ONE_+bdy3d_sponge_size-i)/(_ONE_+bdy3d_sponge_size))**2
            LEVEL4 "sp(",i,")=",real(sp(i))
         end do
      end if
   end if

#ifdef _FABM_
   if (fabm_calc) then
      npel=size(model%info%state_variables)
   end if
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving init_bdy_3d()'
   write(debug,*)
#endif
   return
   end subroutine init_bdy_3d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_bdy_3d  - updating 3D boundary conditions
! \label{sec-do-bdy-3d}
!
! !INTERFACE:
   subroutine do_bdy_3d(tag,field)
!
! !DESCRIPTION:
!
! Here, the boundary conditions for salinity and temperature are
! copied to the boundary points and relaxed to the near boundary points
! by means of the flow relaxation scheme by \cite{MARTINSENea87}.
!
! As an extention to the flow relaxation scheme, it is possible
! to relax the boundary point values to the specified boundary
! condition in time, thus giving more realistic situations
! especially for outgoing flow conditions. This nudging is implemented
! to depend on the local (3D) current velocity perpendicular to
! the boundary. For strong outflow, the boundary condition is turned
! off, while for inflows it is given a high impact.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: tag
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout)             :: field(I3DFIELD)
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,kl,l,n,o,ii,jj,kk
   REALTYPE                  :: rat
   REALTYPE                  :: wsum
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_bdy_3d() # ',Ncall
#endif

#if 0
   select case (tag)
      case (1)
!       Lateral zero-gradient boundary condition (north & south)
         do k=1,kmax
            do i=imin,imax
               if (au(i,jmin) .eq. 3) field(i,jmin,k)=field(i,jmin+1,k)
               if (au(i,jmax) .eq. 3) field(i,jmax,k)=field(i,jmax-1,k)
            end do
         end do
      case (2)
!       Lateral zero-gradient boundary conditions (west & east)
         do k=1,kmax
            do j=jmin,jmax
               if (av(imin,j) .eq. 3) field(imin,j,k)=field(imin+1,j,k)
               if (av(imax,j) .eq. 3) field(imax,j,k)=field(imax-1,j,k)
            end do
         end do
      case default
         FATAL 'Non valid tag'
         stop 'do_bdy_3d'
   end select
#endif

#ifndef NO_BAROCLINIC

   l = 0
   do n=1,NWB
      l = l+1
      k = bdy_index(l)
      kl = bdy_index_l(l)
      i = wi(n)
      select case (bdy_3d_type(l))
         case(ZERO_GRADIENT)
            do j=wfj(n),wlj(n)
               S(i,j,:) = S(i+1,j,:)
               T(i,j,:) = T(i+1,j,:)
            end do
         case(CLAMPED)
            do j=wfj(n),wlj(n)
               if (az(i,j) .eq. 2) then
                  if (bdy3d_tmrlx) then
!                    Temporal relaxation: Weight inner (actual) solution near boundary
!                    with boundary condition (outer solution.)
                     wsum = _ZERO_
                     bdyvertS(:) = _ZERO_
                     bdyvertT(:) = _ZERO_
                     do ii=1,bdy3d_sponge_size
!                       Get (weighted avr of) inner near-bdy solution (sponge) cells:
                        if(az(i+ii,j) .ne. 0) then
                           wsum = wsum + sp(ii)
                           bdyvertS(:) = bdyvertS(:) + sp(ii)*S(i+ii,j,:)
                           bdyvertT(:) = bdyvertT(:) + sp(ii)*T(i+ii,j,:)
                        else
                           exit
                        end if
                     end do
                     if (wsum>_ZERO_) then
!                       Local temporal relaxation coeficient depends on
!                       local current just *inside* domain:
                        do kk=1,kmax
                           if (uu(i,j,kk).ge.bdy3d_tmrlx_ucut) then
                              rlxcoef(kk) = bdy3d_tmrlx_max
                           else if (uu(i,j,kk).le.bdy3d_tmrlx_umin) then
                              rlxcoef(kk) = bdy3d_tmrlx_min
                           else
                              rlxcoef(kk) = (bdy3d_tmrlx_max-bdy3d_tmrlx_min)    &
                                   *(uu(i,j,kk)-bdy3d_tmrlx_umin)                &
                                   /(bdy3d_tmrlx_ucut-bdy3d_tmrlx_umin)          &
                                   + bdy3d_tmrlx_min
                           end if
                        end do
!                       Weight inner and outer (bc) solutions for use
!                       in spatial relaxation/sponge
                        bdyvertS(:) = (_ONE_-rlxcoef(:))*bdyvertS(:)/wsum + rlxcoef(:)*bdy_data_S(:,kl)
                        bdyvertT(:) = (_ONE_-rlxcoef(:))*bdyvertT(:)/wsum + rlxcoef(:)*bdy_data_T(:,kl)
                     else
   !                    No near-bdy points. Just clamp bdy temporally:
                        bdyvertS(:) = bdy_data_S(:,kl)
                        bdyvertT(:) = bdy_data_T(:,kl)
                     end if
                  else
!                    No time-relaxation. Just clamp at bondary points.
                     bdyvertS(:) = bdy_data_S(:,kl)
                     bdyvertT(:) = bdy_data_T(:,kl)
                  end if
                  S(i,j,:) = bdyvertS(:)
                  T(i,j,:) = bdyvertT(:)
#ifdef _FABM_
                  if (fabm_calc) then
                     do o=1,npel
                        if (have_bio_bdy_values(o) .eq. 1) then
                           fabm_pel(i,j,:,o) = bio_bdy(:,k,o)
                        else
!                          zero gradient when we don't have bdy values
                           fabm_pel(i,j,:,o) = fabm_pel(i+1,j,:,o)
                        end if
                     end do
                     fabm_ben(i,j,:) = fabm_ben(i+1,j,:)
                  end if
#endif
                  do ii=1,bdy3d_sponge_size
!                    Spatial relaxation / sponge:
                     if (az(i+ii,j) .ne. 0) then
                        S(i+ii,j,:) = sp(ii)*bdyvertS(:)+(_ONE_-sp(ii))*S(i+ii,j,:)
                        T(i+ii,j,:) = sp(ii)*bdyvertT(:)+(_ONE_-sp(ii))*T(i+ii,j,:)
#ifdef _FABM_
                        if (fabm_calc) then
                           do o=1,npel
                              if (have_bio_bdy_values(o) .eq. 1) then
                                 fabm_pel(i+ii,j,:,o) = sp(ii)*bio_bdy(:,k,o) &
                                                     +(_ONE_-sp(ii))*fabm_pel(i+ii,j,:,o)
                              end if
                           end do
                        end if
#endif
                     else
                        exit
                     end if
                  end do
               end if
               k = k+1
               kl = kl + 1
            end do
      end select
   end do

   do n = 1,NNB
      l = l+1
      k = bdy_index(l)
      kl = bdy_index_l(l)
      j = nj(n)
      select case (bdy_3d_type(l))
         case(ZERO_GRADIENT)
            do i = nfi(n),nli(n)
               S(i,j,:) = S(i,j-1,:)
               T(i,j,:) = T(i,j-1,:)
            end do
         case(CLAMPED)
            do i = nfi(n),nli(n)
               if (az(i,j) .eq. 2) then
                  if (bdy3d_tmrlx) then
                     wsum = _ZERO_
                     bdyvertS(:) = _ZERO_
                     bdyvertT(:) = _ZERO_
                     do jj=1,bdy3d_sponge_size
                        if(az(i,j-jj)  .ne. 0) then
                           wsum = wsum + sp(jj)
                           bdyvertS(:) = bdyvertS(:) + sp(jj)*S(i,j-jj,:)
                           bdyvertT(:) = bdyvertT(:) + sp(jj)*T(i,j-jj,:)
                        else
                           exit
                        end if
                     end do
                     if (wsum>_ZERO_) then
                        do kk=1,kmax
                           if (vv(i,j-1,kk).le.-bdy3d_tmrlx_ucut) then
                              rlxcoef(kk) = bdy3d_tmrlx_max
                           else if (vv(i,j-1,kk).ge.-bdy3d_tmrlx_umin) then
                              rlxcoef(kk) = bdy3d_tmrlx_min
                           else
                              rlxcoef(kk) = -(bdy3d_tmrlx_max-bdy3d_tmrlx_min)   &
                                   *(vv(i,j-1,kk)+bdy3d_tmrlx_umin)              &
                                   /(bdy3d_tmrlx_ucut-bdy3d_tmrlx_umin)          &
                                   + bdy3d_tmrlx_min
                           end if
                        end do
                        bdyvertS(:) = (_ONE_-rlxcoef(:))*bdyvertS(:)/wsum + rlxcoef(:)*bdy_data_S(:,kl)
                        bdyvertT(:) = (_ONE_-rlxcoef(:))*bdyvertT(:)/wsum + rlxcoef(:)*bdy_data_T(:,kl)
                     else
                        bdyvertS(:) = bdy_data_S(:,kl)
                        bdyvertT(:) = bdy_data_T(:,kl)
                     end if
                  else
                     bdyvertS(:) = bdy_data_S(:,kl)
                     bdyvertT(:) = bdy_data_T(:,kl)
                  end if
                  S(i,j,:) = bdyvertS(:)
                  T(i,j,:) = bdyvertT(:)
#ifdef _FABM_
                  if (fabm_calc) then
                     do o=1,npel
                        if (have_bio_bdy_values(o) .eq. 1) then
                           fabm_pel(i,j,:,o) = bio_bdy(:,k,o)
                        else
!                          zero gradient when we don't have bdy values
                           fabm_pel(i,j,:,o) = fabm_pel(i,j-1,:,o)
                        end if
                     end do
                     fabm_ben(i,j,:) = fabm_ben(i,j-1,:)
                  end if
#endif
                  do jj=1,bdy3d_sponge_size
                     if (az(i,j-jj) .ne. 0) then
                        S(i,j-jj,:) = sp(jj)*bdyvertS(:)+(_ONE_-sp(jj))*S(i,j-jj,:)
                        T(i,j-jj,:) = sp(jj)*bdyvertT(:)+(_ONE_-sp(jj))*T(i,j-jj,:)
#ifdef _FABM_
                        if (fabm_calc) then
                           do o=1,npel
                              if (have_bio_bdy_values(o) .eq. 1) then
                                 fabm_pel(i,j-jj,:,o) = sp(jj)*bio_bdy(:,k,o) &
                                                     +(_ONE_-sp(jj))*fabm_pel(i,j-jj,:,o)
                              end if
                           end do
                        end if
#endif
                     else
                        exit
                     end if
                  end do
               end if
               k = k+1
               kl = kl + 1
            end do
      end select
   end do

   do n=1,NEB
      l = l+1
      k = bdy_index(l)
      kl = bdy_index_l(l)
      i = ei(n)
      select case (bdy_3d_type(l))
         case(ZERO_GRADIENT)
            do j=efj(n),elj(n)
               S(i,j,:) = S(i-1,j,:)
               T(i,j,:) = T(i-1,j,:)
            end do
         case(CLAMPED)
            do j=efj(n),elj(n)
               if (az(i,j) .eq. 2) then
                  if (bdy3d_tmrlx) then
                     wsum = _ZERO_
                     bdyvertS(:) = _ZERO_
                     bdyvertT(:) = _ZERO_
                     do ii=1,bdy3d_sponge_size
                        if(az(i-ii,j) .ne. 0) then
                           wsum = wsum + sp(ii)
                           bdyvertS(:) = bdyvertS(:) + sp(ii)*S(i-ii,j,:)
                           bdyvertT(:) = bdyvertT(:) + sp(ii)*T(i-ii,j,:)
                        else
                           exit
                        end if
                     end do
                     if (wsum>_ZERO_) then
                        do kk=1,kmax
                           if (uu(i-1,j,kk).le.-bdy3d_tmrlx_ucut) then
                              rlxcoef(kk) = bdy3d_tmrlx_max
                           else if (uu(i-1,j,kk).ge.-bdy3d_tmrlx_umin) then
                              rlxcoef(kk) = bdy3d_tmrlx_min
                           else
                              rlxcoef(kk) = -(bdy3d_tmrlx_max-bdy3d_tmrlx_min)   &
                                   *(uu(i-1,j,kk)+bdy3d_tmrlx_umin)              &
                                   /(bdy3d_tmrlx_ucut-bdy3d_tmrlx_umin)          &
                                   + bdy3d_tmrlx_min
                           end if
                        end do
                        bdyvertS(:) = (_ONE_-rlxcoef(:))*bdyvertS(:)/wsum + rlxcoef(:)*bdy_data_S(:,kl)
                        bdyvertT(:) = (_ONE_-rlxcoef(:))*bdyvertT(:)/wsum + rlxcoef(:)*bdy_data_T(:,kl)
                     else
                        bdyvertS(:) = bdy_data_S(:,kl)
                        bdyvertT(:) = bdy_data_T(:,kl)
                     end if
                  else
                     bdyvertS(:) = bdy_data_S(:,kl)
                     bdyvertT(:) = bdy_data_T(:,kl)
                  end if
                  S(i,j,:) = bdyvertS(:)
                  T(i,j,:) = bdyvertT(:)
#ifdef _FABM_
                  if (fabm_calc) then
                     do o=1,npel
                        if (have_bio_bdy_values(o) .eq. 1) then
                           fabm_pel(i,j,:,o) = bio_bdy(:,k,o)
                        else
!                          zero gradient when we don't have bdy values
                           fabm_pel(i,j,:,o) = fabm_pel(i-1,j,:,o)
                        end if
                     end do
                     fabm_ben(i,j,:) = fabm_ben(i-1,j,:)
                  end if
#endif
                  do ii=1,bdy3d_sponge_size
                     if (az(i-ii,j) .ne. 0) then
                        S(i-ii,j,:) = sp(ii)*bdyvertS(:)+(_ONE_-sp(ii))*S(i-ii,j,:)
                        T(i-ii,j,:) = sp(ii)*bdyvertT(:)+(_ONE_-sp(ii))*T(i-ii,j,:)
#ifdef _FABM_
                        if (fabm_calc) then
                           do o=1,npel
                              if (have_bio_bdy_values(o) .eq. 1) then
                                 fabm_pel(i-ii,j,:,o) = sp(ii)*bio_bdy(:,k,o) &
                                                     +(_ONE_-sp(ii))*fabm_pel(i-ii,j,:,o)
                              end if
                           end do
                        end if
#endif
                     else
                        exit
                     end if
                  end do
               end if
               k = k+1
               kl = kl + 1
            end do
      end select
   end do

   do n = 1,NSB
      l = l+1
      k = bdy_index(l)
      kl = bdy_index_l(l)
      j = sj(n)
      select case (bdy_3d_type(l))
         case(ZERO_GRADIENT)
            do i = sfi(n),sli(n)
               S(i,j,:) = S(i,j+1,:)
               T(i,j,:) = T(i,j+1,:)
            end do
         case(CLAMPED)
            do i = sfi(n),sli(n)
               if (az(i,j) .eq. 2) then
                  if (bdy3d_tmrlx) then
                     wsum = _ZERO_
                     bdyvertS(:) = _ZERO_
                     bdyvertT(:) = _ZERO_
                     do jj=1,bdy3d_sponge_size
                        if(az(i,j+jj) .ne. 0) then
                           wsum = wsum + sp(jj)
                           bdyvertS(:) = bdyvertS(:) + sp(jj)*S(i,j+jj,:)
                           bdyvertS(:) = bdyvertT(:) + sp(jj)*T(i,j+jj,:)
                        else
                           exit
                        end if
                     end do
                     if (wsum>_ZERO_) then
                        do kk=1,kmax
                           if (vv(i,j,kk).ge.bdy3d_tmrlx_ucut) then
                              rlxcoef(kk) = bdy3d_tmrlx_max
                           else if (vv(i,j,kk).le.bdy3d_tmrlx_umin) then
                              rlxcoef(kk) = bdy3d_tmrlx_min
                           else
                              rlxcoef(kk) = (bdy3d_tmrlx_max-bdy3d_tmrlx_min)    &
                                   *(vv(i,j,kk)-bdy3d_tmrlx_umin)                &
                                   /(bdy3d_tmrlx_ucut-bdy3d_tmrlx_umin)          &
                                   + bdy3d_tmrlx_min
                           end if
                        end do
                        bdyvertS(:) = (_ONE_-rlxcoef(:))*bdyvertS(:)/wsum + rlxcoef(:)*bdy_data_S(:,kl)
                        bdyvertT(:) = (_ONE_-rlxcoef(:))*bdyvertT(:)/wsum + rlxcoef(:)*bdy_data_T(:,kl)
                     else
                        bdyvertS(:) = bdy_data_S(:,kl)
                        bdyvertT(:) = bdy_data_T(:,kl)
                     end if
                  else
                     bdyvertS(:) = bdy_data_S(:,kl)
                     bdyvertT(:) = bdy_data_T(:,kl)
                  end if
                  S(i,j,:) = bdyvertS(:)
                  T(i,j,:) = bdyvertT(:)
#ifdef _FABM_
                  if (fabm_calc) then
                     do o=1,npel
                        if (have_bio_bdy_values(o) .eq. 1) then
                           fabm_pel(i,j,:,o) = bio_bdy(:,k,o)
                        else
!                          zero gradient when we don't have bdy values
                           fabm_pel(i,j,:,o) = fabm_pel(i,j+1,:,o)
                        end if
                     end do
                     fabm_ben(i,j,:) = fabm_ben(i,j+1,:)
                  end if
#endif
                  do jj=1,bdy3d_sponge_size
                     if (az(i,j+jj) .ne. 0) then
                        S(i,j+jj,:) = sp(jj)*bdyvertS(:)+(_ONE_-sp(jj))*S(i,j+jj,:)
                        T(i,j+jj,:) = sp(jj)*bdyvertT(:)+(_ONE_-sp(jj))*T(i,j+jj,:)
#ifdef _FABM_
                        if (fabm_calc) then
                           do o=1,npel
                              if (have_bio_bdy_values(o) .eq. 1) then
                                 fabm_pel(i,j+jj,:,o) = sp(jj)*bio_bdy(:,k,o) &
                                                     +(_ONE_-sp(jj))*fabm_pel(i,j+jj,:,o)
                              end if
                           end do
                        end if
#endif
                     else
                        exit
                     end if
                  end do
               end if
               k = k+1
               kl = kl + 1
            end do
      end select
   end do

#ifdef _FABM_
   if ( allocated(fabm_pel) ) then
      do n=1,size(fabm_pel,4)
         call mirror_bdy_3d(fabm_pel(:,:,:,n),H_TAG)
      end do
   end if
  if ( allocated(fabm_ben) ) then
      do n=1, size(fabm_ben,3)
         call mirror_bdy_3d(fabm_ben(:,:,  n),H_TAG)
      end do
  end if
#endif
#endif

#ifdef DEBUG
   write(debug,*) 'leaving do_bdy_3d()'
   write(debug,*)
#endif
   return
   end subroutine do_bdy_3d
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  LOGICAL function bdy3d_active -
!
! !INTERFACE:
   logical function bdy3d_active(type_3d)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)  :: type_3d
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!EOP
!-----------------------------------------------------------------------
!BOC

   select case (type_3d)
      case (CONSTANT)
         bdy3d_active = .false.
      case (CLAMPED)
         bdy3d_active = .true.
      case (ZERO_GRADIENT)
         bdy3d_active = .false.
      case default
         bdy3d_active = .false.
   end select

   return
   end function bdy3d_active
!EOC
!-----------------------------------------------------------------------

   end module bdy_3d

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Karsten Bolding and Hans Burchard               !
!-----------------------------------------------------------------------
