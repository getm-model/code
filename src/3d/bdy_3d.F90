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
   use domain, only: nsbv,NWB,NNB,NEB,NSB,bdy_index
   use domain, only: wi,wfj,wlj,nj,nfi,nli,ei,efj,elj,sj,sfi,sli
   use variables_3d
#ifdef _FABM_
   use getm_fabm, only: fabm_calc,model,fabm_pel,fabm_ben
#endif
   IMPLICIT NONE
!
   private
!
! !PUBLIC DATA MEMBERS:
   public init_bdy_3d, do_bdy_3d
   REALTYPE, public, allocatable       :: S_bdy(:,:),T_bdy(:,:)
#ifdef _FABM_
   REALTYPE, public, allocatable       :: bio_bdy(:,:,:)
   integer, public, allocatable        :: have_bio_bdy_values(:)
#endif
   logical,  public                    :: bdy3d_tmrlx=.false.
   REALTYPE, public                    :: bdy3d_tmrlx_ucut=_ONE_/50
   REALTYPE, public                    :: bdy3d_tmrlx_max=_ONE_/4
   REALTYPE, public                    :: bdy3d_tmrlx_min=_ZERO_
!
! !PRIVATE DATA MEMBERS:
   REALTYPE,         allocatable       :: bdyvertS(:), bdyvertT(:)
   REALTYPE,         allocatable       :: rlxcoef(:)
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
   subroutine init_bdy_3d()
!
! !DESCRIPTION:
!
! Here, the necessary fields {\tt S\_bdy} and {\tt T\_bdy} for
! salinity and temperature, respectively, are allocated.
!
! !USES:
   IMPLICIT NONE
!
! !LOCAL VARIABLES:
   integer                   :: rc,i,j,k,n
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_bdy_3d() # ',Ncall
#endif

   LEVEL2 'init_bdy_3d()'
   allocate(S_bdy(0:kmax,nsbv),stat=rc)
   if (rc /= 0) stop 'init_init_bdy_3d: Error allocating memory (S_bdy)'

   allocate(T_bdy(0:kmax,nsbv),stat=rc)
   if (rc /= 0) stop 'init_init_bdy_3d: Error allocating memory (T_bdy)'

   allocate(bdyvertS(0:kmax),stat=rc)
   if (rc /= 0) stop 'init_init_bdy_3d: Error allocating memory (bdyvertS)'

   allocate(bdyvertT(0:kmax),stat=rc)
   if (rc /= 0) stop 'init_init_bdy_3d: Error allocating memory (bdyvertT)'

   allocate(rlxcoef(0:kmax),stat=rc)
   if (rc /= 0) stop 'init_init_bdy_3d: Error allocating memory (rlxcoef)'

#ifdef _FABM_
   npel=size(model%info%state_variables)
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
   integer                   :: i,j,k,l,n,o,ii,jj,kk
   REALTYPE                  :: sp(1:4),rat
   REALTYPE                  :: bdy3d_tmrlx_umin
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
! Sponge layer factors according to Martinsen and Engedahl, 1987.
   sp(1)=1.0
   sp(2)=0.5625
   sp(3)=0.25
   sp(4)=0.0625
! Hardcoding of lower limit of velocity cut-off for temporal relaxation.
! Linear variation between bdy3d_tmrlx_umin and bdy3d_tmrlx_ucut.
   bdy3d_tmrlx_umin = -_QUART_*bdy3d_tmrlx_ucut

   l = 0
   k = 0
   do n=1,NWB
      l = l+1
      k = bdy_index(l)
      i = wi(n)
      do j=wfj(n),wlj(n)
         if (bdy3d_tmrlx) then
            if (au(i,j).gt.0) then
!                Local temporal relaxation coeficient depends on
!                local current just *inside* domain:
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
            else
               rlxcoef(:)=bdy3d_tmrlx_max
            end if
!             Temporal relaxation: Weight inner (actual) solution near boundary
!             with boundary condition (outer solution.)
            wsum= MIN(az(i-1+2,j),1)*sp(2)                              &
                 +MIN(az(i-1+3,j),1)*sp(3)                              &
                 +MIN(az(i-1+4,j),1)*sp(4)
            if (wsum>_ZERO_) then
!                Get (weighted avr of) inner near-bdy solution (sponge) cells:
               bdyvertS(:) = (_ONE_/wsum)*(                              &
                               MIN(az(i-1+2,j),1) * sp(2) * S(i-1+2,j,:) &
                              +MIN(az(i-1+3,j),1) * sp(3) * S(i-1+3,j,:) &
                              +MIN(az(i-1+4,j),1) * sp(4) * S(i-1+4,j,:) &
                                          )
               bdyvertT(:) = (_ONE_/wsum)*(                              &
                               MIN(az(i-1+2,j),1) * sp(2) * T(i-1+2,j,:) &
                              +MIN(az(i-1+3,j),1) * sp(3) * T(i-1+3,j,:) &
                              +MIN(az(i-1+4,j),1) * sp(4) * T(i-1+4,j,:) &
                                          )
!                Weight inner and outer (bc) solutions for use
!                in spatial relaxation/sponge
               bdyvertS(:) = (_ONE_-rlxcoef(:))*bdyvertS(:) + rlxcoef(:)*S_bdy(:,k)
               bdyvertT(:) = (_ONE_-rlxcoef(:))*bdyvertT(:) + rlxcoef(:)*T_bdy(:,k)
            else
!                No near-bdy points. Just clamp bdy temporally:
               bdyvertS(:) = S_bdy(:,k)
               bdyvertT(:) = T_bdy(:,k)
            end if
         else
!            No time-relaxation. Just clamp at bondary points.
            bdyvertS(:) = S_bdy(:,k)
            bdyvertT(:) = T_bdy(:,k)
         end if
         do ii=1,4
!           Spatial relaxation / sponge:
            if (az(i-1+ii,j).gt.0) then
               S(i-1+ii,j,:) = sp(ii)*bdyvertS(:)+(_ONE_-sp(ii))*S(i-1+ii,j,:)
               T(i-1+ii,j,:) = sp(ii)*bdyvertT(:)+(_ONE_-sp(ii))*T(i-1+ii,j,:)
#ifdef _FABM_
               if (fabm_calc) then
                  do o=1,npel
                     if (have_bio_bdy_values(o) .eq. 1) then
                        fabm_pel(i-1+ii,j,:,o) = sp(ii)*bio_bdy(:,k,o) &
                                            +(1.-sp(ii))*fabm_pel(i-1+ii,j,:,o)
                     end if
                  end do
               end if
#endif
            end if
         end do
#ifdef _FABM_
!        zero gradient when we don't have bdy values
         if (fabm_calc) then
            do o=1,npel
               if (have_bio_bdy_values(o) .ne. 1) then
                  fabm_pel(i,j,:,o) = fabm_pel(i+1,j,:,o)
               end if
            end do
            fabm_ben(i,j,:) = fabm_ben(i+1,j,:)
         end if
#endif
         k = k+1
      end do
   end do

   do n = 1,NNB
      l = l+1
      k = bdy_index(l)
      j = nj(n)
      do i = nfi(n),nli(n)
         if (bdy3d_tmrlx) then
            if (av(i,j-1).gt.0) then
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
            else
               rlxcoef(:)=bdy3d_tmrlx_max
            end if
            wsum= MIN(az(i,j+1-2),1)*sp(2)                              &
                 +MIN(az(i,j+1-3),1)*sp(3)                              &
                 +MIN(az(i,j+1-4),1)*sp(4)
            if (wsum>_ZERO_) then
               bdyvertS(:) = (_ONE_/wsum)*(                              &
                               MIN(az(i,j+1-2),1) * sp(2) * S(i,j+1-2,:) &
                              +MIN(az(i,j+1-3),1) * sp(3) * S(i,j+1-3,:) &
                              +MIN(az(i,j+1-4),1) * sp(4) * S(i,j+1-4,:) &
                                          )
               bdyvertT(:) = (_ONE_/wsum)*(                              &
                               MIN(az(i,j+1-2),1) * sp(2) * T(i,j+1-2,:) &
                              +MIN(az(i,j+1-3),1) * sp(3) * T(i,j+1-3,:) &
                              +MIN(az(i,j+1-4),1) * sp(4) * T(i,j+1-4,:) &
                                          )
               bdyvertS(:) = (_ONE_-rlxcoef(:))*bdyvertS(:) + rlxcoef(:)*S_bdy(:,k)
               bdyvertT(:) = (_ONE_-rlxcoef(:))*bdyvertT(:) + rlxcoef(:)*T_bdy(:,k)
            else
               bdyvertS(:) = S_bdy(:,k)
               bdyvertT(:) = T_bdy(:,k)
            end if
         else
            bdyvertS(:) = S_bdy(:,k)
            bdyvertT(:) = T_bdy(:,k)
         end if
         do jj=1,4
            if (az(i,j+1-jj).gt.0) then
               S(i,j+1-jj,:) = sp(jj)*bdyvertS(:)+(_ONE_-sp(jj))*S(i,j+1-jj,:)
               T(i,j+1-jj,:) = sp(jj)*bdyvertT(:)+(_ONE_-sp(jj))*T(i,j+1-jj,:)
#ifdef _FABM_
               if (fabm_calc) then
                  do o=1,npel
                     if (have_bio_bdy_values(o) .eq. 1) then
                        fabm_pel(i,j+1-jj,:,o) = sp(jj)*bio_bdy(:,k,o) &
                                            +(1.-sp(jj))*fabm_pel(i,j+1-jj,:,o)
                     end if
                  end do
               end if
#endif
            end if
         end do
#ifdef _FABM_
!        zero gradient when we don't have bdy values
         if (fabm_calc) then
            do o=1,npel
               if (have_bio_bdy_values(o) .ne. 1) then
                  fabm_pel(i,j,:,o) = fabm_pel(i,j-1,:,o)
               end if
            end do
            fabm_ben(i,j,:) = fabm_ben(i,j-1,:)
         end if
#endif
         k = k+1
      end do
   end do

   do n=1,NEB
      l = l+1
      k = bdy_index(l)
      i = ei(n)
      do j=efj(n),elj(n)
         if (bdy3d_tmrlx) then
            if (au(i-1,j).gt.0) then
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
            else
               rlxcoef(:)=bdy3d_tmrlx_max
            end if
            wsum= MIN(az(i+1-2,j),1)*sp(2)                              &
                 +MIN(az(i+1-3,j),1)*sp(3)                              &
                 +MIN(az(i+1-4,j),1)*sp(4)
            if (wsum>_ZERO_) then
               bdyvertS(:) = (_ONE_/wsum)*(                              &
                               MIN(az(i+1-2,j),1) * sp(2) * S(i+1-2,j,:) &
                              +MIN(az(i+1-3,j),1) * sp(3) * S(i+1-3,j,:) &
                              +MIN(az(i+1-4,j),1) * sp(4) * S(i+1-4,j,:) &
                                          )
               bdyvertT(:) = (_ONE_/wsum)*(                              &
                               MIN(az(i+1-2,j),1) * sp(2) * T(i+1-2,j,:) &
                              +MIN(az(i+1-3,j),1) * sp(3) * T(i+1-3,j,:) &
                              +MIN(az(i+1-4,j),1) * sp(4) * T(i+1-4,j,:) &
                                          )
               bdyvertS(:) = (_ONE_-rlxcoef(:))*bdyvertS(:) + rlxcoef(:)*S_bdy(:,k)
               bdyvertT(:) = (_ONE_-rlxcoef(:))*bdyvertT(:) + rlxcoef(:)*T_bdy(:,k)
            else
               bdyvertS(:) = S_bdy(:,k)
               bdyvertT(:) = T_bdy(:,k)
            end if
         else
            bdyvertS(:) = S_bdy(:,k)
            bdyvertT(:) = T_bdy(:,k)
         end if
         do ii=1,4
            if (az(i+1-ii,j).gt.0) then
               S(i+1-ii,j,:) = sp(ii)*bdyvertS(:)+(_ONE_-sp(ii))*S(i+1-ii,j,:)
               T(i+1-ii,j,:) = sp(ii)*bdyvertT(:)+(_ONE_-sp(ii))*T(i+1-ii,j,:)
#ifdef _FABM_
               if (fabm_calc) then
                  do o=1,npel
                     if (have_bio_bdy_values(o) .eq. 1) then
                        fabm_pel(i+1-ii,j,:,o) = sp(ii)*bio_bdy(:,k,o) &
                                            +(1.-sp(ii))*fabm_pel(i+1-ii,j,:,o)
                     end if
                  end do
               end if
#endif
            end if
         end do
#ifdef _FABM_
!        zero gradient when we don't have bdy values
         if (fabm_calc) then
            do o=1,npel
               if (have_bio_bdy_values(o) .ne. 1) then
                  fabm_pel(i,j,:,o) = fabm_pel(i-1,j,:,o)
               end if
            end do
            fabm_ben(i,j,:) = fabm_ben(i-1,j,:)
         end if
#endif
         k = k+1
      end do
   end do

   do n = 1,NSB
      l = l+1
      k = bdy_index(l)
      j = sj(n)
      do i = sfi(n),sli(n)
         if (av(i,j).gt.0) then
            if (bdy3d_tmrlx) then
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
            else
               rlxcoef(:)=bdy3d_tmrlx_max
            end if
            wsum= MIN(az(i,j-1+2),1)*sp(2)                              &
                 +MIN(az(i,j-1+3),1)*sp(3)                              &
                 +MIN(az(i,j-1+4),1)*sp(4)
            if (wsum>_ZERO_) then
               bdyvertS(:) = (_ONE_/wsum)*(                              &
                               MIN(az(i,j-1+2),1) * sp(2) * S(i,j-1+2,:) &
                              +MIN(az(i,j-1+3),1) * sp(3) * S(i,j-1+3,:) &
                              +MIN(az(i,j-1+4),1) * sp(4) * S(i,j-1+4,:) &
                                          )
               bdyvertT(:) = (_ONE_/wsum)*(                              &
                               MIN(az(i,j-1+2),1) * sp(2) * T(i,j-1+2,:) &
                              +MIN(az(i,j-1+3),1) * sp(3) * T(i,j-1+3,:) &
                              +MIN(az(i,j-1+4),1) * sp(4) * T(i,j-1+4,:) &
                                          )
               bdyvertS(:) = (_ONE_-rlxcoef(:))*bdyvertS(:) + rlxcoef(:)*S_bdy(:,k)
               bdyvertT(:) = (_ONE_-rlxcoef(:))*bdyvertT(:) + rlxcoef(:)*T_bdy(:,k)
            else
               bdyvertS(:) = S_bdy(:,k)
               bdyvertT(:) = T_bdy(:,k)
            end if
         else
            bdyvertS(:) = S_bdy(:,k)
            bdyvertT(:) = T_bdy(:,k)
         end if
         do jj=1,4
            if (az(i,j-1+jj).gt.0) then
               S(i,j-1+jj,:) = sp(jj)*bdyvertS(:)+(_ONE_-sp(jj))*S(i,j-1+jj,:)
               T(i,j-1+jj,:) = sp(jj)*bdyvertT(:)+(_ONE_-sp(jj))*T(i,j-1+jj,:)
#ifdef _FABM_
               if (fabm_calc) then
                  do o=1,npel
                     if (have_bio_bdy_values(o) .eq. 1) then
                        fabm_pel(i,j-1+jj,:,o) = sp(jj)*bio_bdy(:,k,o) &
                                            +(1.-sp(jj))*fabm_pel(i,j-1+jj,:,o)
                     end if
                  end do
               end if
#endif
            end if
         end do
#ifdef _FABM_
!        zero gradient when we don't have bdy values
         if (fabm_calc) then
            do o=1,npel
               if (have_bio_bdy_values(o) .ne. 1) then
                  fabm_pel(i,j,:,o) = fabm_pel(i,j+1,:,o)
               end if
            end do
            fabm_ben(i,j,:) = fabm_ben(i,j+1,:)
         end if
#endif
         k = k+1
      end do
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

   end module bdy_3d

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Karsten Bolding and Hans Burchard               !
!-----------------------------------------------------------------------
