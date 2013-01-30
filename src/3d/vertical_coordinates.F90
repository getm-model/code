#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: vertical_coordinates
!
! !INTERFACE:
   module vertical_coordinates
!
! !DESCRIPTION:
!
! !USES:
#ifdef SLICE_MODEL
   use domain, only: imin,imax,jmin,jmax
   use variables_3d, only: kvmin,hvo,hvn
#endif
   use domain, only: kmax,vert_cord,maxdepth
   use exceptions
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   public coordinates
   logical,public  :: restart_with_ho=.false.
   logical,public  :: restart_with_hn=.false.
   REALTYPE,public :: cord_relax=_ZERO_
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister & Knut Klingbeil
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:  coordinates - defines the vertical coordinate
! \label{sec-coordinates}
!
! !INTERFACE:
   subroutine coordinates(hotstart)
!
! !DESCRIPTION:
!
! Here, the vertical layer distribution in T-, U- and V-points is updated
! during every macro time step. This is done for the old and the new
! layer thicknesses at every point. Calculation of the layer distribution
! in the U- and V-points is done indepently from the calculation in the
! T-points, since different methods for the calculation of the
! bathymetry values in the U- and V-points are possible, see routine
! {\tt uv\_depths} described on page \pageref{sec-uv-depth}.
!
! The different methods for the vertical layer distribution
! are initialised and called to be chosen by the namelist paramter {\tt vert\_cord}:\\
! \\
! {\tt vert\_cord=1}: sigma coordinates (section~\ref{sec-sigma-coordinates}) \\
! {\tt vert\_cord=2}: z-level (not coded yet) \\
! {\tt vert\_cord=3}: general vertical coordinates (gvc, section~\ref{sec-general-coordinates})
! \\
! {\tt vert\_cord=5}: adaptive vertical coordinates (section~\ref{sec-adaptive-coordinates}) \\
! \\
!
!
! !USES:
   use getm_timers, only: tic, toc,TIM_COORDS
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!   integer, intent(in)                 :: cord_type
!   REALTYPE, intent(in)                :: cord_relax
!   REALTYPE, intent(in)                :: maxdepth
   logical, intent(in)                 :: hotstart
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   logical, save   :: first=.true.
   integer         :: ii
!   integer         :: preadapt=0
#ifdef SLICE_MODEL
   integer          :: i,j,k
#endif
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'coordinates() # ',Ncall
#endif
   call tic(TIM_COORDS)

   if (first) then
      select case (vert_cord)
         case (_SIGMA_COORDS_) ! sigma coordinates
            LEVEL2 'using ',kmax,' sigma layers'
            call sigma_coordinates(.true.)
         case (_Z_COORDS_) ! z-level
            call getm_error("coordinates()","z-levels not implemented yet")
         case (_GENERAL_COORDS_) ! general vertical coordinates
            LEVEL2 'using ',kmax,' gvc layers'
            call general_coordinates(.true.,cord_relax,maxdepth)
         case (_HYBRID_COORDS_) ! hybrid vertical coordinates
            LEVEL2 'using ',kmax,' hybrid layers'
            call hybrid_coordinates(.true.)
STDERR 'coordinates(): hybrid_coordinates not coded yet'
stop
         case (_ADAPTIVE_COORDS_) ! adaptive vertical coordinates
            LEVEL2 'using ',kmax,' adaptive layers'
            call adaptive_coordinates(.true.,hotstart)
         case default
      end select
      first = .false.
   else
      select case (vert_cord)
         case (_SIGMA_COORDS_) ! sigma coordinates
            call sigma_coordinates(.false.)
         case (_Z_COORDS_) ! z-level
         case (_GENERAL_COORDS_) ! general vertical coordinates
            call general_coordinates(.false.,cord_relax,maxdepth)
         case (_HYBRID_COORDS_) ! hybrid vertical coordinates
            call hybrid_coordinates(.false.)
         case (_ADAPTIVE_COORDS_) ! adaptive vertical coordinates
            call adaptive_coordinates(.false.,hotstart)
         case default
      end select
   end if ! first

#ifdef SLICE_MODEL
   do i=imin,imax
      do k=kvmin(i,2),kmax
         hvo(i,1,k)=hvo(i,2,k)
         hvo(i,3,k)=hvo(i,2,k)
         hvn(i,1,k)=hvn(i,2,k)
         hvn(i,3,k)=hvn(i,2,k)
      end do
   end do
#endif

   call toc(TIM_COORDS)
#ifdef DEBUG
   write(debug,*) 'Leaving coordinates()'
   write(debug,*)
#endif
   return
   end subroutine coordinates
!EOC
!-----------------------------------------------------------------------

   end module vertical_coordinates

!-----------------------------------------------------------------------
! Copyright (C) 2012 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
