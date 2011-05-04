!$Id: les.F90,v 1.23 2009-09-30 11:28:45 bjb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: les
!
! !INTERFACE:
   module les
!
! !DESCRIPTION:
!
!
! !USES:
   use variables_les
   use exceptions

   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   public init_les, do_les_2d
#ifndef NO_3D
   public do_les_3d
#endif
   integer  :: les_method=1
   REALTYPE :: smag_const=0.28d0
!
! !PRIVATE DATA MEMBERS:
   integer, private, parameter         :: SMAGORINSKY=1

!  explicit interface needed due to optional arguments
   interface
      subroutine les_smagorinsky(U,V,DU,DV,Am,AmX,AmU,AmV)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,intent(in)           :: U(E2DFIELD),V(E2DFIELD)
         REALTYPE,intent(in)           :: DU(E2DFIELD),DV(E2DFIELD)
         REALTYPE,intent(out)          :: Am(E2DFIELD),AmX(E2DFIELD)
         REALTYPE,intent(out),optional :: AmU(E2DFIELD),AmV(E2DFIELD)
      end subroutine les_smagorinsky
   end interface
!
! !REVISION HISTORY:
!  Original author(s): Knut Klingbeil
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_les - initialising LES
! \label{sec-init-les}
!
! !INTERFACE:
   subroutine init_les(runtype,Am_method,Am_const)

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: runtype
   integer, intent(in)                 :: Am_method
   REALTYPE, intent(in)                :: Am_const
!
! !DESCRIPTION:
!
! Here, some necessary memory is allocated (in case of the compiler option
! {\tt STATIC}), and information is written to the log-file of
! the simulation.
!
! !LOCAL VARIABLES
   namelist /les/ les_method,smag_const
   integer         :: rc
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_les() # ',Ncall
#endif

   if (Am_method .eq. 2) LEVEL1 'init_les()'

!  Allocates memory for the public data members - if not static
   call init_variables_les(runtype,Am_method,Am_const)

   if (Am_method .eq. 2) then
      read(NAMLST,nml=les)
      select case (les_method)
         case(SMAGORINSKY)
            LEVEL2 'Smagorinsky (1963) parameterisation'
            LEVEL3 'Smagorinsky constant: ',smag_const
         case default
            FATAL 'No valid les_method specified'
            stop 'init_internal_pressure()'
      end select
   end if

   return
   end subroutine init_les
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_les_2d - 2D LES parameterisation
! \label{sec-do-les_2d}
!
! !INTERFACE:
   subroutine do_les_2d(U,V,DU,DV)
!
! !DESCRIPTION:
!
! !USES:
   use getm_timers, only: tic,toc,TIM_LES2D

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in) :: U(E2DFIELD),V(E2DFIELD)
   REALTYPE, intent(in) :: DU(E2DFIELD),DV(E2DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !LOCAL VARIABLES:

!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_les_2d() # ',Ncall
#endif

   call tic(TIM_LES2D)
   select case (les_method)
      case(SMAGORINSKY)
         call les_smagorinsky(U,V,DU,DV,Am2d,AmX2d)
      case default
         FATAL 'No valid les_method specified'
         stop 'do_les_2d()'
   end select
   call toc(TIM_LES2D)

#ifdef DEBUG
   write(debug,*) 'Leaving do_les_2d()'
   write(debug,*)
#endif
   return
   end subroutine do_les_2d
!EOC
!-----------------------------------------------------------------------
#ifndef NO_3D
!BOP
!
! !IROUTINE:  do_les_3d - 3D LES parameterisation
! \label{sec-do-les_3d}
!
! !INTERFACE:
   subroutine do_les_3d(uu,vv,hu,hv)
!
! !DESCRIPTION:
!
! !USES:
   use getm_timers, only: tic,toc,TIM_LES2D

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in) :: uu(I3DFIELD),vv(I3DFIELD)
   REALTYPE, intent(in) :: hu(I3DFIELD),hv(I3DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !LOCAL VARIABLES:
   integer :: k

!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_les_3d() # ',Ncall
#endif

   call tic(TIM_LES3D)
   select case (les_method)
      case(SMAGORINSKY)
         do k=1,kmax
            call les_smagorinsky(uu(:,:,k),vv(:,:,k),hu(:,:,k),hv(:,:,k), &
                                 Am3d(:,:,k),AmX3d(:,:,k), &
                                 AmU3d(:,:,k),AmV3d(:,:,k))
         end do
      case default
         FATAL 'No valid les_method specified'
         stop 'do_les_3d()'
   end select
   call toc(TIM_LES3D)

#ifdef DEBUG
   write(debug,*) 'Leaving do_les_3d()'
   write(debug,*)
#endif
   return
   end subroutine do_les_3d
!EOC

!-----------------------------------------------------------------------
#endif

   end module les

!-----------------------------------------------------------------------
! Copyright (C) 2011 - Knut Klingbeil                                  !
!-----------------------------------------------------------------------
