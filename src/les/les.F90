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
   REALTYPE :: smag_const=0.28d0
!
! !PRIVATE DATA MEMBERS:
   integer,private,parameter :: SMAG_2D=1
   integer,private           :: les_method=SMAG_2D

!  explicit interface needed due to optional arguments
   interface
      subroutine les_smagorinsky(dudxC,dudxV,   &
#ifndef SLICE_MODEL
                                 dvdyC,dvdyU,   &
#endif
                                 shearX,shearU, &
                                 AmC,AmX,AmU,AmV)
         use domain, only: imin,imax,jmin,jmax
         IMPLICIT NONE
         REALTYPE,dimension(E2DFIELD),intent(in) :: dudxC,dudxV
#ifndef SLICE_MODEL
         REALTYPE,dimension(E2DFIELD),intent(in) :: dvdyC,dvdyU
#endif
         REALTYPE,dimension(E2DFIELD),intent(in) :: shearX,shearU
         REALTYPE,dimension(E2DFIELD),intent(out),optional :: AmC,AmX,AmU,AmV
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
   subroutine init_les(runtype)

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in) :: runtype
!
! !DESCRIPTION:
!
! Here, some necessary memory is allocated (in case of the compiler option
! {\tt STATIC}), and information is written to the log-file of
! the simulation.
!
! !LOCAL VARIABLES
   namelist /les/ les_method,smag_const
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_les() # ',Ncall
#endif

   LEVEL1 'init_les'

   call init_variables_les(runtype)

   if (les_mode .ne. NO_LES) then
      read(NAMLST,les)
      select case (les_method)
         case(SMAG_2D)
            LEVEL2 'Smagorinsky (1963) parameterisation'
            LEVEL3 'Smagorinsky constant: ',real(smag_const)
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
   subroutine do_les_2d(dudxC,dudxV, &
#ifndef SLICE_MODEL
                        dvdyC,dvdyU, &
#endif
                        shearX,shearU)
!
! !DESCRIPTION:
!
! !USES:
   use getm_timers, only: tic,toc,TIM_LES2D

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,dimension(E2DFIELD),intent(in) :: dudxC,dudxV
#ifndef SLICE_MODEL
   REALTYPE,dimension(E2DFIELD),intent(in) :: dvdyC,dvdyU
#endif
   REALTYPE,dimension(E2DFIELD),intent(in) :: shearX,shearU
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
      case(SMAG_2D)
         call les_smagorinsky(dudxC,dudxV,         &
#ifndef SLICE_MODEL
                              dvdyC,dvdyU,         &
#endif
                              shearX,shearU,       &
                              AmC=AmC_2d,AmX=AmX_2d)
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
   subroutine do_les_3d(dudxC,dudxV, &
#ifndef SLICE_MODEL
                        dvdyC,dvdyU, &
#endif
                        shearX,shearU)
!
! !DESCRIPTION:
!
!
! !USES:
   use getm_timers, only: tic,toc,TIM_LES3D

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE,dimension(I3DFIELD),intent(in) :: dudxC,dudxV
#ifndef SLICE_MODEL
   REALTYPE,dimension(I3DFIELD),intent(in) :: dvdyC,dvdyU
#endif
   REALTYPE,dimension(I3DFIELD),intent(in) :: shearX,shearU
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
      case(SMAG_2D)
         select case (les_mode)
            case(LES_MOMENTUM)
               do k=1,kmax
                  call les_smagorinsky(dudxC(:,:,k),dudxV(:,:,k),         &
#ifndef SLICE_MODEL
                                       dvdyC(:,:,k),dvdyU(:,:,k),         &
#endif
                                       shearX(:,:,k),shearU(:,:,k),       &
                                       AmC=AmC_3d(:,:,k),AmX=AmX_3d(:,:,k))
               end do
            case(LES_TRACER)
               do k=1,kmax
                  call les_smagorinsky(dudxC(:,:,k),dudxV(:,:,k),   &
#ifndef SLICE_MODEL
                                       dvdyC(:,:,k),dvdyU(:,:,k),   &
#endif
                                       shearX(:,:,k),shearU(:,:,k), &
                                       AmU=AmU_3d(:,:,k)            &
#ifndef SLICE_MODEL
                                       ,AmV=AmV_3d(:,:,k)           &
#endif
                                                                    )
               end do
            case(LES_BOTH)
               do k=1,kmax
                  call les_smagorinsky(dudxC(:,:,k),dudxV(:,:,k),           &
#ifndef SLICE_MODEL
                                       dvdyC(:,:,k),dvdyU(:,:,k),           &
#endif
                                       shearX(:,:,k),shearU(:,:,k),         &
                                       AmC=AmC_3d(:,:,k),AmX=AmX_3d(:,:,k), &
                                       AmU=AmU_3d(:,:,k)                    &
#ifndef SLICE_MODEL
                                       ,AmV=AmV_3d(:,:,k)                   &
#endif
                                                                            )
               end do
            case default
               FATAL 'No valid les_mode specified'
               stop 'do_les_3d()'
         end select
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
