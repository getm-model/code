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

!
! !USES:
   use domain    , only: ill,ihl,ilg,ihg,jll,jhl,jlg,jhg
   use domain    , only: imin,imax,jmin,jmax
   use domain    , only: az
   use halo_zones, only: update_2d_halo,wait_halo,H_TAG

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
   integer                 :: i,j
   REALTYPE                :: smag_const=0.28d0
   character(len=PATH_MAX) :: smag_file
   namelist /les/ les_method,smag_method,smag_const,smag_file
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_les() # ',Ncall
#endif

   LEVEL1 'init_les'

   if (les_mode .eq. NO_LES) return

   read(NAMLST,les)
   call init_variables_les(runtype)

   select case(les_method)

      case (SMAG_2D)

         LEVEL2 'Smagorinsky (1963) parameterisation'

         select case(smag_method)

            case (SMAG_CONSTANT)

               LEVEL3 'Smagorinsky parameter: ',real(smag_const)
               SmagC2_2d = smag_const*smag_const
               SmagX2_2d => SmagC2_2d
               SmagU2_2d => SmagC2_2d
               SmagV2_2d => SmagC2_2d

            case (SMAG_FROMFILE)

               LEVEL3 'getting Smagorinsky parameter from ',trim(smag_file)
               call get_2d_field(trim(smag_file),"smag2d",ilg,ihg,jlg,jhg,.true.,SmagC2_2d(ill:ihl,jll:jhl))
               SmagC2_2d = SmagC2_2d*SmagC2_2d
               call update_2d_halo(SmagC2_2d,SmagC2_2d,az,imin,jmin,imax,jmax,H_TAG)
               call wait_halo(H_TAG)

               if (les_mode.eq.LES_MOMENTUM .or. les_mode.eq.LES_BOTH) then
                  do j=jmin-HALO,jmax+HALO-1
                     do i=imin-HALO,imax+HALO-1
                        SmagX2_2d(i,j) = _QUART_ * (  SmagC2_2d(i,j  ) + SmagC2_2d(i+1,j  ) &
                                                    + SmagC2_2d(i,j+1) + SmagC2_2d(i+1,j+1) )
                     end do
                  end do
               end if
               if (les_mode.eq.LES_TRACER .or. les_mode.eq.LES_BOTH) then
                  do j=jmin-HALO,jmax+HALO
                     do i=imin-HALO,imax+HALO-1
                        SmagU2_2d(i,j) = _HALF_ * ( SmagC2_2d(i,j) + SmagC2_2d(i+1,j  ) )
                     end do
                  end do
                  do j=jmin-HALO,jmax+HALO-1
                     do i=imin-HALO,imax+HALO
                        SmagV2_2d(i,j) = _HALF_ * ( SmagC2_2d(i,j) + SmagC2_2d(i  ,j+1) )
                     end do
                  end do
               end if

            case default
               stop 'init_les(): no valid smag_method specified'

         end select

      case default
         stop 'init_les(): no valid les_method specified'

   end select

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
! Copyright (C) 2011 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
