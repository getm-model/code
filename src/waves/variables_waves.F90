#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: variables_waves - global variables for waves
!
! !INTERFACE:
   module variables_waves
!
! !DESCRIPTION:
!  The module contains public subroutines to initialise and cleanup.
!
! !USES:
   use domain, only: imin,imax,jmin,jmax,kmax

   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   REALTYPE,dimension(:,:)  ,allocatable :: waveDir,coswavedir,sinwavedir
   REALTYPE,dimension(:,:)  ,allocatable :: waveH,waveL,waveT,waveK,waveE
   REALTYPE,dimension(:,:)  ,allocatable :: SJ,SJJ
   REALTYPE,dimension(:,:)  ,allocatable :: kDveln,sinh2kDvelnm1
   logical ,dimension(:,:)  ,allocatable :: is_deepwave
   REALTYPE,dimension(:,:,:),allocatable :: khab,layerratios
   REALTYPE,dimension(:,:)  ,allocatable :: UStokesC,VStokesC
   REALTYPE,dimension(:,:)  ,allocatable :: UStokes,VStokes
   REALTYPE,dimension(:,:)  ,allocatable :: UStokesCint,VStokesCint
   REALTYPE,dimension(:,:)  ,allocatable :: UStokesCadv,VStokesCadv
   REALTYPE,dimension(:,:,:),allocatable :: uuStokesC,vvStokesC
   REALTYPE,dimension(:,:,:),allocatable :: uuStokes,vvStokes
!
! !REVISION HISTORY:
!  Original author(s): Ulf Graewe
!                      Saeed Moghimi
!                      Knut Klingbeil
!
! !LOCAL VARIABLES:
   integer :: rc
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_variables_waves - initialise waves related stuff.
!
! !INTERFACE:
   subroutine init_variables_waves(runtype)
!
! !DESCRIPTION:
!  Allocates memory (unless {\tt STATIC} is set) for les related fields,
!  by an include statement.

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in) :: runtype
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_variables_waves() # ',Ncall
#endif

   LEVEL2 'init_variables_waves'

   allocate(waveDir(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_variables_waves: Error allocating memory (waveDir)'

   allocate(coswavedir(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_variables_waves: Error allocating memory (coswavedir)'

   allocate(sinwavedir(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_variables_waves: Error allocating memory (sinwavedir)'

   allocate(waveH(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_variables_waves: Error allocating memory (waveH)'

   allocate(waveL(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_variables_waves: Error allocating memory (waveL)'

   allocate(waveT(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_variables_waves: Error allocating memory (waveT)'

   allocate(waveK(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_variables_waves: Error allocating memory (waveK)'

   allocate(waveE(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_variables_waves: Error allocating memory (waveE)'

   allocate(SJ(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_variables_waves: Error allocating memory (SJ)'

   allocate(UStokesC(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_variables_waves: Error allocating memory (UStokesC)'
   UStokesC = _ZERO_

   allocate(VStokesC(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_variables_waves: Error allocating memory (VStokesC)'
   VStokesC = _ZERO_

   allocate(UStokes(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_variables_waves: Error allocating memory (UStokes)'
   UStokes = _ZERO_

   allocate(VStokes(E2DFIELD),stat=rc)
   if (rc /= 0) stop 'init_variables_waves: Error allocating memory (VStokes)'
   VStokes = _ZERO_


   if (runtype .gt. 1) then

      allocate(SJJ(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_variables_waves: Error allocating memory (SJJ)'

      allocate(kDveln(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_variables_waves: Error allocating memory (kDveln)'

      allocate(sinh2kDvelnm1(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_variables_waves: Error allocating memory (sinh2kDvelnm1)'

      allocate(is_deepwave(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_variables_waves: Error allocating memory (is_deepwave)'

      allocate(khab(I3DFIELD),stat=rc)
      if (rc /= 0) stop 'init_variables_waves: Error allocating memory (khab)'

      allocate(layerratios(I3DFIELD),stat=rc)
      if (rc /= 0) stop 'init_variables_waves: Error allocating memory (layerratios)'

      allocate(UStokesCint(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_variables_waves: Error allocating memory (UStokesCint)'
      UStokesCint = _ZERO_

      allocate(VStokesCint(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_variables_waves: Error allocating memory (VStokesCint)'
      VStokesCint = _ZERO_

      allocate(UStokesCadv(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_variables_waves: Error allocating memory (UStokesCadv)'
      UStokesCadv = _ZERO_

      allocate(VStokesCadv(I2DFIELD),stat=rc)
      if (rc /= 0) stop 'init_variables_waves: Error allocating memory (VStokesCadv)'
      VStokesCadv = _ZERO_

      allocate(uuStokesC(I3DFIELD),stat=rc)
      if (rc /= 0) stop 'init_variables_waves: Error allocating memory (uuStokesC)'
      uuStokesC = _ZERO_

      allocate(vvStokesC(I3DFIELD),stat=rc)
      if (rc /= 0) stop 'init_variables_waves: Error allocating memory (vvStokesC)'
      vvStokesC = _ZERO_

      allocate(uuStokes(I3DFIELD),stat=rc)
      if (rc /= 0) stop 'init_variables_waves: Error allocating memory (uuStokes)'
      uuStokes = _ZERO_

      allocate(vvStokes(I3DFIELD),stat=rc)
      if (rc /= 0) stop 'init_variables_waves: Error allocating memory (vvStokes)'
      vvStokes = _ZERO_

   end if


#ifdef DEBUG
   write(debug,*) 'Leaving init_variables_waves()'
   write(debug,*)
#endif
   return
   end subroutine init_variables_waves
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clean_variables_waves - cleanup after run.
!
! !INTERFACE:
   subroutine clean_variables_waves()
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This routine is currently empty.
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'clean_variables_waves() # ',Ncall
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving clean_variables_waves()'
   write(debug,*)
#endif
   return
   end subroutine clean_variables_waves
!EOC
!-----------------------------------------------------------------------

   end module variables_waves

!-----------------------------------------------------------------------
! Copyright (C) 2013 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
