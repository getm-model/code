!$Id: nesting.F90,v 1.3 2003-04-23 12:07:12 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: nesting - interface for nesting routines
!
! !INTERFACE:
   module nesting
!
! !DESCRIPTION:
!  Still missing
!
! !USES:
   use domain, only: imin,imax,jmin,jmax
   use domain, only: iimin,iimax,jjmin,jjmax,kmax
   use time, only: julianday,secondsofday
   use domain, only: H,HU,HV
#ifndef NO_3D
   use variables_3d, only: hun,uu,hvn,vv
#ifndef NO_BAROCLINIC
   use variables_3d, only: S,T
#endif
#endif
#ifndef NO_BAROCLINIC
#endif
   IMPLICIT NONE
!
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_nesting, nesting_file, clean_nesting
!
! !PUBLIC DATA MEMBERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: nesting.F90,v $
!  Revision 1.3  2003-04-23 12:07:12  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 12:32:58  kbk
!  parallel support + NO_3D, NO_BAROCLINIC
!
!  Revision 1.1.1.1  2002/05/02 14:01:53  gotm
!  recovering after CVS crash
!
!
! !PRIVATE DATA MEMBERS
   integer                   :: NW_w,NN_w,NE_w,NS_w,nonp_w=0
   REALTYPE, allocatable     :: wrk_2d_w(:),wrk_3d_w(:,:)
   integer, dimension(:), allocatable :: wi_w,wfj_w,wlj_w
   integer, dimension(:), allocatable :: nj_w,nfi_w,nli_w
   integer, dimension(:), allocatable :: ei_w,efj_w,elj_w
   integer, dimension(:), allocatable :: sj_w,sfi_w,sli_w

   integer                   :: NW_r,NN_r,NE_r,NS_r,nonp_r=0
   REALTYPE, allocatable     :: wrk_2d_r(:),wrk_3d_r(:,:)
   integer, dimension(:), allocatable :: wi_r,wfj_r,wlj_r
   integer, dimension(:), allocatable :: nj_r,nfi_r,nli_r
   integer, dimension(:), allocatable :: ei_r,efj_r,elj_r
   integer, dimension(:), allocatable :: sj_r,sfi_r,sli_r

   integer, parameter        :: TGRID=1,UGRID=2,VGRID=3
   integer, parameter        :: uout=60,uin=61 ! kbk
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_nesting
!
! !INTERFACE:
   subroutine init_nesting()
!
! !DESCRIPTION:
!  Description still missing
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer                   :: iunit=25 ! kbk
   integer                   :: dum
   integer                   :: rc
   integer                   :: k,n
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_nesting() # ',Ncall
#endif
   LEVEL2 'init_nesting()'

   LEVEL2 'Reading nesting information from ','nesting.info'
   open(unit=iunit,file='nesting.info',status='unknown',err=90)

   read(iunit,*,end=91,err=92) dum

   do k=1,1

!  Western
   read(iunit,*,end=91,err=92) NW_w
   if (NW_w .ge. 1) then   
      allocate(wi_w(NW_w),stat=rc)
      if (rc /= 0) stop 'init_nesting: Error allocating memory (wi_w)'
      allocate(wfj_w(NW_w),stat=rc)
      if (rc /= 0) stop 'init_nesting: Error allocating memory (wfj_w)'
      allocate(wlj_w(NW_w),stat=rc)
      if (rc /= 0) stop 'init_nesting: Error allocating memory (wlj_w)'
      do n = 1,NW_w
         read(iunit,*,end=91,err=92) wi_w(n),wfj_w(n),wlj_w(n)
         nonp_w = nonp_w + (wlj_w(n)-wfj_w(n)+1)
      end do
   end if

!  Northen
   read(iunit,*,end=91,err=92) NN_w
   if (NN_w .ge. 1) then   
      allocate(nj_w(NN_w),stat=rc)
      if (rc /= 0) stop 'init_nesting: Error allocating memory (nj_w)'
      allocate(nfi_w(NN_w),stat=rc)
      if (rc /= 0) stop 'init_nesting: Error allocating memory (nfi_w)'
      allocate(nli_w(NN_w),stat=rc)
      if (rc /= 0) stop 'init_nesting: Error allocating memory (nli_w)'
      do n = 1,NN_w
         read(iunit,*,end=91,err=92) nj_w(n),nfi_w(n),nli_w(n)
         nonp_w = nonp_w + (nli_w(n)-nfi_w(n)+1)
      end do
   end if

!  Eastern
   read(iunit,*,end=91,err=92) NE_w
   if (NE_w .ge. 1) then
      allocate(ei_w(NE_w),stat=rc)
      if (rc /= 0) stop 'init_nesting: Error allocating memory (ei_w)'
      allocate(efj_w(NE_w),stat=rc)
      if (rc /= 0) stop 'init_nesting: Error allocating memory (efj_w)'
      allocate(elj_w(NE_w),stat=rc)
      if (rc /= 0) stop 'init_nesting: Error allocating memory (elj_w)'
      do n = 1,NE_w
         read(iunit,*,end=91,err=92) ei_w(n),efj_w(n),elj_w(n)
         nonp_w = nonp_w + (elj_w(n)-efj_w(n)+1)
      end do
   end if

!  Southern
   read(iunit,*,end=91,err=92) NS_w
   if (NS_w .ge. 1) then
      allocate(sj_w(NS_w),stat=rc)
      if (rc /= 0) stop 'init_nesting: Error allocating memory (sj_w)'
      allocate(sfi_w(NS_w),stat=rc)
      if (rc /= 0) stop 'init_nesting: Error allocating memory (sfi_w)'
      allocate(sli_w(NS_w),stat=rc)
      if (rc /= 0) stop 'init_nesting: Error allocating memory (sli_w)'
      do n = 1,NS_w
         read(iunit,*,end=91,err=92) sj_w(n),sfi_w(n),sli_w(n)
         nonp_w = nonp_w + (sli_w(n)-sfi_w(n)+1)
      end do
   end if

   end do

   close(iunit)

   allocate(wrk_2d_w(nonp_w),stat=rc)
   if (rc /= 0) stop 'init_nesting: Error allocating memory (wrk_2d_w)'
   allocate(wrk_3d_w(nonp_w,0:kmax),stat=rc)
   if (rc /= 0) stop 'init_nesting: Error allocating memory (wrk_3d_w)'

   return
90 FATAL 'can not open ','nesting.info'
   stop
91 STDERR 'EOF ','nesting.info'
   stop
92 STDERR 'Error reading ','nesting.info'
   stop
   end subroutine init_nesting
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: nesting_file
!
! !INTERFACE:
   subroutine nesting_file(mode)
!
! !DESCRIPTION:
!  Description still missing
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: mode
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,n
   integer                   :: rc
   logical                   :: first_w=.true.,first_r=.true.
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'nesting_file() # ',Ncall
#endif
   if (mode .eq. WRITING) then
      if(first_w) then
         first_w=.false.
         open(unit=uout,file='nesting.out',  &
              status='unknown',form='unformatted')
         write(uout) NW_w,NN_w,NE_w,NS_w,nonp_w
         if(NW_w .gt. 0) then
            write(uout) wi_w,wfj_w,wlj_w
         end if
         if(NN_w .gt. 0) then
            write(uout) nj_w,nfi_w,nli_w
         end if
         if(NE_w .gt. 0) then
            write(uout) ei_w,efj_w,elj_w
         end if
         if(NS_w .gt. 0) then
            write(uout) sj_w,sfi_w,sli_w
         end if
         call collect_2d(TGRID,H)
         write(uout) wrk_2d_w
         call collect_2d(UGRID,HU)
         write(uout) wrk_2d_w
         call collect_2d(VGRID,HV)
         write(uout) wrk_2d_w
      end if
#if 0
      write(uout) julianday,secondsofday
      call collect_3d(UGRID,hun)
      write(uout) wrk_3d_w
      call collect_3d(UGRID,uu)
      write(uout) wrk_3d_w
      call collect_3d(VGRID,hvn)
      write(uout) wrk_3d_w
      call collect_3d(VGRID,vv)
      write(uout) wrk_3d_w
      call collect_3d(TGRID,T)
      write(uout) wrk_3d_w
      call collect_3d(TGRID,S)
      write(uout) wrk_3d_w
#endif
   end if

   if (mode .eq. READING) then
      if(first_r) then
         first_r=.false.
         open(unit=uin,file='nesting.in',  &
              status='unknown',form='unformatted')
         read(uin) NW_r,NN_r,NE_r,NS_r,nonp_r
         if(NW_r .gt. 0) then
            read(uin) wi_r,wfj_r,wlj_r
         end if
         if(NN_r .gt. 0) then
            read(uin) nj_r,nfi_r,nli_r
         end if
         if(NE_r .gt. 0) then
            read(uin) ei_r,efj_r,elj_r
         end if
         if(NS_r .gt. 0) then
            read(uin) sj_r,sfi_r,sli_r
         end if
         read(uin) wrk_2d_r
      end if
      allocate(wrk_2d_r(nonp_r),stat=rc)
      if (rc /= 0) stop 'init_nesting: Error allocating memory (wrk_2d_r)'
      allocate(wrk_3d_w(nonp_r,0:kmax),stat=rc)
      if (rc /= 0) stop 'init_nesting: Error allocating memory (wrk_3d_r)'
      read(uin) julianday,secondsofday
      read(uin) wrk_3d_r
      read(uin) wrk_3d_r
      read(uin) wrk_3d_r
      read(uin) wrk_3d_r
      read(uin) wrk_3d_r
      read(uin) wrk_3d_r
   end if

   return
   end subroutine nesting_file
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clean_nesting
!
! !INTERFACE:
   subroutine clean_nesting()
!
! !DESCRIPTION:
!  Description still missing
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'clean_nesting() # ',Ncall
#endif
   LEVEL2 'clean_nesting()'
!close(uin)
!close(uout)
   return
   end subroutine clean_nesting
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: collect_2d
!
! !INTERFACE:
   subroutine collect_2d(source_grid,f)
!
! !DESCRIPTION:
!  Description still missing
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer                             :: source_grid
   REALTYPE, intent(in)                :: f(E2DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,n
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'collect_2d() # ',Ncall
#endif
   k = 0
   do n=1,NW_w
      i = wi_w(n)
      do j=wfj_w(n),wlj_w(n)
         k = k+1
         select case (source_grid)
            case (TGRID)
               wrk_2d_w(k) = 0.5*(f(i,j)+f(i+1,j))
            case (UGRID)
               wrk_2d_w(k) = f(i,j)
            case (VGRID)
               wrk_2d_w(k) = 0.25*(f(i,j)+f(i,j-1)+f(i+1,j-1)+f(i+1,j))
         end select
      end do
   end do
   do n=1,NN_w
      j = nj_w(n)
      do i=nfi_w(n),nli_w(n)
         k = k+1
         select case (source_grid)
            case (TGRID)
               wrk_2d_w(k) = 0.5*(f(i,j)+f(i,j+1))
            case (UGRID)
               wrk_2d_w(k) = 0.25*(f(i,j)+f(i,j+1)+f(i-1,j)+f(i-1,j+1))
            case (VGRID)
               wrk_2d_w(k) = f(i,j)
         end select
      end do
   end do
   do n=1,NE_w
      i = ei_w(n)
      do j=efj_w(n),elj_w(n)
         k = k+1
         select case (source_grid)
            case (TGRID)
               wrk_2d_w(k) = 0.5*(f(i,j)+f(i+1,j))
            case (UGRID)
               wrk_2d_w(k) = f(i,j)
            case (VGRID)
               wrk_2d_w(k) = 0.25*(f(i,j)+f(i,j-1)+f(i+1,j-1)+f(i+1,j))
         end select
      end do
   end do
   do n=1,NS_w
      j = sj_w(n)
      do i=sfi_w(n),sli_w(n)
         k = k+1
         select case (source_grid)
            case (TGRID)
               wrk_2d_w(k) = 0.5*(f(i,j)+f(i,j+1))
            case (UGRID)
               wrk_2d_w(k) = 0.25*(f(i,j)+f(i,j+1)+f(i-1,j)+f(i-1,j+1))
            case (VGRID)
               wrk_2d_w(k) = f(i,j)
         end select
      end do
   end do
   return
   end subroutine collect_2d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: collect_3d
!
! !INTERFACE:
   subroutine collect_3d(source_grid,f)
!
! !DESCRIPTION:
!  Description still missing
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer                             :: source_grid
   REALTYPE, intent(in)                :: f(I3DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,n
   REALTYPE                  :: col(0:kmax)
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'collect_3d() # ',Ncall
#endif
   k = 0
   do n=1,NW_w
      i = wi_w(n)
      do j=wfj_w(n),wlj_w(n)
         k = k+1
         select case (source_grid)
            case (TGRID)
               wrk_3d_w(k,:) = 0.5*(f(i,j,:)+f(i+1,j,:))
            case (UGRID)
               wrk_3d_w(k,:) = f(i,j,:)
            case (VGRID)
               wrk_3d_w(k,:) = &
                  0.25*(f(i,j,:)+f(i,j-1,:)+f(i+1,j-1,:)+f(i+1,j,:))
         end select
      end do
   end do
   do n=1,NN_w
      j = nj_w(n)
      do i=nfi_w(n),nli_w(n)
         k = k+1
         select case (source_grid)
            case (TGRID)
               wrk_3d_w(k,:) = 0.5*(f(i,j,:)+f(i,j+1,:))
            case (UGRID)
               wrk_3d_w(k,:) = &
                  0.25*(f(i,j,:)+f(i,j+1,:)+f(i-1,j,:)+f(i-1,j+1,:))
            case (VGRID)
               wrk_3d_w(k,:) = f(i,j,:)
         end select
      end do
   end do
   do n=1,NE_w
      i = ei_w(n)
      do j=efj_w(n),elj_w(n)
         k = k+1
         select case (source_grid)
            case (TGRID)
               wrk_3d_w(k,:) = 0.5*(f(i,j,:)+f(i+1,j,:))
            case (UGRID)
               wrk_3d_w(k,:) = f(i,j,:)
            case (VGRID)
               wrk_3d_w(k,:) = &
                  0.25*(f(i,j,:)+f(i,j-1,:)+f(i+1,j-1,:)+f(i+1,j,:))
         end select
      end do
   end do
   do n=1,NS_w
      j = sj_w(n)
      do i=sfi_w(n),sli_w(n)
         k = k+1
         select case (source_grid)
            case (TGRID)
               wrk_3d_w(k,:) = 0.5*(f(i,j,:)+f(i,j+1,:))
            case (UGRID)
               wrk_3d_w(k,:) = &
                  0.25*(f(i,j,:)+f(i,j+1,:)+f(i-1,j,:)+f(i-1,j+1,:))
            case (VGRID)
               wrk_3d_w(k,:) = f(i,j,:)
         end select
      end do
   end do

   return
   end subroutine collect_3d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: spread_3d
!
! !INTERFACE:
   subroutine spread_3d(source_grid,f)
!
! !DESCRIPTION:
!  Description still missing
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer                             :: source_grid
   REALTYPE, intent(in)                :: f(I3DFIELD)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,n
   REALTYPE                  :: col(0:kmax)
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'spread_3d() # ',Ncall
#endif
   LEVEL2 'spread_3d()'
#if 0
   k = 0
   do n=1,NW_r
      i = wi_r(n)
      do j=wfj_r(n),wlj_r(n)
         k = k+1
         select case (source_grid)
            case (TGRID)
               wrk_3d_r(k,:) = 0.5*(f(i,j,:)+f(i+1,j,k))
            case (UGRID)
               wrk_3d_r(k,:) = f(i,j,:)
            case (VGRID)
               wrk_3d_r(k,:) = &
                  0.25*(f(i,j,:)+f(i,j-1,:)+f(i+1,j-1,:)+f(i+1,j,:))
         end select
      end do
   end do
   do n=1,NN_r
      j = nj_r(n)
      do i=nfi_r(n),nli_r(n)
         k = k+1
         select case (source_grid)
            case (TGRID)
               wrk_3d_r(k,:) = 0.5*(f(i,j,:)+f(i,j+1,k))
            case (UGRID)
               wrk_3d_r(k,:) = &
                  0.25*(f(i,j,:)+f(i,j+1,:)+f(i-1,j,:)+f(i-1,j+1,:))
            case (VGRID)
               wrk_3d_r(k,:) = f(i,j,:)
         end select
      end do
   end do
   do n=1,NE_r
      i = ei_r(n)
      do j=efj_r(n),elj_r(n)
         k = k+1
         select case (source_grid)
            case (TGRID)
               wrk_3d_r(k,:) = 0.5*(f(i,j,:)+f(i+1,j,k))
            case (UGRID)
               wrk_3d_r(k,:) = f(i,j,:)
            case (VGRID)
               wrk_3d_r(k,:) = &
                  0.25*(f(i,j,:)+f(i,j-1,:)+f(i+1,j-1,:)+f(i+1,j,:))
         end select
      end do
   end do
   do n=1,NS_r
      j = sj_r(n)
      do i=sfi_r(n),sli_r(n)
         k = k+1
         select case (source_grid)
            case (TGRID)
               wrk_3d_r(k,:) = 0.5*(f(i,j,:)+f(i,j+1,k))
            case (UGRID)
               wrk_3d_r(k,:) = &
                  0.25*(f(i,j,:)+f(i,j+1,:)+f(i-1,j,:)+f(i-1,j+1,:))
            case (VGRID)
               wrk_3d_r(k,:) = f(i,j,:)
         end select
      end do
   end do
#endif
   return
   end subroutine spread_3d
!EOC

!-----------------------------------------------------------------------

   end module nesting

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Karsten Bolding and Hans Burchard (BBH)         !
!-----------------------------------------------------------------------
