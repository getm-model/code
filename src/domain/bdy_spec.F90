!$Id: bdy_spec.F90,v 1.5 2005-04-29 12:55:31 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: bdy_spec() - defines the open boundaries.
!
! !INTERFACE:
   subroutine bdy_spec(FName)
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: NWB,NNB,NEB,NSB,NOB
   use domain, only: wi,wfj,wlj,nj,nfi,nli,ei,efj,elj,sj,sfi,sli
   use domain, only: bdy_index,bdy_map,nsbv
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: FName
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: bdy_spec.F90,v $
!  Revision 1.5  2005-04-29 12:55:31  kbk
!  removing print statement
!
!  Revision 1.4  2003/08/03 09:52:11  kbk
!  nicer print statements
!
!  Revision 1.3  2003/04/23 11:59:39  kbk
!  update_2d_halo on spherical variables + TABS to spaces
!
!  Revision 1.2  2003/04/07 15:20:53  kbk
!  added bdy_index and bdy_map
!
!  Revision 1.1.1.1  2002/05/02 14:01:11  gotm
!  recovering after CVS crash
!
!  Revision 1.2  2001/09/01 17:12:13  bbh
!  Removed a STDERR
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,l,n,rc
!
!EOP
!-----------------------------------------------------------------------
!BOC

#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'bdy_spec() # ',Ncall
#endif

   LEVEL2 'Reading boundary location information from:'
   LEVEL3 TRIM(FName)
   open(BDYINFO,file=FName,status='unknown',ERR=90)

!  Western boundary info
   read(BDYINFO,*,END=91,ERR=92) NWB
   if (NWB .ge. 1) then
      allocate(wi(NWB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (wi)'
      allocate(wfj(NWB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (wfj)'
      allocate(wlj(NWB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (wlj)'
      do n = 1,NWB
         read(BDYINFO,*,END=91,ERR=92) wi(n),wfj(n),wlj(n)
         nsbv = nsbv + (wlj(n)-wfj(n)+1)
      end do
   end if

!  Northern boundary info
   read(BDYINFO,*,END=91,ERR=92) NNB
   if (NNB .ge. 1 ) then
      allocate(nj(NNB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (nj)'
      allocate(nfi(NNB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (nfi)'
      allocate(nli(NNB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (nli)'
      do n = 1,NNB
         read(BDYINFO,*,END=91,ERR=92) nj(n),nfi(n),nli(n)
         nsbv = nsbv + (nli(n)-nfi(n)+1)
      end do
   end if

!  Eastern boundary info
   read(BDYINFO,*,END=91,ERR=92) NEB
   if (NEB .ge. 1 ) then
      allocate(ei(NEB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (ei)'
      allocate(efj(NEB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (efj)'
      allocate(elj(NEB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (elj)'
      do n = 1,NEB
         read(BDYINFO,*,END=91,ERR=92) ei(n),efj(n),elj(n)
         nsbv = nsbv + (elj(n)-efj(n)+1)
      end do
   end if

!  Southern boundary info
   read(BDYINFO,*,END=91,ERR=92) NSB
   if (NSB .ge. 1 ) then
      allocate(sj(NSB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (sj)'
      allocate(sfi(NSB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (sfi)'
      allocate(sli(NSB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (sli)'
      do n = 1,NSB
         read(BDYINFO,*,END=91,ERR=92) sj(n),sfi(n),sli(n)
         nsbv = nsbv + (sli(n)-sfi(n)+1)
      end do
   end if

   close(BDYINFO)

   NOB = NWB+NNB+NEB+NSB

   LEVEL2 '# of open boundaries ',NOB
   LEVEL2 '# of open bdy points ',nsbv

   if(NOB .gt. 0) then
      allocate(bdy_index(NOB),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (bdy_index)'

      allocate(bdy_map(nsbv,2),stat=rc)
      if (rc /= 0) stop 'bdy_spec: Error allocating memory (bdy_map)'

      k = 1
      l = 1
      do n=1,NWB
         bdy_index(l) = k
         l = l+1
         i = wi(n)
         do j=wfj(n),wlj(n)
            bdy_map(k,1) = i
            bdy_map(k,2) = j
            k = k+1
         end do
      end do
      do n=1,NNB
         bdy_index(l) = k
         l = l+1
         j = nj(n)
         do i=nfi(n),nli(n)
            bdy_map(k,1) = i
            bdy_map(k,2) = j
            k = k+1
         end do
      end do
      do n=1,NEB
         bdy_index(l) = k
         l = l+1
         i = ei(n)
         do j=efj(n),elj(n)
            bdy_map(k,1) = i
            bdy_map(k,2) = j
            k = k+1
         end do
      end do
      do n=1,NSB
         bdy_index(l) = k
         l = l+1
         j = sj(n)
         do i=sfi(n),sli(n)
            bdy_map(k,1) = i
            bdy_map(k,2) = j
            k = k+1
         end do
      end do
   end if

#ifdef DEBUG
   write(debug,*) 'Leaving bdy_spec()'
   write(debug,*)
#endif

   return
90 FATAL 'can not open ',TRIM(FName)
   stop
91 STDERR 'EOF ',TRIM(FName)
92 STDERR 'Error reading ',TRIM(FName)

   return
   end subroutine bdy_spec
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
