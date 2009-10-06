!$Id: parser.F90,v 1.1 2009-10-06 11:42:49 kb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 
!
! !INTERFACE:
   module parser
!
! !DESCRIPTION:
!  This module define a type $output_set$ defining variables - i.e.
!  hyperslab and variable list information - for allowing a very flexible
!  outut scheme in GETM. The actual full output specifications is in a
!  list of $output_set$'s. Information about output-definitions is read
!  from an ASCII-file with a relative simple format.
!  The file contains a number of lines with a given format. Lines starting
!  with \# or \! and empty lines are discarded. All non-discarded lines must
!  have a forma like:
!  \begin{verbatim}
!   all      |n=50:end:25|elev,temp,salt
!   test2d   |i=50:60, j=40:70:2, k=10:end,n=10:end:10|metrics,2d
!   test3d   |i=50:60, j=40:70:2, k=10:end,n=10:end:10|3d,U,V
!   all      |n=50:end:25|airp
!  \end{verbatim}
!  Each line has 3 fields separated by "|". The first field is the
!  id for an output set. The second field is the hyperslab - in indices
!  i,j,k,n and the last field is the variable list.\newline
!  Note that an id can be repeated. In this case the first hyperslab
!  definition is used and the variable list is the union. For the hyperslap
!  field it is not necessary to specify for all indices. If a given index 
!  is not defined the minimum and maximum values are used with a stride of
!  1. Note that the keywords $begin$ and $end$ can be used. The 3 field
!  contains a comma-separated list of variables to include iin the 
!  output-set. The variable names must match the names of variables given
!  in $variable_info.F90$
!  \newline
!  All implementation details not decided yet.
!  \newline
!  Work in progress and some issues still remaining.
!  A very first version of this routine was made by Alexander Barth.
!
! !USES:
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:

! !DEFINED PARAMETERS
   character, parameter :: field_sep = '|'
   character, parameter :: var_sep = ','
   character, parameter :: range_sep = ','
   character, parameter :: index_sep = ':'
   integer, parameter   :: max_length = 256
   integer, parameter   :: max_name   = 16
   character(len=max_name) :: metrics(8) = [ "dxc", "dyc", "dxu", "dyu", &
                                             "dxv", "dyv", "dxx", "dyx"  ]
   character(len=max_name) :: a2d(3) = [ "elev", "U", "V" ]
#ifdef NO_BAROCLINIC
   character(len=max_name) :: a3d(3) = [ "elev", "uu", "vv" ]
#else
   character(len=max_name) :: a3d(5) = [ "elev", "uu", "vv", "S", "T" ]
#endif
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  $Log: parser.F90,v $
!  Revision 1.1  2009-10-06 11:42:49  kb
!  first steps towards imore flexible output
!
!
!  subset of the model domain for output
   type output_set
      character(len=max_length) :: id
!      integer :: nfid

      character(len=max_name), allocatable :: var_list(:)
      integer, pointer :: varinfo_index(:)
      integer, pointer :: varids(:)

      integer :: il,ih,is
      integer :: jl,jh,js
      integer :: kl,kh,ks
      integer :: nl,nh,ns
   end type output_set

 ! list of all subsets
   type(output_set), pointer :: output_list(:)

!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: parse_output_list() - 
!
! !INTERFACE:
   subroutine parse_output_list(fn,imin,imax,jmin,jmax,kmin,kmax,nmin,nmax,ol)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  This subroutine opens and parses the file with output-set definitions.
!  On exit of the subroutine a list of output definitions are in the
!  array of $output_set$-type -  $ol$.
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)   :: fn
   integer, intent(in)            :: imin,imax
   integer, intent(in)            :: jmin,jmax
   integer, intent(in)            :: kmin,kmax
   integer, intent(in)            :: nmin,nmax
!
! !OUTPUT PARAMETERS:
   type(output_set), pointer :: ol(:)
!
! !LOCAL VARIABLES:
   character(len=max_length) :: line
   integer                   :: i,j,k=1,l,m,n
   integer                   :: n_output_set=-1
   integer                   :: nb_vars
   integer                   :: Nlines
   integer                   :: iostat
   integer                   :: iunit = 10
   character(len=max_length), allocatable :: id_list(:)
   character(len=max_length), allocatable :: lines(:)
   character(len=max_name),   allocatable :: var_names(:)
!EOP
!-------------------------------------------------------------------------
!BOC
   open(unit=iunit,file=fn,status='old')
   read(iunit,'(A)',iostat=iostat) line

! loop over all lines in file and get the number of 'valid' lines
   Nlines = 0
   do while (iostat == 0)
     if (line(1:1) == '#' .or. line(1:1) == '!' .or. &
         len_trim(line) .eq. 0) then
       ! comments (ignore)
      else
         Nlines=Nlines+1
      end if
      read(iunit,'(A)',iostat=iostat) line
   end do
   allocate(id_list(Nlines))
   allocate(lines(Nlines))
!write(*,*) 'Nlines= ',Nlines

!  Get full list of id's and read in lines with information
   rewind (iunit)
   read(iunit,'(A)',iostat=iostat) line
   n=0
   do while (iostat == 0)
      if (line(1:1) == '#' .or. line(1:1) == '!' .or. &
          len_trim(line) .eq. 0) then
        ! comments (ignore)
      else
         i = 1
         j = next_index(line,field_sep,i)-1
         n=n+1
         id_list(n) = trim(line(i:j))
         lines(n) = line
      end if
      read(iunit,'(A)',iostat=iostat) line
   end do
   close(iunit)

!  Get unique list of id's
   n_output_set = n
   i = n_output_set
   do n=1,n_output_set
      if (len_trim(id_list(n)) .gt. 0) then
         do k=n+1,n_output_set
            if (len_trim(id_list(k)) .gt. 0) then
               if (trim(id_list(n)) == trim(id_list(k))) then
                  id_list(k) = ''
                  i=i-1
               end if
            end if
         end do 
      end if
   end do
!write(*,*) 'Nuniq=  ',i
   n_output_set=i

!  allocate output list and copy unique id's
   allocate(ol(n_output_set))
   do n=1,n_output_set
      ol(n)%id = trim(id_list(n))
!write(*,*) trim(ol(n)%id)
   end do

!  get indices range
   do n=1,n_output_set
!     default: output everything
      ol(n)%il = imin; ol(n)%ih = imax; ol(n)%is = 1
      ol(n)%jl = jmin; ol(n)%jh = jmax; ol(n)%js = 1
      ol(n)%kl = kmin; ol(n)%kh = kmax; ol(n)%ks = 1
      ol(n)%nl = nmin; ol(n)%nh = nmax; ol(n)%ns = 1

!     i and j are start and end index of every token on a line
      k = 1
      i = 1
      j = next_index(lines(k),field_sep,i)-1
      do while (trim(ol(n)%id) .ne. trim(lines(k)(i:j)))
         k=k+1
      end do
!write(*,*) 'bb ',n,k,trim(lines(k)(i:j))
      j = next_index(lines(k),field_sep,i)-1
!      ol(k)%id = trim(line(i:j))

!     get indices range
      i=j+2
      j = next_index(lines(k),field_sep,i)-1

      m = next_index(lines(k),range_sep,i,j)-1
!write (*,*) 'aa ',i,j,m
!write (*,*) 'aa ',trim(lines(k))
!write (*,*) 'aa ',lines(k)(i:j)
!write (*,*) 'aa ',lines(k)(i:m)

      do while (i < j)
         l = next_index(lines(k),'=',i,m)-1

!write(*,*) 'cc ',n,k
!write(*,*) 'cc ',i,j,l,m
!write(*,*) 'cc ',lines(k)(i:l)
!write(*,*) 'cc ',lines(k)(i:m)
         if (lines(k)(l:l) == 'i') then
            i=l+2
            call parse_range(lines(k)(i:m),ol(n)%il,ol(n)%ih,ol(n)%is,imin,imax)   
         else if (lines(k)(l:l) == 'j') then
            i=l+2
            call parse_range(lines(k)(i:m),ol(n)%jl,ol(n)%jh,ol(n)%js,jmin,jmax)   
         else if (lines(k)(l:l) == 'k') then
            i=l+2
            call parse_range(lines(k)(i:m),ol(n)%kl,ol(n)%kh,ol(n)%ks,kmin,kmax)   
         else if (lines(k)(l:l) == 'n') then
            i=l+2
!write(*,*) 'dd ',lines(k)(i:m)
            call parse_range(lines(k)(i:m),ol(n)%nl,ol(n)%nh,ol(n)%ns,nmin,nmax)   
         end if

         i = m+2
         m = next_index(lines(k),range_sep,i,j)-1
      end do
   end do

!  concatenate variable lists from dublicated id's
   do n=1,n_output_set

!     i and j are start and end index of every token on a line
      do k=1,Nlines
         i = 1
         j = next_index(lines(k),field_sep,i)-1
         if (trim(ol(n)%id) == trim(lines(k)(i:j))) then
!           skip to variable list
            j = next_index(lines(k),field_sep,i)
            i=j+1
            j = next_index(lines(k),field_sep,i)
            i=j+1
!write(*,*) n,k,i,j
!write(*,*) lines(k)(i:)
            if (k .eq. n) then
               line = lines(k)(i:)
            else
               line = trim(line)//','//lines(k)(i:)
            end if
         end if
      end do
!write(*,*) 'vars= ',trim(line)

!     get rid of duplicate variable names
      nb_vars = count_chars(line,var_sep) + 1
      allocate(var_names(nb_vars))
      i = 1
      j = next_index(line,var_sep,i)-1
      do l=1,nb_vars
         var_names(l)=trim(line(i:j)) 
         i=j+2
         j = next_index(line,var_sep,i)-1
      end do
#if 0
do l=1,nb_vars
  write(*,*) 'aa ',l,trim(var_names(l))
end do
#endif
      m=nb_vars
      do l=1,nb_vars
         do k=l+1,nb_vars
            if (len_trim(var_names(l)) .gt. 0) then
               if(trim(var_names(l)) == trim(var_names(k))) then
                  var_names(k) = ''
                  m=m-1
               end if
            end if
         end do
      end do
      nb_vars=m
#if 0
i = 0 
do l=1,size(var_names)
   if (len_trim(var_names(l)) .gt. 0) then
      i = i+1
      write(*,*) 'bb ',l,trim(var_names(i))
   end if
end do
#endif
!     expand predefined definitions of variable lists
      m=nb_vars
      do l=1,nb_vars
         if (trim(var_names(l)) .eq. 'metrics') m=m+size(metrics)-1
         if (trim(var_names(l)) .eq. '2d')      m=m+size(a2d)-1
         if (trim(var_names(l)) .eq. '3d')      m=m+size(a3d)-1
         if (trim(var_names(l)) .eq. 'turb')    m=m+5-1
      end do
      nb_vars=m
!write(*,*) 'm= ',m
!write(*,*) 'n= ',n

!     allocate(ol(k)%varinfo_index(nb_vars))
      allocate(ol(n)%var_list(nb_vars))
      i = 0
      do l=1,size(var_names)
         if (len_trim(var_names(l)) .gt. 0) then
            if (trim(var_names(l)) == "metrics") then
               do j=1,size(metrics)
                  i = i+1
                  ol(n)%var_list(i) = metrics(j)
               end do
            else if (trim(var_names(l)) == "2d") then
               do j=1,size(a2d)
                  i = i+1
                  ol(n)%var_list(i) = a2d(j)
               end do
            else if (trim(var_names(l)) == "3d") then
               do j=1,size(a3d)
                  i = i+1
                  ol(n)%var_list(i) = a3d(j)
               end do
            else
               i = i+1
               ol(n)%var_list(i) = trim(var_names(l))
            end if
          end if
      end do
      deallocate(var_names)
   end do
   end subroutine parse_output_list
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: print_output_list - 
!
! !DESCRIPTION:
!
! !INTERFACE:
   subroutine print_output_list(ol)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type(output_set), pointer      :: ol(:)
!
! !LOCAL VARIABLES:
  integer          :: k,i
!EOP
!-------------------------------------------------------------------------
!BOC
  do k=1,size(ol)
    write(6,'(A,A)') 'id: ',trim(ol(k)%id)
    do i=1,size(ol(k)%var_list)
!      write(6,'(A,A)') '  variable: ', trim(varinfo_list(ol(k)%varinfo_index(i))%name)
      write(6,'(A,A)') '  variable: ', trim((ol(k)%var_list(i)))
    end do
    write(6,'(A,3I10)') '  i range: ',ol(k)%il,ol(k)%ih,ol(k)%is
    write(6,'(A,3I10)') '  j range: ',ol(k)%jl,ol(k)%jh,ol(k)%js
    write(6,'(A,3I10)') '  k range: ',ol(k)%kl,ol(k)%kh,ol(k)%ks
    write(6,'(A,3I10)') '  n range: ',ol(k)%nl,ol(k)%nh,ol(k)%ns
  end do
 end subroutine print_output_list
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: deallocate_output_list - 
!
! !DESCRIPTION:
!
! !INTERFACE:
   subroutine deallocate_output_list(ol)
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
  type(output_set), pointer       :: ol(:)
!
! !LOCAL VARIABLES:
   integer         :: k
!EOP
!-------------------------------------------------------------------------
!BOC
  do k=1,size(ol)
    deallocate(ol(k)%varinfo_index)
    deallocate(ol(k)%varids)
  end do
  deallocate(ol)
  end subroutine deallocate_output_list
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: next_index - 
!
! !INTERFACE:
   function next_index(s,sep,start,finish)
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)   :: s
   character, intent(in)          :: sep
   integer, intent(in)            :: start
   integer, intent(in), optional  :: finish
!
! !LOCAL VARIABLES:
  integer          :: next_index
   integer         :: i,f
!EOP
!-------------------------------------------------------------------------
!BOC
   if (present(finish)) then
      f = finish
   else
      f = len_trim(s)
   end if
   do i=start,f
      if (s(i:i) == sep) then
         next_index = i
         return
      end if
   end do
   next_index = f+1
   end function next_index
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: count_chars() - 
!
! !INTERFACE:
   function count_chars(s,sep)
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)   :: s
   character, intent(in)          :: sep
!
! !LOCAL VARIABLES:
   integer         :: i
   integer         :: count_chars
!EOP
!-------------------------------------------------------------------------
!BOC
   count_chars = 0
   do i=1,len_trim(s)
     if (s(i:i) == sep) then
       count_chars = count_chars+1
     end if
   end do
   end function count_chars
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: parse_index() - 
!
! !DESCRIPTION:
!
! !INTERFACE:
   function parse_index(s,min,max)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)   :: s
   integer, intent(in)            :: min,max
!
! !LOCAL VARIABLES:
   integer parse_index
!EOP
!-------------------------------------------------------------------------
!BOC
   if (s == 'end') then
     parse_index = max
   elseif (s == 'begin') then
     parse_index = min
   else
     read(s,*) parse_index
   end if
   end function parse_index
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: parse_range() - 
!
! !INTERFACE:
   subroutine parse_range(s,first,last,step,min,max)
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)   :: s
   integer, intent(in)            :: min,max
!
! !OUTPUT PARAMETERS:
   integer, intent(out)           :: first,last,step
!
! !LOCAL VARIABLES:
   integer         :: i,j  
!EOP
!-------------------------------------------------------------------------
!BOC
   i = 1
   j = next_index(s,index_sep,i)-1
   first = parse_index(s(i:j),min,max)
   i=j+2
   j = next_index(s,index_sep,i)-1
   if (len_trim(s(i:j)) == 0) then
      last = first
      step = 1
   else
      last = parse_index(s(i:j),min,max)

      i=j+2
      j = next_index(s,index_sep,i)-1

      if (len_trim(s(i:j)) == 0) then
         step = 1
      else
         step = parse_index(s(i:j),min,max)
      end if
   end if
   end subroutine parse_range
!EOC

!-----------------------------------------------------------------------

   end module parser

!-----------------------------------------------------------------------
! Copyright (C) 2009 - Karsten Bolding (BB)                            !
!-----------------------------------------------------------------------
