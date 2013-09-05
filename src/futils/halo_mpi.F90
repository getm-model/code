#include "cppdefs.h"
#ifndef HALO
#define HALO 0
#endif
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: halo_mpi - mpi interface to getm
!
! !INTERFACE:
   module halo_mpi

! !DESCRIPTION:
!  This module provides necessary routines to enable parallel execution
!  of 'getm' using the Message Passage Interface (MPI) method.
!
! !USES:
   use mpi
   IMPLICIT NONE
!
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public                              :: init_mpi,postinit_mpi
   public                              :: print_MPI_info,barrier
   public                              :: set_active_communicator
   public                              :: update_2d_halo_mpi
   public                              :: update_3d_halo_mpi
   public                              :: wait_halo_mpi
   public                              :: set_flag_mpi
   public:: part_domain_mpi
   integer, public, parameter          :: H_TAG=10,HU_TAG=11,HV_TAG=12
   integer, public, parameter          :: D_TAG=20,DU_TAG=21,DV_TAG=22
   integer, public, parameter          :: z_TAG=30,U_TAG=31,V_TAG=32

   private                             :: MPI_data_types

! !PUBLIC DATA MEMBERS:
   integer, public                     :: myid=-1, nprocs=1
   integer, public                     :: comm_hd=MPI_COMM_WORLD
   LONGINT, public                     :: all_2d_exchange, all_3d_exchange
!   integer, public                    :: comm_wave=MPI_COMM_WORLD
!   integer, public                    :: comm_biology=MPI_COMM_WORLD
!
!  !DEFINED PARAMETERS:
   integer, parameter                  :: nneighbours=8
!  Different mesh specification methods
   integer, parameter                  :: ONE_CELL=-1
   integer, parameter                  :: ONED_MESH=0
   integer, parameter                  :: TWOD_MESH=1
   integer, parameter                  :: MESH_FROM_FILE=2
!  Methods of communication
   integer, parameter                  :: ONE_PROCESS=-1
   integer, parameter                  :: ONED_SENDRECV=0
   integer, parameter                  :: ONED_NONBLOCKING=1
   integer, parameter                  :: TWOD_SENDRECV=2
   integer, parameter                  :: TWOD_NONBLOCKING=3
!  Direction in case of ONED_? communications
   integer, parameter                  :: RIGHT_LEFT=1
   integer, parameter                  :: DOWN_UP=2
!  Last action
   integer, parameter                  :: SENDING=1
   integer, parameter                  :: WAITING=2
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   character(LEN = 256), private:: pname
   integer                   :: active_comm=MPI_COMM_WORLD
   integer                   :: mesh_method=ONE_CELL
   logical                   :: re_order=.false.
   integer                   :: comm_method=ONE_PROCESS
   integer                   :: len
   integer                   :: ierr
   integer                   :: x_line,x_lines,y_line,y_lines
   integer                   :: halo_line,halo_square
!KBK   integer                   :: i1_slice
   integer                   :: xz_slice,xz_slices
   integer                   :: yz_slice,yz_slices
   integer                   :: z_column
   integer                   :: halo_columns
!
! There is a discrepancy between new and older versions of MPI(2).
! In particular, the use of MPI_TYPE_EXTENT is deprecated, but _MAY_ still
! not be available in all mpi-lib implementations. 
! Also, mpif.h _should_  contain a definition of MPI_ADDRESS_KIND, but that 
! may also not always be the case. If x_size etc. is defined wrong, then the 
! data exchanged between subdomains may be wrong.
! The "old" behaviour (without MPI_ADDRESS_KIND and MPI_TYPE_GET_EXTENT) can
! be tried by defining _MPI_TYPE_EXTENT_ as a directive at compile time, e.g. by 
! adding -D_MPI_TYPE_EXTENT_ to the compilation.
#ifdef _MPI_TYPE_EXTENT_
   integer                   :: x_size,y_size,z_size
   integer                   :: xy_size,xz_size,yz_size,xyz_size
#else
   integer(KIND=MPI_ADDRESS_KIND) :: x_size,y_size,z_size
   integer(KIND=MPI_ADDRESS_KIND) :: xy_size,xz_size,yz_size,xyz_size
#endif
   integer                   :: com_direction
   integer                   :: req(2*nneighbours)
   integer                   :: status_array(MPI_STATUS_SIZE,2*nneighbours)
   integer                   :: dims(2),coords(2)
   logical                   :: periods(2)
   integer                   :: up,down,left,right
   integer                   :: ur,ul,ll,lr
   integer                   :: status(MPI_STATUS_SIZE)
   integer                   :: last_action=WAITING
   integer                   :: size_point
   integer                   :: size_left, size_ul,size_up,  size_ur
   integer                   :: size_right,size_lr,size_down,size_ll
!EOP
!-----------------------------------------------------------------------
!BOC

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_mpi - initialize the basic MPI environment
!
! !INTERFACE:
   subroutine init_mpi
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Initialize MPI parallel environment, i.e. getting process id etc.
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!  Revised by: Bjarne Buchmann, 2006
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   write(debug,*) 'init_mpi'
#endif

!  Initialize the MPI environment
   call MPI_INIT(ierr)
   if(ierr .ne. MPI_SUCCESS) then
      STDERR 'Fatal error: unable to initialize MPI.'
      call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
   end if

!  Get number of processes
   call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
   if (ierr .ne. MPI_SUCCESS) THEN
      STDERR 'Fatal error: unable to get number of processes.'
      call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
   end if

!  Get rank of current process
   call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
   if (ierr .ne. MPI_SUCCESS) THEN
      STDERR 'Fatal error: unable to get MYID.'
      call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
   end if

!  Get the processor names
   call MPI_GET_PROCESSOR_NAME(pname,len,ierr)
   if(ierr .ne. MPI_SUCCESS) THEN
      STDERR 'Fatal error: unable to get processor name.'
      call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
   end if

   all_2d_exchange = 0
   all_3d_exchange = 0


#ifdef DEBUG
   write(debug,*) 'Leaving init_mpi()'
   write(debug,*)
#endif
   return
   end subroutine init_mpi
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: postinit_mpi - more MPI initialization
!
! !INTERFACE:
   subroutine postinit_mpi(input_dir)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Initialization that requires read namelist.
!
! !INPUT PARAMTERS:
   character(len=*)                    :: input_dir
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!  Revised by: Bjarne Buchmann, 2006
!
! !LOCAL VARIABLES:
   integer                   :: MeshMethod,MsgMethod
   logical                   :: reorder
   namelist /nampar/ &
             MeshMethod,reorder,MsgMethod
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   write(debug,*) 'postinit_mpi'
#endif

!  Read parallel/MPI specific things from the namelist.
   open(11,file=trim(input_dir) // 'parallel.inp')
   read(11,nampar)
   close(11)

   call set_com_method(MsgMethod)
   call set_mesh_method(MeshMethod,reorder)

#ifdef DEBUG
   write(debug,*) 'Leaving postinit_mpi()'
   write(debug,*)
#endif
   return
   end subroutine postinit_mpi
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: barrier - stop until all processes reach this point
!
! !INTERFACE:
   subroutine barrier()
   IMPLICIT NONE
!
! !DESCRIPTION:
!  When this subroutine is called all processes wait for all others to
!  reach this point in the execution. Use with care - slows down.
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-------------------------------------------------------------------------
!BOC
   STDERR 'barrier(): we are waiting'
   call MPI_BARRIER(active_comm,ierr)
   return
   end subroutine barrier
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: print_MPI_info - write various MPI related info.
!
! !INTERFACE:
   SUBROUTINE print_MPI_info()
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Used to print information concerning the actual process. Id, name of
!  processor etc..
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer         :: ver,subver
   character(len=8) :: vstr
!EOP
!-------------------------------------------------------------------------
!BOC
!  Get the MPI version
   call MPI_GET_VERSION(ver,subver,ierr)
   if (ierr .ne. MPI_SUCCESS) THEN
      STDERR 'Fatal error: unable to get MPI version information.'
      call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
   end if
   write(vstr,'(a1,I1,a1,I1)') 'v',ver,'.',subver
   LEVEL0 "MPI is initialised - ",trim(vstr)

   LEVEL0 'Process ',myid,' of ',nprocs,' is alive on ',pname(1:len)
   end subroutine print_MPI_info
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_mesh_method - defines how calculation meshes are made
!
! !INTERFACE:
   SUBROUTINE set_mesh_method(method,reorder)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Sets global variables for calculation mesh and communication method.
!
! !INPUT PARAMTERS:
   integer, intent(in)                 :: method
   logical, intent(in)                 :: reorder
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-------------------------------------------------------------------------
!BOC
   mesh_method=method
   re_order=reorder
   end subroutine set_mesh_method
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_active_communicator - sets the active communicator to use.
!
! !INTERFACE:
   SUBROUTINE set_active_communicator(comm)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Sets the active communicator used in successive operations by setting
!  active\_comm to comm.
!
! !INPUT PARAMTERS:
   integer, intent(in)                 :: comm
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-------------------------------------------------------------------------
!BOC
   active_comm=comm
   end subroutine set_active_communicator
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_com_method - sets the communication method.
!
! !INTERFACE:
   SUBROUTINE set_com_method(method)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Set the communication method to method by setting private member
!  comm\_method to method.
!
! !INPUT PARAMTERS:
    integer, intent(in)                :: method
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-------------------------------------------------------------------------
!BOC
   if(nprocs .eq. 1) then
      comm_method=ONE_PROCESS
   else
      comm_method=method
   end if
   end subroutine set_com_method
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: part_domain_mpi - partition the calculation domain.
!
! !INTERFACE:
   SUBROUTINE part_domain_mpi(iextr,jextr,kmax,imin,imax,jmin,jmax,ioff,joff)
!$ use omp_lib
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Partition the calculation domain.
!
! !INPUT PARAMTERS:
   integer, intent(in)                 :: iextr,jextr,kmax
#ifdef STATIC
   integer, intent(in)                 :: imin,imax,jmin,jmax
#endif
!
! !OUTPUT PARAMTERS:
#ifndef STATIC
   integer, intent(out)                :: imin,imax,jmin,jmax
#endif
   integer, intent(out)                :: ioff,joff
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: i,j,zz(2),neighbours(nneighbours),numthreads
   character(len=*),parameter:: par_setup="par_setup.dat"
!EOP
!-------------------------------------------------------------------------
!BOC
   periods(1) = .false.
   periods(2) = .false.

   if(nprocs .eq. 1) then

#ifndef STATIC
      imin=1 ; imax=iextr ; jmin=1 ; jmax=jextr
#endif
#if 0
   STDERR 'AAAAkbk'
   STDERR imin,imax,jmin,jmax
   stop
#endif
      ioff=0 ; joff=0
      return
   end if

!  Set up the cartesian vitual topology
   select case (mesh_method)
      case(ONE_CELL)
#ifndef STATIC
         imin=1 ; imax = iextr ; jmin=1;  jmax=jextr
#endif
      case(ONED_MESH)
         if(iextr .gt. jextr) then
            dims(1) = 1
            dims(2) = nprocs
            com_direction = RIGHT_LEFT
         else
            dims(1) = nprocs
            dims(2) = 1
            com_direction = DOWN_UP
         end if
      case(TWOD_MESH)
#ifdef STATIC
         dims(1) = jextr/(jmax-jmin+1)
         if(mod(jextr,(jmax-jmin+1)) .ne. 0) dims(1) = dims(1)+1
         dims(2) = iextr/(imax-imin+1)
         if(mod(iextr,(imax-imin+1)) .ne. 0) dims(2) = dims(2)+1
#else
         dims(1) = nprocs
         dims(2) = 1
         dims(1)=0 ; dims(2)=0
         call MPI_DIMS_CREATE(nprocs, 2, dims, ierr)
!         STDERR 'dynamic dims ',dims
#endif
      case(MESH_FROM_FILE)
         call read_par_setup(par_setup,nprocs,myid,imax,jmax,iextr,jextr, &
                             ioff,joff,neighbours,numthreads)
         left  = neighbours(1) ; if (left  .eq. -1) left  = MPI_PROC_NULL
         ul    = neighbours(2) ; if (ul    .eq. -1) ul    = MPI_PROC_NULL
         up    = neighbours(3) ; if (up    .eq. -1) up    = MPI_PROC_NULL
         ur    = neighbours(4) ; if (ur    .eq. -1) ur    = MPI_PROC_NULL
         right = neighbours(5) ; if (right .eq. -1) right = MPI_PROC_NULL
         lr    = neighbours(6) ; if (lr    .eq. -1) lr    = MPI_PROC_NULL
         down  = neighbours(7) ; if (down  .eq. -1) down  = MPI_PROC_NULL
         ll    = neighbours(8) ; if (ll    .eq. -1) ll    = MPI_PROC_NULL

! IF we use OMP and IF the number of read threads is sensible (>0), then set #threads:
!$       if (numthreads>0) then
!$          LEVEL1 'Setting number of threads to ',numthreads
!$          call omp_set_num_threads(numthreads)
!$       end if
!$       LEVEL1 'Number of threads is ',omp_get_max_threads()
      case default
         FATAL 'A non valid partitioning method has been chosen'
         call MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
   end select

   if (mesh_method .eq. ONED_MESH .or. mesh_method .eq. TWOD_MESH) then
      call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periods,re_order,comm_hd,ierr)
      call MPI_CART_COORDS(comm_hd,myid,2,coords,ierr)
      call MPI_CART_SHIFT(comm_hd,0,1,down,up,ierr)
      call MPI_CART_SHIFT(comm_hd,1,1,left,right,ierr)

      zz(1) = coords(1)+1
      zz(2) = coords(2)+1
      if (zz(1) .gt. dims(1)-1 .or. zz(2) .gt. dims(2)-1) then
         ur = MPI_PROC_NULL
      else
         call MPI_CART_RANK(comm_hd,zz,ur,ierr)
      end if

      zz(1) = coords(1)+1
      zz(2) = coords(2)-1
      if (zz(1) .gt. dims(1)-1 .or. zz(2) .lt. 0) then
         ul  = MPI_PROC_NULL
      else
         call MPI_CART_RANK(comm_hd,zz,ul,ierr)
      end if

      zz(1) = coords(1)-1
      zz(2) = coords(2)-1
      if (zz(1) .lt. 0 .or. zz(2) .lt. 0) then
         ll  = MPI_PROC_NULL
      else
         call MPI_CART_RANK(comm_hd,zz,ll,ierr)
      end if

      zz(1) = coords(1)-1
      zz(2) = coords(2)+1
      if (zz(1) .lt. 0 .or. zz(2) .gt. dims(2)-1) then
         lr  = MPI_PROC_NULL
      else
         call MPI_CART_RANK(comm_hd,zz,lr,ierr)
      end if

#ifdef STATIC
      if(dims(2)*imax .lt. iextr) then
         FATAL 'problems - i - part_domain_mpi'
         stop
      end if
      if(dims(1)*jmax .lt. jextr) then
         FATAL 'problems - j - part_domain_mpi'
         stop
      end if
#else
      if(mod(iextr,dims(2)) .eq. 0) then
         imin=1;imax=iextr/dims(2)
      else
         imin=1;imax=iextr/dims(2)+1
      end if
      if(mod(jextr,dims(1)) .eq. 0) then
         jmin=1;jmax=jextr/dims(1)
      else
         jmin=1;jmax=jextr/dims(1)+1
      end if
#endif
      ioff=coords(2)*imax
      joff=coords(1)*jmax
   end if

   call MPI_data_types(imin,imax,jmin,jmax,kmax)

   STDERR LINE
   LEVEL2 'My id, coordinates, off-set and neighbours ....'
   STDERR LINE
   write(0,*) 'id = ',myid,coords,ioff,joff
   write(0,'(3(i3))') ul,up,ur
   write(0,'(3(i3))') left,myid,right
   write(0,'(3(i3))') ll,down,lr
   STDERR LINE

#ifdef DEBUG
   call MPI_BARRIER(comm_hd,ierr)
   if (myid .lt. 10) STDERR 'id = ',myid,dims,coords,left,right,down,up
   call MPI_BARRIER(comm_hd,ierr)
#endif
   return
   end subroutine part_domain_mpi
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: decomp_1d - decompose an array over m processors
!
! !INTERFACE:
   SUBROUTINE decomp_1d(n,np,m,s,e)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Decompoes an array over m processes
!
! !INPUT PARAMTERS:
  integer, intent(in)                  :: n,np,m
!
! !OUTPUT PARAMTERS:
  integer, intent(out)                 :: s,e
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: nlocal,deficit
!EOP
!-------------------------------------------------------------------------
!BOC
   nlocal  = n / np
   s       = m * nlocal + 1
   deficit = mod(n,np)
   s       = s + min(m,deficit)
   if (m .lt. deficit) THEN
       nlocal = nlocal + 1
   end if
   e = s + nlocal - 1
   if (e .gt. n .or. m .eq. np-1) e = n

   end subroutine decomp_1d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: MPI_data_types - defines misc. data types.
!
! !INTERFACE:
   SUBROUTINE MPI_data_types(imin,imax,jmin,jmax,kmax)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Sets up a number of convinient data types for passing various sub-sections
!  of 2D and 3D fields - when MPI-2 is out one should probably use
!  MPI\_TYPE\_CREATE\_SUBARRAY instead.
!
! !INPUT PARAMTERS:
   integer, intent(in)                 :: imin,imax,jmin,jmax,kmax
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: m,n,o
   integer                   :: real_extent
#ifdef _MPI_TYPE_EXTENT_
   INTEGER                   :: lower_bound,sizeof_realtype
#else
   INTEGER(KIND=MPI_ADDRESS_KIND) :: idum1,idum2
   INTEGER(KIND=MPI_ADDRESS_KIND) :: lower_bound,sizeof_realtype
#endif
!EOP
!-------------------------------------------------------------------------
!BOC
!kbk whatch out for this
   m = imax - imin + 1
   n = jmax - jmin + 1
   o = kmax + 1

   x_size = m+2*HALO
   y_size = n+2*HALO
   z_size = kmax+1
   xy_size = x_size*y_size
   xz_size = x_size*z_size
   yz_size = y_size*z_size
   xyz_size = x_size*y_size*z_size

!  Set up different data types
   call MPI_TYPE_SIZE(MPI_REALTYPE,real_extent,ierr)
#ifdef _MPI_TYPE_EXTENT_
!  deprecated - use MPI_TYPE_GET_EXTENT instead - does not work yet
   call MPI_TYPE_EXTENT(MPI_REALTYPE,sizeof_realtype,ierr)
#else
   call MPI_TYPE_GET_EXTENT(MPI_REALTYPE,idum1,idum2,ierr)
   lower_bound     = idum1
   sizeof_realtype = idum2
#endif

!  1 x-line
   call MPI_TYPE_VECTOR(m,1,1,MPI_REALTYPE,x_line,ierr)
   call MPI_TYPE_COMMIT(x_line,ierr)

!  HALO x-lines
   call MPI_TYPE_HVECTOR(HALO,1,x_size*sizeof_realtype,x_line,x_lines,ierr)
   call MPI_TYPE_COMMIT(x_lines,ierr)

!  1 y-line
   call MPI_TYPE_VECTOR(n,1,x_size,MPI_REALTYPE,y_line,ierr)
   call MPI_TYPE_COMMIT(y_line,ierr)

!  HALO y-lines
   call MPI_TYPE_HVECTOR(HALO,1,1*sizeof_realtype,y_line,y_lines,ierr)
   call MPI_TYPE_COMMIT(y_lines,ierr)

!  1 HALO-line
   call MPI_TYPE_VECTOR(HALO,1,1,MPI_REALTYPE,halo_line,ierr)
   call MPI_TYPE_COMMIT(halo_line,ierr)

!  HALO square
   call MPI_TYPE_HVECTOR(HALO,1,x_size*sizeof_realtype,halo_line,halo_square,ierr)
   call MPI_TYPE_COMMIT(halo_square,ierr)

!  1 xz-slice
   call MPI_TYPE_HVECTOR(o,1,xy_size*sizeof_realtype,x_line,xz_slice,ierr)
   call MPI_TYPE_COMMIT(xz_slice,ierr)

!  HALO xz-slices
   call MPI_TYPE_HVECTOR(o,1,xy_size*sizeof_realtype,x_lines,xz_slices,ierr)
   call MPI_TYPE_COMMIT(xz_slices,ierr)

#if 1
!  1 yz slice
   call MPI_TYPE_HVECTOR(o,1,xy_size*sizeof_realtype,y_line,yz_slice,ierr)
   call MPI_TYPE_COMMIT(yz_slice,ierr)

!  HALO yz-slices
   call MPI_TYPE_HVECTOR(o,1,xy_size*sizeof_realtype,y_lines,yz_slices,ierr)
   call MPI_TYPE_COMMIT(yz_slices,ierr)

#else
!  Using MPI: p. 292
!  1 yz slice
!  first data type for a(i,:sy:ey:k)
   call MPI_TYPE_VECTOR(n,1,x_size,MPI_REALTYPE,i1_slice,ierr)
   call MPI_TYPE_COMMIT(i1_slice,ierr)

!  vector if i1_slices
   call MPI_TYPE_HVECTOR(o,1,xy_size*sizeof_realtype,i1_slice,yz_slice,ierr)
   call MPI_TYPE_COMMIT(yz_slice,ierr)

   call MPI_TYPE_HVECTOR(HALO,1,1*sizeof_realtype,yz_slice,yz_slices,ierr)
   call MPI_TYPE_COMMIT(yz_slices,ierr)
#endif

!  a vertical column
   call MPI_TYPE_VECTOR(o,1,xy_size,MPI_REALTYPE,z_column,ierr)
   call MPI_TYPE_COMMIT(z_column,ierr)

!  HALO square vertical columns
   call MPI_TYPE_HVECTOR(o,1,xy_size*sizeof_realtype,halo_square, &
                         halo_columns,ierr)
   call MPI_TYPE_COMMIT(halo_columns,ierr)

   size_point = sizeof_realtype
   if (left  .eq. MPI_PROC_NULL) then
      size_left = 0
   else
      size_left = n*sizeof_realtype
   end if

   if (ul    .eq. MPI_PROC_NULL) then
      size_ul = 0
   else
      size_ul = sizeof_realtype
   end if

   if (up    .eq. MPI_PROC_NULL) then
      size_up = 0
   else
      size_up = m*sizeof_realtype
   end if

   if (ur    .eq. MPI_PROC_NULL) then
      size_ur = 0
   else
      size_ur = sizeof_realtype
   end if

   if (right .eq. MPI_PROC_NULL) then
      size_right = 0
   else
      size_right = n*sizeof_realtype
   end if

   if (lr    .eq. MPI_PROC_NULL) then
      size_lr = 0
   else
      size_lr = sizeof_realtype
   end if

   if (down  .eq. MPI_PROC_NULL) then
      size_down = 0
   else
      size_down = m*sizeof_realtype
   end if

   if (ll    .eq. MPI_PROC_NULL) then
      size_ll = 0
   else
      size_ll = sizeof_realtype
   end if

   return
   end subroutine MPI_data_types
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: update_2d_halo_mpi - updates the halo zones for 2D fields.
!
! !INTERFACE:
   subroutine update_2d_halo_mpi(f1,f2,imin,jmin,imax,jmax,tag,mirror)
   use getm_timers, only: tic, toc, TIM_HALO2D
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Print information on the MPI environment
!
! !INPUT PARAMTERS:
   integer, intent(in)                 :: imin,jmin,imax,jmax
   integer, intent(in)                 :: tag
   logical, optional, intent(in)       :: mirror
!
! !INPUT/OUTPUT PARAMTERS:
   REALTYPE, intent(inout), dimension(E2DFIELD):: f1,f2
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: il,jl,ih,jh
   logical                   :: do_mirror=.true.
!EOP
!-------------------------------------------------------------------------
!BOC
   if (last_action .ne. WAITING) then
      FATAL 'Last action was not WAITING - not ready for sending (2D)'
      call MPI_ABORT(active_comm,-1,ierr)
   end if
   call tic(TIM_HALO2D)
   il=imin;ih=imax;jl=jmin;jh=jmax

   if ( present(mirror) ) do_mirror = mirror

   select case (comm_method)
      case(ONE_PROCESS)
         if (do_mirror) then
         f1(il-1, : )  = f2(il, :  )
         f1(ih+1, : )  = f2(ih, :  )
         f1( :, jl-1 ) = f2( :, jl )
         f1( :, jh+1 ) = f2( :, jh )
         f1(il-1,jh+1) = f2(il,jh)
         f1(ih+1,jh+1) = f2(ih,jh)
         f1(ih+1,jl-1) = f2(ih,jl)
         f1(il-1,jl-1) = f2(il,jl)
         end if
      case(ONED_SENDRECV)
         if(com_direction .eq. RIGHT_LEFT) then
#ifdef DEBUG
STDERR 'ONED_SENDRECV - y_lines'
#endif
            call MPI_SENDRECV(f1(il,jl),   1, y_line, left , tag, &
                              f2(ih+1,jl), 1, y_line, right, tag, &
                              active_comm, status, ierr)

            call MPI_SENDRECV(f1(ih,jl),   1, y_line, right, tag, &
                              f2(il-1,jl), 1, y_line, left , tag, &
                              active_comm, status, ierr)
         else
#ifdef DEBUG
STDERR 'ONED_SENDRECV - x_lines'
#endif
            call MPI_SENDRECV(f1(il,jh),   1, x_line, up,   tag, &
                              f2(il,jl-1), 1, x_line, down, tag, &
                              active_comm, status, ierr)
            call MPI_SENDRECV(f1(il,jl),   1, x_line, down, tag, &
                              f2(il,jh+1), 1, x_line, up  , tag, &
                              active_comm, status, ierr)
         end if
      case(ONED_NONBLOCKING)
         if(com_direction .eq. RIGHT_LEFT) then
#ifdef DEBUG
STDERR 'ONED_NONBLOCKING - y_lines'
#endif
            call MPI_IRECV(f2(il-HALO,jl),     1, y_lines, left,  tag, &
                              active_comm, req(1), ierr)
            call MPI_IRECV(f2(ih+1,jl),        1, y_lines, right, tag, &
                              active_comm, req(2), ierr)
            call MPI_ISEND(f1(il,jl),          1, y_lines, left,  tag, &
                              active_comm, req(3), ierr)
            call MPI_ISEND(f1(ih-(HALO-1),jl), 1, y_lines, right, tag, &
                              active_comm, req(4), ierr)
         else
#ifdef DEBUG
STDERR 'ONED_NONBLOCKING - x_lines'
#endif
            call MPI_IRECV(f2(il,jl-HALO),     1, x_lines, down, tag, &
                              active_comm, req(1), ierr)
            call MPI_IRECV(f2(il,jh+1),        1, x_lines, up,   tag, &
                              active_comm, req(2), ierr)
            call MPI_ISEND(f1(il,jl),          1, x_lines, down, tag, &
                              active_comm, req(3), ierr)
            call MPI_ISEND(f1(il,jh-(HALO-1)), 1, x_lines, up,   tag, &
                              active_comm, req(4), ierr)
         end if
      case(TWOD_SENDRECV)
#ifdef DEBUG
STDERR 'TWOD_SENDRECV'
#endif
         call MPI_SENDRECV(f1(il,jh),   1, x_line, up  , tag, &
                           f2(il,jl-1), 1, x_line, down, tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(il,jl),   1, x_line, down, tag, &
                           f2(il,jh+1), 1, x_line, up  , tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(ih,jl),   1, y_line, right, tag, &
                           f2(il-1,jl), 1, y_line, left,  tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(il,jl),   1, y_line, left , tag, &
                           f2(ih+1,jl), 1, y_line, right, tag, &
                           active_comm, status, ierr)
!        Corner points
         call MPI_SENDRECV(f1(il,jl),    1, MPI_REALTYPE, ll, tag, &
                           f2(ih+1,jh+1),1, MPI_REALTYPE, ur,tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(ih,jl),    1, MPI_REALTYPE, lr,tag, &
                           f2(il-1,jh+1),1, MPI_REALTYPE, ul, tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(ih,jh),    1, MPI_REALTYPE, ur,tag, &
                           f2(il-1,jl-1),1, MPI_REALTYPE, ll, tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(il,jh),    1, MPI_REALTYPE, ul,tag, &
                           f2(ih+1,jl-1),1, MPI_REALTYPE, lr,tag, &
                           active_comm, status, ierr)
      case(TWOD_NONBLOCKING)
!        Recieving x_lines
         call MPI_IRECV(f2(il,jl-HALO), 1, x_lines, down,  tag, &
                           active_comm, req(1), ierr)
         call MPI_IRECV(f2(il,jh+1),    1, x_lines, up,    tag, &
                           active_comm, req(2), ierr)

!        Recieving y_lines
         call MPI_IRECV(f2(il-HALO,jl), 1, y_lines, left,  tag, &
                           active_comm, req(3), ierr)
         call MPI_IRECV(f2(ih+1,jl),    1, y_lines, right, tag, &
                           active_comm, req(4), ierr)

!        Recieving HALOxHALO corner squares
         call MPI_IRECV(f2(il-HALO,jl-HALO), 1, halo_square,ll,tag, &
                           active_comm, req(5), ierr)
         call MPI_IRECV(f2(ih+1,jl-HALO), 1, halo_square,lr,tag, &
                           active_comm, req(6), ierr)
         call MPI_IRECV(f2(ih+1,jh+1), 1, halo_square,ur,tag, &
                           active_comm, req(7), ierr)
         call MPI_IRECV(f2(il-HALO,jh+1), 1, halo_square,ul,tag, &
                           active_comm, req(8), ierr)

!        Sending x_lines
         call MPI_ISEND(f1(il,jl),          1, x_lines, down,  tag, &
                           active_comm, req(9), ierr)
         call MPI_ISEND(f1(il,jh-(HALO-1)), 1, x_lines, up,    tag, &
                           active_comm, req(10), ierr)

!        Sending y_lines
         call MPI_ISEND(f1(il,jl),          1, y_lines, left,  tag, &
                           active_comm, req(11), ierr)
         call MPI_ISEND(f1(ih-(HALO-1),jl), 1, y_lines, right, tag, &
                           active_comm, req(12), ierr)

!        Sending HALOxHALO corner squares
         call MPI_ISEND(f1(ih-(HALO-1),jh-(HALO-1)), 1, halo_square, ur,tag, &
                           active_comm, req(13), ierr)

         call MPI_ISEND(f1(il,jh-(HALO-1)), 1, halo_square, ul,tag, &
                           active_comm, req(14), ierr)

         call MPI_ISEND(f1(il,jl), 1, halo_square, ll,tag, &
                           active_comm, req(15), ierr)

         call MPI_ISEND(f1(ih-(HALO-1),jl), 1, halo_square, lr,tag, &
                           active_comm, req(16), ierr)

         all_2d_exchange =  all_2d_exchange      &
                          + HALO*size_left       &
                          + HALO*size_up         &
                          + HALO*size_right      &
                          + HALO*size_down       &
                          + HALO*HALO*size_point &
                          * (size_ul+size_ur+size_lr+size_ll)

      case default
         FATAL 'A non valid communication method has been chosen'
         stop 'update_2d_halo_mpi'
   end select

   if (do_mirror) then
      if ( comm_method .ne. ONE_PROCESS ) then
         if(left  .eq. MPI_PROC_NULL) f1(il-1, : )  = f1(il, : )
         if(right .eq. MPI_PROC_NULL) f1(ih+1, : )  = f1(ih, : )
         if(down  .eq. MPI_PROC_NULL) f1( :, jl-1)  = f1( :, jl)
         if(up    .eq. MPI_PROC_NULL) f1( :, jh+1)  = f1( :, jh)
         if(ul    .eq. MPI_PROC_NULL) f1(il-1,jh+1) = f1(il,jh)
         if(ur    .eq. MPI_PROC_NULL) f1(ih+1,jh+1) = f1(ih,jh)
         if(lr    .eq. MPI_PROC_NULL) f1(ih+1,jl-1) = f1(ih,jl)
         if(ll    .eq. MPI_PROC_NULL) f1(il-1,jl-1) = f1(il,jl)
      end if
   end if

   last_action = SENDING

   call toc(TIM_HALO2D)
   return
   end subroutine update_2d_halo_mpi
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: update_3d_halo_mpi - updates the halo zones for 3D fields.
!
! !INTERFACE:
   SUBROUTINE update_3d_halo_mpi(f1,f2,imin,jmin,imax,jmax,kmax,tag,mirror)
   use getm_timers, only: tic, toc, TIM_HALO3D
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Updates the halo zones for 3D fields.
!
! !INPUT PARAMTERS:
   integer, intent(in)                 :: imin,jmin,imax,jmax,kmax
   integer, intent(in)                 :: tag
   logical, intent(in)                 :: mirror
!
! !INPUT/OUTPUT PARAMTERS:
   REALTYPE, intent(inout), DIMENSION(I3DFIELD) :: f1,f2
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: il,jl,ih,jh
!EOP
!-------------------------------------------------------------------------
!BOC
   if (last_action .ne. WAITING) then
      FATAL 'Last action was not WAITING - not ready for sending (3D)'
      call MPI_ABORT(active_comm,-1,ierr)
   end if
   call tic(TIM_HALO3D)
   il=imin;ih=imax;jl=jmin;jh=jmax
   select case (comm_method)
      case(ONE_PROCESS)
         if ( mirror ) then
         f1(il-1, :, : )  = f2(il, :, :  )
         f1(ih+1, :, : )  = f2(ih, :, :  )
         f1( :, jl-1, : ) = f2( :, jl, : )
         f1( :, jh+1, : ) = f2( :, jh, : )
         end if
      case(ONED_SENDRECV)
         if(com_direction .eq. RIGHT_LEFT) then
#ifdef DEBUG
STDERR 'ONED_SENDRECV - yz_slices'
#endif
            call MPI_SENDRECV(f1(il,jl,0),   1, yz_slice, left , tag, &
                              f2(ih+1,jl,0), 1, yz_slice, right, tag, &
                              active_comm, status, ierr)
            call MPI_SENDRECV(f1(ih,jl,0),   1, yz_slice, right, tag, &
                              f2(il-1,jl,0), 1, yz_slice, left , tag, &
                              active_comm, status, ierr)
         else
#ifdef DEBUG
STDERR 'ONED_SENDRECV - xz_slices'
#endif
            call MPI_SENDRECV(f1(il,jl,0),   1, xz_slice, down, tag, &
                              f2(il,jh+1,0), 1, xz_slice, up  , tag, &
                              active_comm, status, ierr)
            call MPI_SENDRECV(f1(il,jh,0),   1, xz_slice, up,   tag, &
                              f2(il,jl-1,0), 1, xz_slice, down, tag, &
                              active_comm, status, ierr)
         end if
      case(ONED_NONBLOCKING)
         if(com_direction .eq. RIGHT_LEFT) then
#ifdef DEBUG
STDERR 'ONED_NONBLOCKING - yz_slices'
#endif
            call MPI_IRECV(f2(il-HALO,jl,0),     1, yz_slices, left,  tag, &
                              active_comm, req(1), ierr)
            call MPI_IRECV(f2(ih+1,jl,0),        1, yz_slices, right, tag, &
                              active_comm, req(2), ierr)
            call MPI_ISEND(f1(il,jl,0),          1, yz_slices, left,  tag, &
                              active_comm, req(3), ierr)
            call MPI_ISEND(f1(ih-(HALO-1),jl,0), 1, yz_slices, right, tag, &
                              active_comm, req(4), ierr)
         else
#ifdef DEBUG
STDERR 'ONED_NONBLOCKING - xz_slices'
#endif
            call MPI_IRECV(f2(il,jl-HALO,0),     1, xz_slices, down, tag, &
                              active_comm, req(1), ierr)
            call MPI_IRECV(f2(il,jh+1,0),        1, xz_slices, up,   tag, &
                              active_comm, req(2), ierr)
            call MPI_ISEND(f1(il,jl,0),          1, xz_slices, down, tag, &
                              active_comm, req(3), ierr)
            call MPI_ISEND(f1(il,jh-(HALO-1),0), 1, xz_slices, up, tag,   &
                              active_comm, req(4), ierr)
         end if
      case(TWOD_SENDRECV)
#ifdef DEBUG
STDERR 'TWOD_SENDRECV'
#endif
         call MPI_SENDRECV(f1(il,jh,0),   1, xz_slice, up  , tag, &
                           f2(il,jl-1,0), 1, xz_slice, down, tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(il,jl,0),   1, xz_slice, down, tag, &
                           f2(il,jh+1,0), 1, xz_slice, up  , tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(ih,jl,0),   1, yz_slice, right, tag, &
                           f2(il-1,jl,0), 1, yz_slice, left,  tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(il,jl,0),   1, yz_slice, left , tag, &
                           f2(ih+1,jl,0), 1, yz_slice, right, tag, &
                           active_comm, status, ierr)
!        Corner points
         call MPI_SENDRECV(f1(il,jl,0),    1, z_column, ll, tag, &
                           f2(ih+1,jh+1,0),1, z_column, ur, tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(ih,jl,0),    1, z_column, lr, tag, &
                           f2(il-1,jh+1,0),1, z_column, ul, tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(ih,jh,0),    1, z_column, ur, tag, &
                           f2(il-1,jl-1,0),1, z_column, ll, tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(il,jh,0),    1, z_column, ul, tag, &
                           f2(ih+1,jl-1,0),1, z_column, lr, tag, &
                           active_comm, status, ierr)
      case(TWOD_NONBLOCKING)
#ifdef DEBUG
STDERR 'TWOD_NONBLOCKING'
#endif
!        Recieving xz_slices
         call MPI_IRECV(f2(il,jl-HALO,0), 1, xz_slices, down, tag, &
                           active_comm, req(1), ierr)
         call MPI_IRECV(f2(il,jh+1,0),    1, xz_slices, up,   tag, &
                           active_comm, req(2), ierr)

!        Recieving yz_slices
         call MPI_IRECV(f2(il-HALO,jl,0), 1, yz_slices, left,  tag, &
                           active_comm, req(3), ierr)
         call MPI_IRECV(f2(ih+1,jl,0),    1, yz_slices, right, tag, &
                           active_comm, req(4), ierr)

!        Recieving corner columns
         call MPI_IRECV(f2(il-HALO,jl-HALO,0), 1, halo_columns,ll,tag, &
                           active_comm, req(5), ierr)
         call MPI_IRECV(f2(ih+1,jl-HALO,0), 1, halo_columns,lr,tag, &
                           active_comm, req(6), ierr)
         call MPI_IRECV(f2(ih+1,jh+1,0), 1, halo_columns,ur,tag, &
                           active_comm, req(7), ierr)
         call MPI_IRECV(f2(il-HALO,jh+1,0), 1, halo_columns,ul,tag, &
                           active_comm, req(8), ierr)

!        Sending xz_slices
         call MPI_ISEND(f1(il,jl,0),          1, xz_slices, down,  tag, &
                           active_comm,  req(9), ierr)
         call MPI_ISEND(f1(il,jh-(HALO-1),0), 1, xz_slices, up,    tag, &
                           active_comm, req(10), ierr)

!        Sending yz_slices
         call MPI_ISEND(f1(il,jl,0),          1, yz_slices, left,  tag, &
                           active_comm, req(11), ierr)
         call MPI_ISEND(f1(ih-(HALO-1),jl,0), 1, yz_slices, right, tag, &
                           active_comm, req(12), ierr)

!        Sending corner columns
         call MPI_ISEND(f1(ih-(HALO-1),jh-(HALO-1),0), 1, halo_columns, ur, tag, &
                           active_comm, req(13), ierr)

         call MPI_ISEND(f1(il,jh-(HALO-1),0), 1, halo_columns, ul, tag, &
                           active_comm, req(14), ierr)

         call MPI_ISEND(f1(il,jl,0), 1, halo_columns, ll, tag, &
                           active_comm, req(15), ierr)

         call MPI_ISEND(f1(ih-(HALO-1),jl,0), 1, halo_columns, lr, tag, &
                           active_comm, req(16), ierr)

         all_3d_exchange =  all_3d_exchange      &
                          + (kmax+1)*            &
                          ( HALO*size_left       &
                          + HALO*size_up         &
                          + HALO*size_right      &
                          + HALO*size_down       &
                          + HALO*HALO*size_point &
                          * (size_ul+size_ur+size_lr+size_ll) &
                          )

      case default
         FATAL 'A non valid communication method has been chosen'
         stop 'update_3d_halo_mpi'
   end select

! Produces an error in some sub-domain layout cases if the following is
! included.
   if ( mirror ) then
   if ( comm_method .ne. ONE_PROCESS ) then
      if(left  .eq. MPI_PROC_NULL) f1(il-1, :, : )   = f1(il, :, : )
      if(right .eq. MPI_PROC_NULL) f1(ih+1, :, : )   = f1(ih, :, : )
      if(down  .eq. MPI_PROC_NULL) f1( :, jl-1, : )  = f1( :, jl, : )
      if(up    .eq. MPI_PROC_NULL) f1( :, jh+1, : )  = f1( :, jh, : )
      if(ul    .eq. MPI_PROC_NULL) f1(il-1,jh+1, : ) = f1(il,jh, : )
      if(ur    .eq. MPI_PROC_NULL) f1(ih+1,jh+1, : ) = f1(ih,jh, : )
      if(lr    .eq. MPI_PROC_NULL) f1(ih+1,jl-1, : ) = f1(ih,jl, : )
      if(ll    .eq. MPI_PROC_NULL) f1(il-1,jl-1, : ) = f1(il,jl, : )
   end if
   end if
   last_action = SENDING
   call toc(TIM_HALO3D)
   return
   end subroutine update_3d_halo_mpi
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: wait_halo_mpi - wait for any un-finished communications
!
! !INTERFACE:
   SUBROUTINE wait_halo_mpi(tag)
   use getm_timers, only: tic, toc, TIM_HALOWAIT
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Call MPI\_WAITALL to wait for any un-finished communications. If SENDRECV
!  communications are used this call has no effect.
!
! !INPUT PARAMTERS:
   integer, intent(in)                 :: tag
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-------------------------------------------------------------------------
!BOC
   if (last_action .ne. SENDING) then
      FATAL 'Last action was not sending - nothing to wait for'
      call MPI_ABORT(active_comm,-1,ierr)
   end if
   call tic(TIM_HALOWAIT)
   select case (comm_method)
      case(ONE_PROCESS)
      case(ONED_SENDRECV)
      case(ONED_NONBLOCKING)
         !Waiting for 2 sends and 2 recieves.
         call MPI_WAITALL (4,req,status_array,ierr)
      case(TWOD_SENDRECV)
      case(TWOD_NONBLOCKING)
         if(last_action .ne. SENDING) then
         end if
         !Waiting for 8 sends and 8 recieves.
         call MPI_WAITALL (16,req,status_array,ierr)
      case default
         FATAL 'A non valid communication method has been chosen'
         stop 'wait_mpi'
   end select
   last_action = WAITING
   call toc(TIM_HALOWAIT)
   return
   end subroutine wait_halo_mpi
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_flag_mpi - sets and checks integer flag
!
! !INTERFACE:
   SUBROUTINE set_flag_mpi(n,flag,flags)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Call MPI\_WAITALL to wait for any un-finished communications. If SENDRECV
!  communications are used this call has no effect.
!
! !INPUT PARAMTERS:
   integer, intent(in)                 :: n,flag
!
! !OUTPUT PARAMTERS:
   integer, intent(out)                :: flags(n)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Peter Holtermann
!
!EOP
!-------------------------------------------------------------------------
!BOC

   CALL MPI_GATHER(flag,1,MPI_INTEGER,flags,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr);
   CALL MPI_BCAST(flags,nprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

   return
   end subroutine set_flag_mpi
!EOC

!-----------------------------------------------------------------------

   end module halo_mpi

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
