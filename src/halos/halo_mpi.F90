!$Id: halo_mpi.F90,v 1.1.1.1 2002-05-02 14:01:30 gotm Exp $
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
!  If \#define STATIC - imin,imax,jmin,jmax are input parameters else output
!
! !USES:
   use mpi
   IMPLICIT NONE
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public	:: init_mpi,print_MPI_info,set_active_communicator
   public	:: update_2d_halo_mpi,update_3d_halo_mpi,wait_halo_mpi
   private	:: part_domain
   private	:: MPI_data_types

! !PUBLIC DATA MEMBERS:
   integer, public	:: myid, nprocs
   integer, public	:: comm_hd=MPI_COMM_WORLD
!   integer, public	:: comm_wave=MPI_COMM_WORLD
!   integer, public	:: comm_biology=MPI_COMM_WORLD
!
!  !DEFINED PARAMETERS:
   integer, parameter		:: neighbors=8
!  Different mesh specification methods
   integer, parameter		:: ONE_CELL=-1
   integer, parameter		:: ONED_MESH=0
   integer, parameter		:: TWOD_MESH=1
   integer, parameter		:: MESH_FROM_FILE=2
!  Methods of communication
   integer, parameter		:: ONE_PROCESS=-1
   integer, parameter		:: ONED_SENDRECV=0
   integer, parameter		:: ONED_NONBLOCKING=1
   integer, parameter		:: TWOD_SENDRECV=2
   integer, parameter		:: TWOD_NONBLOCKING=3
!  Direction in case of ONED_? communications
   integer, parameter		:: RIGHT_LEFT=1
   integer, parameter		:: DOWN_UP=2
!  Last action
!   integer, parameter		:: SEND_1D=1
!   integer, parameter		:: SEND_2D=2
   integer, parameter		:: SENDING=1
   integer, parameter		:: WAITING=2
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: halo_mpi.F90,v $
!  Revision 1.1.1.1  2002-05-02 14:01:30  gotm
!  recovering after CVS crash
!
!
! !LOCAL VARIABLES:
   character(LEN = 256), private	:: pname
   integer		:: active_comm=MPI_COMM_WORLD
   integer		:: mesh_method=ONE_CELL
   logical		:: re_order=.false.
   integer		:: comm_method=ONE_PROCESS
   integer		:: len
   integer		:: ierr
   integer		:: x_line,y_line
   integer		:: xz_slice,yz_slice,i1_slice
   integer		:: z_column
   integer		:: x_size,y_size,z_size
   integer		:: xy_size,xz_size,yz_size,xyz_size
   integer		:: com_direction
   integer 		:: req(2*neighbors)
   integer		:: status_array(MPI_STATUS_SIZE,2*neighbors)
   integer		:: dims(2),coords(2)
   logical		:: periods(2)
   integer		:: up,down,left,right
   integer		:: uright,uleft,lleft,lright
   integer		:: status(MPI_STATUS_SIZE)
   integer		:: last_action=WAITING
!EOP
!-----------------------------------------------------------------------
!BOC

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_mpi - initialize MPI environment
!
! !INTERFACE:
   SUBROUTINE init_mpi(iextr,jextr,imin,imax,jmin,jmax,kmax)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Initialize Parallel environment
!
! !INPUT PARAMTERS:
   integer, intent(in)	:: iextr,jextr,kmax
#ifdef STATIC
   integer, intent(in)	:: imin,imax,jmin,jmax
#endif
!
! !INPUT/OUTPUT PARAMTERS:
!
! !OUTPUT PARAMTERS:
#ifndef STATIC
   integer, intent(out)	:: imin,imax,jmin,jmax
#endif
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer	:: MeshMethod,MsgMethod
   logical	:: reorder
   namelist /nampar/ MeshMethod,reorder,MsgMethod
!
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   write(debug,*) 'init_mpi'
#endif

   LEVEL1 'OK - we are running in parallel'
   LEVEL1 'So we have to do a bit more initialization'

!  Read parallel/MPI specific things from the namelist.
   open(NAMLST,file='/home/kbk/mpi_test/parallel.inp')
   read(NAMLST,nampar)
   close(NAMLST)

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

   call set_com_method(MsgMethod)
   call set_mesh_method(MeshMethod,reorder)

!  Get the processor names
   call MPI_GET_PROCESSOR_NAME(pname,len,ierr)
   if(ierr .ne. MPI_SUCCESS) THEN
     STDERR 'Fatal error: unable to get processor name.'
     call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
   end if

   call part_domain(iextr,jextr,comm_hd,imin,imax,jmin,jmax)

   call MPI_data_types(imin,imax,jmin,jmax,kmax)

   write(*,*) '-------------------------------'
   write(*,*) 'id = ',myid,dims,coords
   write(*,'(3(i3))') uleft,up,uright
   write(*,'(3(i3))') left,myid,right
   write(*,'(3(i3))') lleft,down,lright
   write(*,*) '-------------------------------'

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
! !INPUT PARAMTERS:
!
! !INPUT/OUTPUT PARAMTERS:
!
! !OUTPUT PARAMTERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
!
!EOP
!-------------------------------------------------------------------------
!BOC
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
! !INPUT PARAMTERS:
!
! !INPUT/OUTPUT PARAMTERS:
!
! !OUTPUT PARAMTERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
!EOP
!-------------------------------------------------------------------------
!BOC
   STDERR 'Process ',myid,' of ',nprocs,' is alive on ',pname(1:len)
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
   integer, INTENT(IN)	:: method
   logical, INTENT(IN)	:: reorder
!
! !INPUT/OUTPUT PARAMTERS:
!
! !OUTPUT PARAMTERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
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
   integer, INTENT(IN)		:: comm
!
! !INPUT/OUTPUT PARAMTERS:
!
! !OUTPUT PARAMTERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
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
    integer, INTENT(IN)		:: method
!
! !INPUT/OUTPUT PARAMTERS:
!
! !OUTPUT PARAMTERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
!EOP
!-------------------------------------------------------------------------
!BOC
   if(nprocs .eq. 1) THEN
      comm_method=ONE_PROCESS
   ELSE
      comm_method=method
   end if
   end subroutine set_com_method
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: part_domain - partition the calculation domain.
!
! !INTERFACE:
   SUBROUTINE part_domain(iextr,jextr,comm,imin,imax,jmin,jmax)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Partition the calculation domain.
!
! !INPUT PARAMTERS:
   integer, INTENT(IN)		:: iextr,jextr
#ifdef STATIC
   integer, INTENT(in)		:: imin,imax,jmin,jmax
#endif
!
! !INPUT/OUTPUT PARAMTERS:
!
! !OUTPUT PARAMTERS:
   integer, INTENT(OUT)		:: comm
#ifndef STATIC
   integer, INTENT(OUT)		:: imin,imax,jmin,jmax
#endif
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer	:: i,j,zz(2)
!EOP
!-------------------------------------------------------------------------
!BOC
   periods(1) = .false.
   periods(2) = .false.

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
	 STDERR 'dims ',dims
#else
         dims(1) = nprocs
         dims(2) = 1
         dims(1)=0 ; dims(2)=0
         call MPI_DIMS_CREATE(nprocs, 2, dims, ierr)
#endif
      case(MESH_FROM_FILE)
         call MPI_COMM_DUP(MPI_COMM_WORLD,comm,ierr)
         FATAL 'MESH_FROM_FILE not implemented yet'
         call MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
      case default
         FATAL 'A non valid partitioning method has been chosen'
         call MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
    end select

    if (mesh_method .EQ. ONED_MESH .OR. mesh_method .EQ. TWOD_MESH) THEN
       call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periods,re_order,comm,ierr)
       call MPI_CART_COORDS(comm,myid,2,coords,ierr)
       call MPI_CART_SHIFT(comm,0,1,down,up,ierr)
       call MPI_CART_SHIFT(comm,1,1,left,right,ierr)

       uright = MPI_PROC_NULL
       lright = MPI_PROC_NULL
       uleft  = MPI_PROC_NULL
       lleft  = MPI_PROC_NULL

       zz(1) = coords(1)+1
       zz(2) = coords(2)+1
       call MPI_CART_RANK(comm,zz,uright,ierr)

       zz(1) = coords(1)+1
       zz(2) = coords(2)-1
       call MPI_CART_RANK(comm,zz,uleft,ierr)

       zz(1) = coords(1)-1
       zz(2) = coords(2)-1
       call MPI_CART_RANK(comm,zz,lleft,ierr)

       zz(1) = coords(1)-1
       zz(2) = coords(2)+1
       call MPI_CART_RANK(comm,zz,lright,ierr)
#ifndef STATIC
!       call decomp_1d(iextr,dims(2),coords(2),imin,imax)
!       call decomp_1d(jextr,dims(1),coords(1),jmin,jmax)
#endif
    end if

#ifdef DEBUG
   call MPI_BARRIER(comm,ierr)
   if (myid .lt. 10) STDERR 'id = ',myid,dims,coords,left,right,down,up
   call MPI_BARRIER(comm,ierr)
#endif
   return
   end subroutine part_domain
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
  integer, INTENT(IN)	:: n,np,m
!
! !INPUT/OUTPUT PARAMTERS:
!
! !OUTPUT PARAMTERS:
  integer, INTENT(OUT)	:: s,e
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer		:: nlocal,deficit
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
   integer, INTENT(IN)	:: imin,imax,jmin,jmax,kmax
!
! !INPUT/OUTPUT PARAMTERS:
!
! !OUTPUT PARAMTERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer		:: m,n,o,sizeof
   integer		:: real_extent
!
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

!  left-right strip
   call MPI_TYPE_VECTOR(1,m,1,MPI_REALTYPE,x_line,ierr)
   call MPI_TYPE_COMMIT(x_line,ierr)

!  down-up strip
   call MPI_TYPE_VECTOR(n,1,x_size,MPI_REALTYPE,y_line,ierr)
   call MPI_TYPE_COMMIT(y_line,ierr)

!  down-up slice
   call MPI_TYPE_VECTOR(o,m,xy_size,MPI_REALTYPE,xz_slice,ierr)
   call MPI_TYPE_COMMIT(xz_slice,ierr)

!  Using MPI: p. 292
!  left-right slice
!  first data type for a(i,:sy:ey:k)
   call MPI_TYPE_VECTOR(n,1,x_size,MPI_REALTYPE,i1_slice,ierr)
   call MPI_TYPE_COMMIT(i1_slice,ierr)
!  vector if i1_slices
   call MPI_TYPE_EXTENT(MPI_REALTYPE,sizeof,ierr)
   call MPI_TYPE_HVECTOR(o,1,xy_size*sizeof,i1_slice,yz_slice,ierr)
   call MPI_TYPE_COMMIT(yz_slice,ierr)

!  a vertical column
   call MPI_TYPE_VECTOR(o,1,xy_size,MPI_REALTYPE,z_column,ierr)
   call MPI_TYPE_COMMIT(z_column,ierr)

   return
   end subroutine MPI_data_types
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: update_2d_halo_mpi - updates the halo zones for 2D fields.
!
! !INTERFACE:
   SUBROUTINE update_2d_halo_mpi(f1,f2,imin,jmin,imax,jmax,tag)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Print information on the MPI environment
!
! !INPUT PARAMTERS:
   integer, INTENT(IN)	:: imin,jmin,imax,jmax
   integer, INTENT(IN)	:: tag
!
! !INPUT/OUTPUT PARAMTERS:
   REALTYPE, INTENT(INOUT), DIMENSION(E2DFIELD)	:: f1,f2
!
! !OUTPUT PARAMTERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
!EOP
!-------------------------------------------------------------------------
!BOC
   if (last_action .ne. WAITING) then
      FATAL 'Last action was not WAITING - not ready for sending (2D)'
      call MPI_ABORT(active_comm,-1,ierr)
   end if
   select case (comm_method)
      case(ONE_PROCESS)
         f1(imin-1, : )  = f2(imin, :  )
         f1(imax+1, : )  = f2(imax, :  )
         f1( :, jmin-1 ) = f2( :, jmin )
         f1( :, jmax+1 ) = f2( :, jmax )
      case(ONED_SENDRECV)
         if(com_direction .eq. RIGHT_LEFT) then
#ifdef DEBUG
STDERR 'ONED_SENDRECV - y_line'
#endif
            call MPI_SENDRECV(f1(imin,jmin),   1, y_line, left , tag, &
                              f2(imax+1,jmin), 1, y_line, right, tag, &
                              active_comm, status, ierr)

            call MPI_SENDRECV(f1(imax,jmin),   1, y_line, right, tag, &
                              f2(imin-1,jmin), 1, y_line, left , tag, &
                              active_comm, status, ierr)
         else
#ifdef DEBUG
STDERR 'ONED_SENDRECV - x_line'
#endif
            call MPI_SENDRECV(f1(imin,jmax),   1, x_line, up,   tag, &
                              f2(imin,jmin-1), 1, x_line, down, tag, &
                              active_comm, status, ierr)
            call MPI_SENDRECV(f1(imin,jmin),   1, x_line, down, tag, &
                              f2(imin,jmax+1), 1, x_line, up  , tag, &
                              active_comm, status, ierr)
         end if
      case(ONED_NONBLOCKING)
         if(com_direction .eq. RIGHT_LEFT) then
#ifdef DEBUG
STDERR 'ONED_NONBLOCKING - y_line'
#endif
            call MPI_IRECV(f2(imax+1,jmin), 1, y_line, right, tag, &
                              active_comm, req(2), ierr)
            call MPI_IRECV(f2(imin-1,jmin), 1, y_line, left,  tag, &
                              active_comm, req(1), ierr)
            call MPI_ISEND(f1(imin,jmin),   1, y_line, left,  tag, &
                              active_comm, req(4), ierr)
            call MPI_ISEND(f1(imax,jmin),   1, y_line, right, tag, &
                              active_comm, req(3), ierr)
         else
#ifdef DEBUG
STDERR 'ONED_NONBLOCKING - x_line'
#endif
            call MPI_IRECV(f2(imin,jmin-1), 1, x_line, down,  tag, &
                              active_comm, req(1), ierr)
            call MPI_IRECV(f2(imin,jmax+1), 1, x_line, up,    tag, &
                              active_comm, req(2), ierr)
            call MPI_ISEND(f1(imin,jmin),   1, x_line, down,  tag, &
                              active_comm, req(3), ierr)
            call MPI_ISEND(f1(imin,jmax),   1, x_line, up,    tag, &
                              active_comm, req(4), ierr)
         end if
      case(TWOD_SENDRECV)
#ifdef DEBUG
STDERR 'TWOD_SENDRECV'
#endif
         call MPI_SENDRECV(f1(imin,jmax),   1, x_line, up  , tag, &
                           f2(imin,jmin-1), 1, x_line, down, tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(imin,jmin),   1, x_line, down, tag, &
                           f2(imin,jmax+1), 1, x_line, up  , tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(imax,jmin),   1, y_line, right, tag, &
                           f2(imin-1,jmin), 1, y_line, left,  tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(imin,jmin),   1, y_line, left , tag, &
                           f2(imax+1,jmin), 1, y_line, right, tag, &
                           active_comm, status, ierr)
!        Corner points
         call MPI_SENDRECV(f1(imin,jmin),    1, MPI_REALTYPE, lleft, tag, &
                           f2(imax+1,jmax+1),1, MPI_REALTYPE, uright,tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(imax,jmin),    1, MPI_REALTYPE, lright,tag, &
                           f2(imin-1,jmax+1),1, MPI_REALTYPE, uleft, tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(imax,jmax),    1, MPI_REALTYPE, uright,tag, &
                           f2(imin-1,jmin-1),1, MPI_REALTYPE, lleft, tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(imin,jmax),    1, MPI_REALTYPE, uleft,tag, &
                           f2(imax+1,jmin-1),1, MPI_REALTYPE, lright,tag, &
                           active_comm, status, ierr)
      case(TWOD_NONBLOCKING)
#ifdef DEBUG
STDERR 'TWOD_NONBLOCKING'
#endif
!        Recieving x_lines
         call MPI_IRECV(f2(imin,jmin-1), 1, x_line, down,  tag, &
                           active_comm, req(1), ierr)
         call MPI_IRECV(f2(imin,jmax+1), 1, x_line, up,    tag, &
                           active_comm, req(2), ierr)

!        Recieving y_lines
         call MPI_IRECV(f2(imin-1,jmin), 1, y_line, left,  tag, &
                           active_comm, req(3), ierr)
         call MPI_IRECV(f2(imax+1,jmin), 1, y_line, right, tag, &
                           active_comm, req(4), ierr)

!        Recieving corner points
         call MPI_IRECV(f2(imin-1,jmin-1), 1, MPI_REALTYPE,lleft,tag, &
                           active_comm, req(5), ierr)
         call MPI_IRECV(f2(imax+1,jmin-1), 1, MPI_REALTYPE,lright,tag, &
                           active_comm, req(6), ierr)
         call MPI_IRECV(f2(imax+1,jmax+1), 1, MPI_REALTYPE,uright,tag, &
                           active_comm, req(7), ierr)
         call MPI_IRECV(f2(imin-1,jmax+1), 1, MPI_REALTYPE,uleft,tag, &
                           active_comm, req(8), ierr)

!        Sending x_lines
         call MPI_ISEND(f1(imin,jmin),   1, x_line, down,  tag, &
                           active_comm, req(9), ierr)
         call MPI_ISEND(f1(imin,jmax),   1, x_line, up,    tag, &
                           active_comm, req(10), ierr)

!        Sending y_lines
         call MPI_ISEND(f1(imin,jmin),   1, y_line, left, tag, &
                           active_comm, req(11), ierr)
         call MPI_ISEND(f1(imax,jmin),   1, y_line, right,  tag, &
                           active_comm, req(12), ierr)

!        Sending corner points
         call MPI_ISEND(f1(imax,jmax), 1, MPI_REALTYPE, uright,tag, &
                           active_comm, req(13), ierr)

         call MPI_ISEND(f1(imin,jmax), 1, MPI_REALTYPE, uleft,tag, &
                           active_comm, req(14), ierr)

         call MPI_ISEND(f1(imin,jmin), 1, MPI_REALTYPE, lleft,tag, &
                           active_comm, req(15), ierr)

         call MPI_ISEND(f1(imax,jmin), 1, MPI_REALTYPE, lright,tag, &
                           active_comm, req(16), ierr)
      case default
         FATAL 'A non valid communication method has been chosen'
         stop 'update_2d_halo_mpi'
   end select
   last_action = SENDING
   return
   end subroutine update_2d_halo_mpi
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: update_3d_halo_mpi - updates the halo zones for 3D fields.
!
! !INTERFACE:
   SUBROUTINE update_3d_halo_mpi(f1,f2,iimin,jjmin,iimax,jjmax,kmax,tag)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Updates the halo zones for 3D fields.
!
! !INPUT PARAMTERS:
   integer, INTENT(IN)		:: iimin,jjmin,iimax,jjmax,kmax
   integer, INTENT(IN)		:: tag
!
! !INPUT/OUTPUT PARAMTERS:
   REALTYPE, INTENT(INOUT), DIMENSION(I3DFIELD)	:: f1,f2
!
! !OUTPUT PARAMTERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
!
!EOP
!-------------------------------------------------------------------------
!BOC
   if (last_action .ne. WAITING) then
      FATAL 'Last action was not WAITING - not ready for sending (3D)'
      call MPI_ABORT(active_comm,-1,ierr)
   end if
   select case (comm_method)
      case(ONE_PROCESS)
         f1(iimin-1, :, : )  = f2(iimin, :, :  )
         f1(iimax+1, :, : )  = f2(iimax, :, :  )
         f1( :, jjmin-1, : ) = f2( :, jjmin, : )
         f1( :, jjmax+1, : ) = f2( :, jjmax, : )
      case(ONED_SENDRECV)
         if(com_direction .eq. RIGHT_LEFT) then
#ifdef DEBUG
STDERR 'ONED_SENDRECV - yz_slice'
#endif
            call MPI_SENDRECV(f1(iimin,jjmin,0),   1, yz_slice, left , tag, &
                              f2(iimax+1,jjmin,0), 1, yz_slice, right, tag, &
                              active_comm, status, ierr)
            call MPI_SENDRECV(f1(iimax,jjmin,0),   1, yz_slice, right, tag, &
                              f2(iimin-1,jjmin,0), 1, yz_slice, left , tag, &
                              active_comm, status, ierr)
         else
#ifdef DEBUG
STDERR 'ONED_SENDRECV - xz_slice'
#endif
            call MPI_SENDRECV(f1(iimin,jjmin,0),   1, xz_slice, down, tag, &
                              f2(iimin,jjmax+1,0), 1, xz_slice, up  , tag, &
                              active_comm, status, ierr)
            call MPI_SENDRECV(f1(iimin,jjmax,0),   1, xz_slice, up,   tag, &
                              f2(iimin,jjmin-1,0), 1, xz_slice, down, tag, &
                              active_comm, status, ierr)
         end if
      case(ONED_NONBLOCKING)
         if(com_direction .eq. RIGHT_LEFT) then
#ifdef DEBUG
STDERR 'ONED_NONBLOCKING - yz_slice'
#endif
            call MPI_IRECV(f2(iimax+1,jjmin,0), 1, yz_slice, right, tag, &
                              active_comm, req(2), ierr)
            call MPI_IRECV(f2(iimin-1,jjmin,0), 1, yz_slice, left,  tag, &
                              active_comm, req(1), ierr)
            call MPI_ISEND(f1(iimin,jjmin,0),   1, yz_slice, left,  tag, &
                              active_comm, req(4), ierr)
            call MPI_ISEND(f1(iimax,jjmin,0),   1, yz_slice, right, tag, &
                              active_comm, req(3), ierr)
         else
#ifdef DEBUG
STDERR 'ONED_NONBLOCKING - xz_slice'
#endif
            call MPI_IRECV(f2(iimin,jjmin-1,0), 1, xz_slice, down,  tag, &
                              active_comm, req(1), ierr)
            call MPI_IRECV(f2(iimin,jjmax+1,0), 1, xz_slice, up,    tag, &
                              active_comm, req(2), ierr)
            call MPI_ISEND(f1(iimin,jjmin,0),   1, xz_slice, down,  tag, &
                              active_comm, req(3), ierr)
            call MPI_ISEND(f1(iimin,jjmax,0),   1, xz_slice, up,    tag, &
                              active_comm, req(4), ierr)
         end if
      case(TWOD_SENDRECV)
#ifdef DEBUG
STDERR 'TWOD_SENDRECV'
#endif
         call MPI_SENDRECV(f1(iimin,jjmax,0),   1, xz_slice, up  , tag, &
                           f2(iimin,jjmin-1,0), 1, xz_slice, down, tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(iimin,jjmin,0),   1, xz_slice, down, tag, &
                           f2(iimin,jjmax+1,0), 1, xz_slice, up  , tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(iimax,jjmin,0),   1, yz_slice, right, tag, &
                           f2(iimin-1,jjmin,0), 1, yz_slice, left,  tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(iimin,jjmin,0),   1, yz_slice, left , tag, &
                           f2(iimax+1,jjmin,0), 1, yz_slice, right, tag, &
                           active_comm, status, ierr)
!        Corner points
         call MPI_SENDRECV(f1(iimin,jjmin,0),    1, z_column, lleft, tag, &
                           f2(iimax+1,jjmax+1,0),1, z_column, uright,tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(iimax,jjmin,0),    1, z_column, lright,tag, &
                           f2(iimin-1,jjmax+1,0),1, z_column, uleft, tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(iimax,jjmax,0),    1, z_column, uright,tag, &
                           f2(iimin-1,jjmin-1,0),1, z_column, lleft, tag, &
                           active_comm, status, ierr)
         call MPI_SENDRECV(f1(iimin,jjmax,0),    1, z_column, uleft,tag, &
                           f2(iimax+1,jjmin-1,0),1, z_column, lright,tag, &
                           active_comm, status, ierr)
      case(TWOD_NONBLOCKING)
#ifdef DEBUG
STDERR 'TWOD_NONBLOCKING'
#endif
!        Recieving xz_slice
         call MPI_IRECV(f2(iimin,jjmin-1,0), 1, xz_slice, down,  tag, &
                           active_comm, req(1), ierr)
         call MPI_IRECV(f2(iimin,jjmax+1,0), 1, xz_slice, up,    tag, &
                           active_comm, req(2), ierr)

!        Recieving yz_slice
         call MPI_IRECV(f2(iimin-1,jjmin,0), 1, yz_slice, left,  tag, &
                           active_comm, req(3), ierr)
         call MPI_IRECV(f2(iimax+1,jjmin,0), 1, yz_slice, right, tag, &
                           active_comm, req(4), ierr)

!        Recieving corner columns
         call MPI_IRECV(f2(iimin-1,jjmin-1,0), 1, z_column,lleft,tag, &
                           active_comm, req(5), ierr)
         call MPI_IRECV(f2(iimax+1,jjmin-1,0), 1, z_column,lright,tag, &
                           active_comm, req(6), ierr)
         call MPI_IRECV(f2(iimax+1,jjmax+1,0), 1, z_column,uright,tag, &
                           active_comm, req(7), ierr)
         call MPI_IRECV(f2(iimin-1,jjmax+1,0), 1, z_column,uleft,tag, &
                           active_comm, req(8), ierr)

!        Sending xz_slice
         call MPI_ISEND(f1(iimin,jjmin,0),   1, xz_slice, down,  tag, &
                           active_comm, req(9), ierr)
         call MPI_ISEND(f1(iimin,jjmax,0),   1, xz_slice, up,    tag, &
                           active_comm, req(10), ierr)

!        Sending yz_slice
         call MPI_ISEND(f1(iimin,jjmin,0),   1, yz_slice, left, tag, &
                           active_comm, req(11), ierr)
         call MPI_ISEND(f1(iimax,jjmin,0),   1, yz_slice, right,  tag, &
                           active_comm, req(12), ierr)

!        Sending corner points
         call MPI_ISEND(f1(iimax,jjmax,0), 1, z_column, uright,tag, &
                           active_comm, req(13), ierr)

         call MPI_ISEND(f1(iimin,jjmax,0), 1, z_column, uleft,tag, &
                           active_comm, req(14), ierr)

         call MPI_ISEND(f1(iimin,jjmin,0), 1, z_column, lleft,tag, &
                           active_comm, req(15), ierr)

         call MPI_ISEND(f1(iimax,jjmin,0), 1, z_column, lright,tag, &
                           active_comm, req(16), ierr)
      case default
         FATAL 'A non valid communication method has been chosen'
         stop 'update_3d_halo_mpi'
   end select
   last_action = SENDING
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
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Call MPI\_WAITALL to wait for any un-finished communications. If SENDRECV
!  communications are used this call has no effect.
!
! !INPUT PARAMTERS:
   integer, INTENT(IN)		:: tag
!
! !INPUT/OUTPUT PARAMTERS:
!
! !OUTPUT PARAMTERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
!
!EOP
!-------------------------------------------------------------------------
!BOC
   if (last_action .ne. SENDING) then
      FATAL 'Last action was not sending - nothing to wait for'
      call MPI_ABORT(active_comm,-1,ierr)
   end if
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
   return
   end subroutine wait_halo_mpi
!EOC

!-----------------------------------------------------------------------

   end module halo_mpi

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
