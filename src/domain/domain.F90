!$Id: domain.F90,v 1.5 2003-04-07 14:34:42 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: domain - sets up the calculation domain.
!
! !INTERFACE:
   module domain
!
! !DESCRIPTION:
!
! !USES:
   use halo_zones, only	: update_2d_halo,wait_halo,H_TAG
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   integer 		:: format=NETCDF
   integer		:: grid_type=1,vert_cord=1
   REALTYPE, allocatable, dimension(:) :: ga
   REALTYPE		:: ddu=-_ONE_,ddl=-_ONE_,d_gamma=20.
   logical		:: gamma_surf=.true.
   integer		:: NWB,NNB,NEB,NSB,NOB
   REALTYPE    		:: latitude=0.
   REALTYPE    		:: Hland
   REALTYPE		:: min_depth,crit_depth
   logical		:: openbdy=.false.
   integer		:: calc_points
#ifdef STATIC
#include "static_domain.h"
#else
#include "dynamic_declarations_domain.h"
#endif
   integer		:: nsbv
   integer, parameter	:: INNER= 1
   integer		:: ioff=0,joff=0
   integer, dimension(:), allocatable	:: wi,wfj,wlj
   integer, dimension(:), allocatable	:: nj,nfi,nli
   integer, dimension(:), allocatable	:: ei,efj,elj
   integer, dimension(:), allocatable	:: sj,sfi,sli
   integer, allocatable			:: bdy_index(:),bdy_map(:,:)
!kbk
   REALTYPE		:: cori= _ZERO_
   REALTYPE, parameter	:: rearth=6370.9490e3
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: domain.F90,v $
!  Revision 1.5  2003-04-07 14:34:42  kbk
!  parallel support, proper spherical grid init. support
!
!  Revision 1.1.1.1  2002/05/02 14:01:11  gotm
!  recovering after CVS crash
!
!  Revision 1.19  2001/10/23 14:15:55  bbh
!  Moved ga from coordinates.F90 to domain.F90
!
!  Revision 1.18  2001/10/22 12:10:26  bbh
!  Partly support for SPHERICAL grid is coded
!
!  Revision 1.17  2001/09/26 10:01:41  bbh
!  lat and lon maps now read in ncdf_topo.F90
!
!  Revision 1.16  2001/09/24 07:49:32  bbh
!  Include .h files for memory declaration/allocation
!
!  Revision 1.15  2001/09/21 11:52:47  bbh
!  Minimum depth for specific areas through - set_min_depth()
!
!  Revision 1.14  2001/09/14 12:04:15  bbh
!  Added xc,yc to hold coordinates + cleaning
!
!  Revision 1.13  2001/09/04 08:00:14  bbh
!  Fill coru and corv arrays
!
!  Revision 1.12  2001/09/04 07:36:32  bbh
!  We need ioff and joff in parallel runs
!
!  Revision 1.11  2001/09/03 15:14:22  bbh
!  Bug with Coriolis removed
!
!  Revision 1.10  2001/09/01 17:10:25  bbh
!  Vertical coordinate definition now specified via namelist
!
!  Revision 1.9  2001/08/27 11:55:02  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.8  2001/08/01 08:19:57  bbh
!  Fields for CURVILINEAR - now done
!
!  Revision 1.7  2001/07/26 14:31:43  bbh
!  Manual merge
!
!  Revision 1.6  2001/07/26 14:20:02  bbh
!  Added grid_type, vert_cord, lonmap and latmap
!
!  Revision 1.5  2001/06/22 08:19:10  bbh
!  Compiler options such as USE_MASK and OLD_DRY deleted.
!  Open and passive boundary for z created.
!  Various inconsistencies removed.
!  wait_halo added.
!  Checked loop boundaries
!
!  Revision 1.4  2001/05/14 12:38:58  bbh
!  Set minimum detph to 10. meters if not COAST_TEST - to be fixed later.
!
!  Revision 1.3  2001/05/06 18:51:55  bbh
!  Towards proper implementation of specified 2D bdy.
!
!  Revision 1.2  2001/04/24 08:24:58  bbh
!  Use runtype instead of macro
!
!  Revision 1.1.1.1  2001/04/17 08:43:08  bbh
!  initial import into CVS
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_domain - initialise the computational domain
!
! !INTERFACE:
   subroutine init_domain(input_dir)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
   character(len=*)	:: input_dir
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Initialize the calculation domain - that is obtain information of the
!  extent of the domain - imin,imax,jmin and jmax.
!  This subroutine operates in two principally different ways.
!  In the first case it has been compiled with \#define STATIC which means that
!  most arrays are statically allocated. In this case the information of
!  the extenstion of the domain is obtained from a file called dimensions.h,
!  where imin,imax,jmin and jmax are specified as integer parameters.
!  If STATIC has not been defined during compilation all arrays are allocatable.!  In this case information of the extension of the domain is obtained from
!  the bathymetry file. After imin,imax,jmin and jmax are known arrays are
!  allocated.
!  Next step is to read in the batymetry. Then location of boundaries are read.
!  Hereafter calculation masks are setup - note that the masks are explicitely
!  set to 0 in the halo zones.
!  Finally the halo zones are filled. In the case of only one process the
!  information is simply copied where as in a parallel run the information is
!  communicated from the neighboring processes - this is all done in
!  update\_2d\_halo().
!
! !REVISION HISTORY:
!
!  See log for module
!
! !LOCAL VARIABLES:
   integer			:: rc
   integer			:: np,sz
   integer			:: i,j,n
   integer			:: kdum
   REALTYPE, parameter		:: pi=3.141592654
   REALTYPE, parameter		:: deg2rad=pi/180.
   REALTYPE, parameter		:: omega=2.*pi/86400.
   character(len=PATH_MAX)	:: bathymetry='topo.nc'
   character(len=PATH_MAX)	:: bdyinfofile='bdyinfo.dat'
   character(len=PATH_MAX)	:: min_depth_file='minimum_depth.dat'
   character(len=PATH_MAX)	:: bathymetry_adjust_file='bathymetry.adjust'
   character(len=PATH_MAX)	:: mask_adjust_file='mask.adjust'
   integer			:: il=-1,ih=-1,jl=-1,jh=-1
   namelist /domain/		grid_type,vert_cord,ddu,ddl,	&
   				d_gamma,gamma_surf,	&
                                bathymetry,openbdy,bdyinfofile,latitude, &
                                crit_depth,min_depth,kdum,il,ih,jl,jh
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_domain()'
#endif

   LEVEL1 'init_domain'

!  Read domain specific things from the namelist.
   read(NAMLST,domain)
!kbk   REWIND(NAMLST)

   select case (grid_type)
      case(1)
         LEVEL2 'Using cartesian grid'
      case(2)
         LEVEL2 'Using spherical grid'
      case(3)
         LEVEL2 'Using curvi-linear grid'
      case default
         FATAL 'A non valid grid type has been chosen'
         stop 'init_domain'
   end select

   select case (vert_cord)
      case(1)
         LEVEL2 'Using sigma coordinates'
      case(2)
         LEVEL2 'Using z-level coordinates'
      case(3)
         LEVEL2 'Using general vertical coordinates'
      case default
         FATAL 'A non valid vertical coordinate system has been chosen'
         stop 'init_domain'
   end select

   select case (format)
      case(NETCDF)
#ifdef STATIC
!         call get_dimensions(trim(input_dir) // bathymetry,rc)
         call get_dimensions('topo.nc',rc)
#else
call get_dimensions(trim(input_dir) // bathymetry,iextr,jextr,rc)
         kmax=kdum
#endif
         call part_domain()
         il=imin ; ih=imax ; jl=jmin ; jh=jmax
#ifndef STATIC
#include "dynamic_allocations_domain.h"
#endif
         H = -10.
         HU = -10.
         HV = -10.

lonc = -1000.
latc = -1000.

         call get_bathymetry(H,Hland,iextr,jextr,ioff,joff,imin,imax,jmin,jmax,rc)
         call update_2d_halo(H,H,az,imin,jmin,imax,jmax,H_TAG)
         call wait_halo(H_TAG)

         call update_2d_halo(lonc,lonc,az,imin,jmin,imax,jmax,H_TAG)
         call wait_halo(H_TAG)

         call update_2d_halo(latc,latc,az,imin,jmin,imax,jmax,H_TAG)
         call wait_halo(H_TAG)
      case default
         FATAL 'A non valid input format has been chosen'
         stop 'init_domain'
   end select

#if 0
#ifdef DK_03NM_TEST
where (H .lt. 2.) 
H = -10.
end where
#endif
#endif

   select case (grid_type)
      case(1)
#if ! ( defined SPHERICAL || defined CURVILINEAR)
         xc(1) = (ioff+0.5)*dx
         do i=imin+1,imax
            xc(i) = xc(i-1)+dx
         end do
         yc(1) = (joff+0.5)*dy
         do j=jmin+1,jmax
            yc(j) = yc(j-1)+dy
         end do
         ard1=_ONE_/(dx*dy)
!        The Coriolis parameter
         cori = 2.*omega*sin(deg2rad*latitude)
         coru = cori
         corv = cori
#endif
      case(2)
#ifdef SPHERICAL
          lonu=lonx
          lonv=lonc
          latv=latx
          latu=latc
          do j=jmin,jmax
             do i=imin,imax
                dxc(i,j)=deg2rad*(lonu(i,j)-lonu(i-1,j))*rearth &
		         *cos(deg2rad*latc(i,j))
             end do
          end do
          do j=jmin,jmax
             do i=imin-1,imax
                dxu(i,j)=deg2rad*(lonc(i+1,j)-lonc(i,j))*rearth &
		         *cos(deg2rad*latc(i,j))
             end do
          end do
          do j=jmin-1,jmax
             do i=imin,imax
                dxv(i,j)=deg2rad*(lonx(i,j)-lonx(i-1,j))*rearth &
		         *cos(deg2rad*latv(i,j))
             end do
          end do
          do j=jmin-1,jmax
             do i=imin-1,imax
                dxx(i,j)=deg2rad*(lonv(i+1,j)-lonv(i,j))*rearth &
		         *cos(deg2rad*latx(i,j))
             end do
          end do
          do j=jmin,jmax
             do i=imin,imax
                dyc(i,j)=deg2rad*(latv(i,j)-latv(i,j-1))*rearth
             end do
          end do
          do i=imin-1,imax
             do j=jmin,jmax
                dyu(i,j)=deg2rad*(latx(i,j)-latx(i,j-1))*rearth
             end do
          end do
          do j=jmin-1,jmax
             do i=imin,imax
                dyv(i,j)=deg2rad*(latc(i,j+1)-latc(i,j))*rearth
             end do
          end do
          do j=jmin-1,jmax
             do i=imin-1,imax
                dyx(i,j)=deg2rad*(latu(i,j+1)-latu(i,j))*rearth
             end do
          end do
          do j=jmin,jmax
             do i=imin,imax
                arcd1(i,j)=_ONE_/(dxc(i,j)*dyc(i,j))
                arud1(i,j)=_ONE_/(dxu(i,j)*dyu(i,j))
                arvd1(i,j)=_ONE_/(dxv(i,j)*dyv(i,j))
             end do
          end do
          do j=jmin,jmax
             do i=imin,imax
                coru(i,j)=2.*omega*sin(deg2rad*latu(i,j))
                corv(i,j)=2.*omega*sin(deg2rad*latv(i,j))
             end do
          end do
#endif
      case(3)
#ifdef CURVILINEAR
         do i=imin-1,imax
            do j=jmin,jmax
	       xu(i,j)=0.5*(xx(i,j)+xx(i,j-1))
	       yu(i,j)=0.5*(yx(i,j)+yx(i,j-1))
            end do
	 end do 	
         do i=imin,imax
            do j=jmin-1,jmax
	       xv(i,j)=0.5*(xx(i,j)+xx(i-1,j))
	       yv(i,j)=0.5*(yx(i,j)+yx(i-1,j))
            end do
	 end do 	
         do i=imin,imax
            do j=jmin,jmax
	       if (H(i,j).gt.Hland) then 
	          dxc(i,j)=sqrt((xu(i,j)-xu(i-1,j))**2+(yu(i,j)-yu(i-1,j))**2)
	          dyc(i,j)=sqrt((xv(i,j)-xv(i,j-1))**2+(yv(i,j)-yv(i,j-1))**2)
                  arcd1(i,j)=_ONE_/(dxc(i,j)*dyc(i,j))
	       end if 
            end do
	 end do 	
         do j=jmin,jmax
            do i=imin,imax-1
	       if ((H(i+1,j).gt.Hland).and.(H(i,j).gt.Hland)) then 
	          dxu(i,j)=sqrt((xc(i+1,j)-xc(i,j))**2+(yc(i+1,j)-yc(i,j))**2)
	          dyu(i,j)=sqrt((xx(i,j)-xx(i,j-1))**2+(yx(i,j)-yx(i,j-1))**2) 
                  arud1(i,j)=_ONE_/(dxu(i,j)*dyu(i,j))
	       end if
            end do
	    dxu(imin-1,j)=dxu(imin,j)
	    dxu(imax,j)=dxu(imax-1,j)
         end do 	
         do i=imin,imax
            do j=jmin,jmax-1
	       if ((H(i,j+1).gt.Hland).and.(H(i,j).gt.Hland)) then 
	          dyv(i,j)=sqrt((xc(i,j+1)-xc(i,j))**2+(yc(i,j+1)-yc(i,j))**2)
	          dxv(i,j)=sqrt((xx(i,j)-xx(i-1,j))**2+(yx(i,j)-yx(i-1,j))**2) 
                  arvd1(i,j)=_ONE_/(dxv(i,j)*dyv(i,j))
	       end if
            end do
	    dyv(i,jmin-1)=dyv(i,jmin)
	    dyv(i,jmax)=dyv(i,jmax-1)
         end do 	
	 do j=jmin-1,jmax
	    do i=imin,imax-1
	       if (((H(i,j+1).gt.Hland).and.(H(i,j).gt.Hland)).or.         &
	           ((H(i+1,j).gt.Hland).and.(H(i+1,j+1).gt.Hland))) then 
	          dxx(i,j)=sqrt((xv(i+1,j)-xv(i,j))**2+(yv(i+1,j)-yv(i,j))**2)
	       end if 	  
            end do
	    dxx(imin-1,j)=dxx(imin,j)
	    dxx(imax,j)=dxx(imax-1,j)
	 end do 	
	 do i=imin-1,imax
	    do j=jmin,jmax-1
	       if (((H(i,j+1).gt.Hland).and.(H(i,j).gt.Hland)).or.         &
	           ((H(i+1,j).gt.Hland).and.(H(i+1,j+1).gt.Hland))) then 
	          dyx(i,j)=sqrt((xu(i,j+1)-xu(i,j))**2+(yu(i,j+1)-yu(i,j))**2)
	       end if 	  
            end do
	    dyx(i,jmin-1)=dyx(i,jmin)
	    dyx(i,jmax)=dyx(i,jmax-1)
	 end do 	
	 coru = 2.*omega*sin(deg2rad*latu)
	 corv = 2.*omega*sin(deg2rad*latv)
#endif
      case default
         FATAL 'A non valid grid type has been chosen'
         stop 'init_domain'
   end select

!  Do we want to set a minimum depth for certain regions
   call set_min_depth(trim(input_dir) // min_depth_file)

!  Do we want to do adjust the bathymetry
   call adjust_bathymetry(trim(input_dir) // bathymetry_adjust_file)

!  Reads boundary location information
   if (openbdy) then
      call bdy_spec(trim(input_dir) // bdyinfofile)
      call print_bdy('Global Boundary Information')
call have_bdy()
call print_bdy('Local Boundary Information')
   end if

!  Define calculation masks
   az = 0
   where (H .gt. Hland+SMALL)
      az=1
   end where

#define BOUNDARY_POINT 2
!  western boundary - at present elev only
   do n=1,NWB
      az(wi(n),wfj(n):wlj(n)) = BOUNDARY_POINT
      if(wfj(n) .eq. jmin) az(wi(n),jmin-1) = az(wi(n),jmin)
      if(wlj(n) .eq. jmax) az(wi(n),jmax+1) = az(wi(n),jmax)
   end do
!  northern boundary - at present elev only
   do n=1,NNB
      az(nfi(n):nli(n),nj(n)) = BOUNDARY_POINT
      if(nfi(n) .eq. imin) az(imin-1,nj(n)) = az(imin,nj(n))
      if(nli(n) .eq. imax) az(imax+1,nj(n)) = az(imax,nj(n))
   end do
!  easter boundary - at present elev only
   do n=1,NEB
      az(ei(n),efj(n):elj(n)) = BOUNDARY_POINT
      if(efj(n) .eq. jmin) az(ei(n),jmin-1) = az(ei(n),jmin)
      if(elj(n) .eq. jmax) az(ei(n),jmax+1) = az(ei(n),jmax)
   end do
!  southern boundary - at present elev only
   do n=1,NSB
      az(sfi(n):sli(n),sj(n)) = BOUNDARY_POINT
      if(sfi(n) .eq. imin) az(imin-1,sj(n)) = az(imin,sj(n))
      if(sli(n) .eq. imax) az(imax+1,sj(n)) = az(imax,sj(n))
   end do
#undef BOUNDARY_POINT

!  Do we want to do adjust the mask
   call adjust_mask(trim(input_dir) // mask_adjust_file)

   au=0
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1 .and. az(i+1,j) .eq. 1) then
            au(i,j)=1
         end if
         if ((az(i,j) .eq. 1 .and. az(i+1,j) .eq. 2).or.    &
             (az(i,j) .eq. 2 .and. az(i+1,j) .eq. 1)) then
            au(i,j)=2
         end if
         if (az(i,j) .eq. 2 .and. az(i+1,j) .eq. 2) then
            au(i,j)=3
         end if
      end do
   end do

   av=0
   do j=jmin,jmax
      do i=imin,imax
         if (az(i,j) .eq. 1 .and. az(i,j+1) .eq. 1) then
            av(i,j)=1
         end if
         if ((az(i,j) .eq. 1 .and. az(i,j+1) .eq. 2).or.    &
             (az(i,j) .eq. 2 .and. az(i,j+1) .eq. 1)) then
            av(i,j)=2
         end if
         if (az(i,j) .eq. 2 .and. az(i,j+1) .eq. 2) then
	    av(i,j)=3
	 end if
      end do
   end do

   ax=0
   do j=jmin,jmax
      do i=imin,imax
         if (az(i  ,j) .eq. 1 .and. az(i  ,j+1) .eq. 1 .and.    &
             az(i+1,j) .eq. 1 .and. az(i+1,j+1) .eq. 1) then
	        ax(i,j)=1
	 end if 	
      end do
   end do

#ifdef DEBUG
   STDERR 'az'
   call print_mask(az)
   STDERR 'au'
   call print_mask(au)
   STDERR 'av'
   call print_mask(av)
#endif

   np = count(az(1:imax,1:jmax) .gt. 0)
   sz = (imax-imin+1)*(jmax-jmin+1)
   LEVEL2 'Dimensions: ',imin,':',imax,',',jmin,':',jmax,',',0,':',kmax
   LEVEL2 '# waterpoints = ',np,' of ',sz

   calc_points = np

#ifdef DEBUG
   write(debug,*) 'Leaving init_domain()'
   write(debug,*)
#endif
   return
   end subroutine init_domain
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_par_setup - reads domain partition
!
! !INTERFACE:
   subroutine read_par_setup(myid)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)			::	myid
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Reads the partitioning of the domain in parallel run
!
! !REVISION HISTORY:
!
!  22Apr99   Karsten Bolding & Hans Burchard  Initial code.
!
! !LOCAL VARIABLES:
   integer				:: id
!
!EOP
!------------------------------------------------------------------------
!BOC

#ifndef STATIC
!   open(PARSETUP,file=input_dir // 'par_setup')
   read(PARSETUP,*)

100  read(PARSETUP,*,ERR=110,END=111) id,imin,imax,jmin,jmax


   if(id .eq. myid ) then 
      close(PARSETUP)
      ioff=imin-1 ; joff=jmin-1
      imax=imax-imin+1 ; imin=1
      jmax=jmax-jmin+1 ; jmin=1
      LEVEL2 'From read_par_setup ',id,ioff,imin,imax,joff,jmin,jmax
      return
   end if

   goto 100

110   FATAL 'reading domain partition information.'
      stop
111   FATAL 'End of file reached - (read_par_setup).'
      stop

   stop 'read_par_setup'
#endif
   return
   end subroutine read_par_setup
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_min_depth -
!
! !INTERFACE:
   subroutine set_min_depth(fn)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)	:: fn
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Read mask adjustments from file.
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
   integer	:: unit = 25 ! kbk
   integer	:: i,j,k,n
   integer	:: il,jl,ih,jh
   integer	:: i1,j1
   REALTYPE	:: dmin
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Should read in to a buffer at some time to allow for #
   open(unit,file=fn,action='read',status='old',err=90)
   read(unit,*,end=91,err=92) n
   if(n .ge. 1) then
      LEVEL2 'Setting minimum depths according to ',trim(fn)
   end if
   do k=1,n
      read(unit,*,end=91,err=92) il,jl,ih,jh,dmin
      LEVEL3 'setting min depth in ',il,jl,ih,jh,' to ',dmin
      do j=jl,jh
         do i=il,ih
            if(imin+ioff .le. i .and. i .le. imax+ioff .and. &
               jmin+joff .le. j .and. j .le. jmax+joff ) then
               i1 = i-ioff
               j1 = j-joff
               if(H(i1,j1) .gt. -9. .and. H(i1,j1) .lt. dmin) then
                  H(i1,j1) = dmin
               end if
            end if
         end do
      end do
   end do
   close(unit)

   return

90 LEVEL2 'could not open ',trim(fn),' no minimum depth adjustments done'
   return
91 FATAL 'eof: ',trim(fn)
   stop 'set_min_depth'
92 FATAL 'error reading: ',trim(fn)
   stop 'set_min_depth'
   end subroutine set_min_depth
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: adjust_bathymetry - read mask adjustments from file.
!
! !INTERFACE:
   subroutine adjust_bathymetry(fn)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)	:: fn
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Read bathymetry adjustments from file.
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
   integer	:: unit = 25 ! kbk
   integer	:: i,j,k,n
   integer	:: il,jl,ih,jh
   integer	:: i1,j1
   REALTYPE	:: x
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Should read in to a buffer at some time to allow for #
   open(unit,file=fn,action='read',status='old',err=90)
   read(unit,*,end=91,err=92) n
   if(n .gt. 1) then
      LEVEL2 'Adjusting bathymetry according to ',trim(fn)
   end if
   do k=1,n
      read(unit,*,end=91,err=92) il,jl,ih,jh,x
      LEVEL3 'setting bathymetry in ',il,jl,ih,jh,' to ',x
      do j=jl,jh
         do i=il,ih
            if(imin+ioff .le. i .and. i .le. imax+ioff .and. &
               jmin+joff .le. j .and. j .le. jmax+joff ) then
               i1 = i-ioff
               j1 = j-joff
	       H(i1,j1) = x
            end if
         end do
      end do
   end do
   close(unit)

   return

90 LEVEL2 'could not open ',trim(fn),' no bathymetry adjustments done'
   return
91 FATAL 'eof: ',trim(fn)
   stop 'adjust_bathymetry'
92 FATAL 'error reading: ',trim(fn)
   stop 'adjust_bathymetry'
   end subroutine adjust_bathymetry
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: adjust_mask - read mask adjustments from file.
!
! !INTERFACE:
   subroutine adjust_mask(fn)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)	:: fn
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Read mask adjustments from file.
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
   integer	:: unit = 25 ! kbk
   integer	:: i,j,k,n
   integer	:: il,jl,ih,jh
   integer	:: i1,j1
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Should read in to a buffer at some time to allow for #
   open(unit,file=fn,action='read',status='old',err=90)
   read(unit,*,end=91,err=92) n
   if(n .gt. 1) then
      LEVEL2 'Adjusting mask according to ',trim(fn)
   end if
   do k=1,n
      read(unit,*,end=91,err=92) il,jl,ih,jh
      LEVEL3 'masking area ',il,jl,ih,jh
      do j=jl,jh
         do i=il,ih
            if(imin+ioff .le. i .and. i .le. imax+ioff .and. &
               jmin+joff .le. j .and. j .le. jmax+joff ) then
               i1 = i-ioff
               j1 = j-joff
               az(i1,j1) = 0
            end if
         end do
      end do
   end do
   close(unit)

   return

90 LEVEL2 'could not open ',trim(fn),' no mask adjustments done'
   return
91 FATAL 'eof: ',trim(fn)
   stop 'adjust_mask'
92 FATAL 'error reading: ',trim(fn)
   stop 'adjust_mask'
   end subroutine adjust_mask
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: print_mask - prints a mask in readable format
!
! !INTERFACE:
   subroutine print_mask(mask)
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in), dimension(E2DFIELD) 	:: mask
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Prints a integer mask in a human readable form.
!
! !REVISION HISTORY:
!
!  22Apr99   Karsten Bolding & Hans Burchard  Initial code.
!
! !LOCAL VARIABLES:
   integer	:: i,j
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
#endif

   do j=jmax,jmin,-1
!      write(0,'(5000(i1,x))') (mask(i,j), i=imin,imax)
      write(0,'(5000(i1))') (mask(i,j), i=imin,imax)
   end do

   return
   end subroutine print_mask
!EOC

!-----------------------------------------------------------------------

   end module domain

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
