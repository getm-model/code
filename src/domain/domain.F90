!$Id: domain.F90,v 1.3 2003-03-17 15:00:20 gotm Exp $
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
   use commhalo, only	: myid,nprocs,comm_hd,update_2d_halo,H_TAG
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   integer 		:: format=NETCDF
   integer		:: grid_type=1,vert_cord=1
   REALTYPE, allocatable, dimension(:) :: ga
   REALTYPE		:: ddu=-_ONE_,ddl=-_ONE_,d_gamma=20.
   logical		:: gamma_surf=.true.
   integer		:: NWB,NNB,NEB,NSB
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
!kbk   integer, parameter			:: OPENBDY=2
!kbk  integer, parameter			:: HALO=3
   integer		:: ioff=0,joff=0
   integer, dimension(:), allocatable	:: wi,wfj,wlj
   integer, dimension(:), allocatable	:: nj,nfi,nli
   integer, dimension(:), allocatable	:: ei,efj,elj
   integer, dimension(:), allocatable	:: sj,sfi,sli
!kbk
   REALTYPE		:: cori=0.
   REALTYPE, parameter	:: rearth=6370.9490e3
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: domain.F90,v $
!  Revision 1.3  2003-03-17 15:00:20  gotm
!  setting lonmap and latmap
!
!  Revision 1.2  2002/05/29 13:37:50  gotm
!  New naming of .h files
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
   subroutine init_domain()
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
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
   write(debug,*) 'init_domain(): id =  ',myid
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
         call get_dimensions(bathymetry,rc)
#else
         call get_dimensions(bathymetry,iextr,jextr,rc)
#endif

         if ( myid .lt. 0 .or. nprocs .eq. 1) then
#ifndef STATIC
!kbk          imin = 1 ; imax = iextr ; jmin = 1 ; jmax = jextr;
            if(il .eq. -1 .or. ih .eq. -1 .or. jl .eq. -1 .or. jh .eq. -1) then
               imin = 1 ; imax = iextr ; jmin = 1 ; jmax = jextr;
               il = imin ; il = imax ; jl = jmin ; jh = jmax
            else
               imin = 1 ; imax = ih-il+1 ; jmin = 1 ; jmax = jh-jl+1;
            end if
#endif
         else
            call read_par_setup(myid)
         end if
         il = imin ; ih = imax ; jl = jmin ; jh = jmax
#ifndef STATIC
         iimin = imin ; iimax = imax ; jjmin = jmin ; jjmax = jmax; kmax = kdum
#include "dynamic_allocations_domain.h"
#endif
         call get_bathymetry(H,Hland,il,ih,jl,jh,rc)
      case default
         FATAL 'A non valid input format has been chosen'
         stop 'init_domain'
   end select

#ifdef DK_03NM_TEST
where (H .lt. 2.) 
H = -10.
end where
#endif

! Define calculation masks
   az = 0
   az(imin-HALO:imin, : ) = 0
   az(imax:imax+HALO, : ) = 0
   az( :, jmin-HALO:jmin) = 0
   az( :, jmax:jmax+HALO) = 0

   az(imin-HALO:imin-1, : ) = 0
   az(imax+1:imax+HALO, : ) = 0
   az( :, jmin-HALO:jmin-1) = 0
   az( :, jmax+1:jmax+HALO) = 0

   where (H .gt. Hland+SMALL)
      az=1
   end where

   select case (grid_type)
      case(1)
#if ! ( defined SPHERICAL || defined CURVILINEAR)
         xc(1) = 0.5*dx
         do i=imin+1,imax
            xc(i) = xc(i-1)+dx
         end do
         yc(1) = 0.5*dy
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
#if 0
#ifdef NS_NOMADS_TEST
       xx(0) = -5.083333-1./6.
       do i=imin,imax+1
          xx(i)=xx(0)+i*1./3.
       end do
       yx(0) = 48.8833-0.1
       do j=jmin,jmax+1
          yx(j)=yx(0)+j*0.2
       end do
#endif
#endif
#ifdef MED_15X15MINS_TEST
       do i=imin,imax
          lonmap(i,:)=-9.25+(i-1)*0.25
       end do

       do j=jmin,jmax
          latmap(:,j)=30.25+(j-1)*0.25
       end do
#endif
       do i=imin,imax+1
          xc(i)=0.5*(xx(i)+xx(i-1))
       end do
       xc(imin-1)=2.*xc(imin)-xc(imin+1)
!KBK       xc(imax+1)=xc(imax) + 0.5*(xx(imax+1)-xx(imax))
       xu=xx
       xv=xc
       do j=jmin,jmax+1
          yc(j)=0.5*(yx(j)+yx(j-1))
       end do
       yc(jmin-1)=2.*yc(jmin)-yc(jmin+1)
!KBK       yc(jmax+1)=yc(jmax) + 0.5*(yx(jmax+1)-yx(jmax))
       yv=yx
       yu=yc
       do j=jmin,jmax
          do i=imin,imax
             dxc(i,j)=deg2rad*(xu(i)-xu(i-1))*rearth*cos(deg2rad*yc(j))
          end do
       end do
       do j=jmin,jmax
          do i=imin-1,imax
             dxu(i,j)=deg2rad*(xc(i+1)-xc(i))*rearth*cos(deg2rad*yc(j))
          end do
       end do
       do j=jmin-1,jmax
          do i=imin,imax
             dxv(i,j)=deg2rad*(xx(i)-xx(i-1))*rearth*cos(deg2rad*yv(j))
          end do
       end do
       do j=jmin-1,jmax
          do i=imin-1,imax
             dxx(i,j)=deg2rad*(xv(i+1)-xv(i))*rearth*cos(deg2rad*yx(j))
          end do
       end do
       do j=jmin,jmax
          do i=imin,imax
             dyc(i,j)=deg2rad*(yv(j)-yv(j-1))*rearth
          end do
       end do
       do i=imin-1,imax
          do j=jmin,jmax
             dyu(i,j)=deg2rad*(yx(j)-yx(j-1))*rearth
          end do
       end do
       do j=jmin-1,jmax
          do i=imin,imax
             dyv(i,j)=deg2rad*(yc(j+1)-yc(j))*rearth
          end do
       end do
       do j=jmin-1,jmax
          do i=imin-1,imax
             dyx(i,j)=deg2rad*(yu(j+1)-yu(j))*rearth
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
             coru(i,j)=2.*omega*sin(deg2rad*yu(j))
             corv(i,j)=2.*omega*sin(deg2rad*yv(j))
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
lonmap = lonc
latmap = latc
#endif
      case default
         FATAL 'A non valid grid type has been chosen'
         stop 'init_domain'
   end select

!  Do we want to set a minimum depth for certain regions
   call set_min_depth(min_depth_file)

!  Do we want to do adjust the bathymetry
   call adjust_bathymetry(bathymetry_adjust_file)

!  Reads boundary location information
   if (openbdy) then
      call bdy_spec(bdyinfofile)
      call print_bdy('Global Boundary Information')
   end if

! Define calculation masks
   az = 0
   az(imin-HALO:imin, : ) = 0
   az(imax:imax+HALO, : ) = 0
   az( :, jmin-HALO:jmin) = 0
   az( :, jmax:jmax+HALO) = 0

   az(imin-HALO:imin-1, : ) = 0
   az(imax+1:imax+HALO, : ) = 0
   az( :, jmin-HALO:jmin-1) = 0
   az( :, jmax+1:jmax+HALO) = 0

   where (H .gt. Hland+SMALL)
      az=1
   end where

#define BOUNDARY_POINT 2
!  western boundary - at present elev only
   do n=1,NWB
#if 0
      i=wi(n)
      do j=wfj(n),wlj(n)
         az(i,j) = BOUNDARY_POINT
      end do
#else
      az(wi(n),wfj(n):wlj(n)) = BOUNDARY_POINT
#endif
   end do
!  northern boundary - at present elev only
   do n=1,NNB
#if 0
      j=nj(n)
      do i=nfi(n),nli(n)
         az(i,j) = BOUNDARY_POINT
      end do
#else
      az(nfi(n):nli(n),nj(n)) = BOUNDARY_POINT
#endif
   end do
!  easter boundary - at present elev only
   do n=1,NEB
#if 0
      i=ei(n)
      do j=efj(n),elj(n)
         az(i,j) = BOUNDARY_POINT
      end do
   end do
#else
      az(ei(n),efj(n):elj(n)) = BOUNDARY_POINT
#endif
   end do
!  southern boundary - at present elev only
   do n=1,NSB
#if 0
      j=sj(n)
      do i=sfi(n),sli(n)
         az(i,j) = BOUNDARY_POINT
      end do
#else
      az(sfi(n):sli(n),sj(n)) = BOUNDARY_POINT
#endif
   end do
#undef BOUNDARY_POINT

!  Do we want to do adjust the mask
   call adjust_mask(mask_adjust_file)

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
   np = count(az .gt. 0)
   sz = (imax-imin+1)*(jmax-jmin+1)
   LEVEL2 'Dimensions: ',imin,':',imax,',',jmin,':',jmax,',',0,':',kmax
   LEVEL2 '# waterpoints = ',np,' of ',sz

   calc_Points = np

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
   open(PARSETUP,file='par_setup')

   read(PARSETUP,*)

100  read(PARSETUP,*,ERR=110,END=111) id,imin,imax,jmin,jmax

   close(PARSETUP)

   LEVEL2 'From read_par_setup ',id,imin,imax,jmin,jmax

   if(id .eq. myid ) return

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
   REALTYPE	:: dmin
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Should read in to a buffer at some time to allow for #
   open(unit,file=fn,action='read',status='old',err=90)
   read(unit,*,end=91,err=92) n
   if(n .gt. 1) then
      LEVEL2 'Setting minimum depths according to ',trim(fn)
   end if
   do k=1,n
      read(unit,*,end=91,err=92) il,jl,ih,jh,dmin
      LEVEL3 'setting min depth in ',il,jl,ih,jh,' to ',dmin
      do j=jl,jh
         do i=il,ih
            if(H(i,j) .gt. -9. .and. H(i,j) .lt. dmin) then
               H(i,j) = dmin
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
	    H(i,j) = x
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
            az(i,j) = 0
#if 0
            H(i,j) = -99.
#endif
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
      write(0,'(1000(i1))') (mask(i,j), i=imin,imax)
   end do

   return
   end subroutine print_mask
!EOC

!-----------------------------------------------------------------------

   end module domain

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
