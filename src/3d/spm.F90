!$Id: spm.F90,v 1.1.1.1 2002-05-02 14:00:59 gotm Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  suspended_matter
!
! !INTERFACE:
   module suspended_matter
!
! !DESCRIPTION:
!  Description still missing
!
! !USES:
   use commhalo, only: myid,update_3d_halo,wait_halo,D_TAG
   use domain, only: imin,jmin,imax,jmax,H,az
   use domain, only: iimin,jjmin,iimax,jjmax,kmax
   use parameters, only: rho_0
   use variables_3d, only: hn,spm,spm_ws,spm_pool,taub
   IMPLICIT NONE
!
   private
!
! !PUBLIC DATA MEMBERS:
   public init_spm, do_spm
!
! !PRIVATE DATA MEMBERS:
   integer		:: spm_method=1
   integer              :: spm_init_method=1, spm_format=2
   character(len=32)	:: spm_file="spm.nc",spm_name='spm'
   integer		:: spm_hor_adv=1,spm_ver_adv=1,spm_strang=0
   REALTYPE             :: spm_AH = -_ONE_
   REALTYPE             :: spm_const=_ZERO_  
   REALTYPE             :: spm_init= _ZERO_
   integer              :: spm_ws_method = 0
   REALTYPE             :: spm_ws_const=0.001
   REALTYPE             :: spm_erosion_const, spm_tauc_sedimentation
   REALTYPE             :: spm_tauc_erosion, spm_pool_init
!  For erosion-sedimentation flux
   REALTYPE             :: Erosion_flux , Sedimentation_flux
   logical              :: erosed_flux =.false.
!
! !REVISION HISTORY:
!  Original author(s): Manuel Ruiz Villarreal, Karsten Bolding and Hans Burchard
!
!  $Log: spm.F90,v $
!  Revision 1.1.1.1  2002-05-02 14:00:59  gotm
!  recovering after CVS crash
!
!
!  Revision 1.2  2001/10/23 07:30:19  bbh
!  module spm renamed to suspended_matter to allow for variable named spm
!
!  Revision 1.1  2001/10/22 11:16:09  bbh
!  Added support for particulate suspended matter - no IO yet
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_spm
!
! !INTERFACE:
   subroutine init_spm(adv_method)
!
! !DESCRIPTION:
!  Description still missing
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)	:: adv_method 
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  See the log for the module
!
! !LOCAL VARIABLES:
   integer	:: i,j,k,n
   integer      :: rc
   NAMELIST /spm_nml/	spm_method,spm_init_method,  &
                        spm_const,spm_format,spm_file,spm_name, &
                        spm_hor_adv, spm_ver_adv,spm_AH,spm_strang, &
                        spm_ws_method,spm_ws_const,                   &
			spm_erosion_const, spm_tauc_sedimentation, &
			spm_tauc_erosion, spm_pool_init
   integer, parameter	:: nmax=10000
   REALTYPE		:: zlev(nmax),prof(nmax)			
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_spm() # ',Ncall
#endif

   LEVEL2 'init_spm()'
   read(NAMLST,spm_nml)

   LEVEL3 'spm_hor_adv= ',spm_hor_adv
   LEVEL3 'spm_ver_adv= ',spm_ver_adv
   LEVEL3 'spm_strang=  ',spm_strang

   spm_ws = _ZERO_
       
!  Compute settling velocity
   select case(spm_ws_method)
      case(0) !constant
         do k=1,kmax-1
	    spm_ws(:,:,k) = spm_ws_const ! positive downwards
         end do
!     Zanke formula for fall velocity of sediment 
      case(1) 
         FATAL 'Zanke method for settling velocity not yet coded'
         stop  'init_spm'
   end select 
   
   select case (spm_init_method)
      case(0)
         LEVEL3 'getting initial fields from hotstart'
      case(1)
         LEVEL3 'setting to constant value'
            forall(i=iimin:iimax,j=jjmin:jjmax, az(i,j) .ne. 0) &
	           spm(i,j,:) = spm_const
      case(2)
         LEVEL3 'using profile'
         call read_profile(spm_file,nmax,zlev,prof,n)
         call ver_interpol(n,zlev,prof,imin,jmin,imax,jmax,az,H,	&
	                   iimin,jjmin,iimax,jjmax,kmax,hn,spm)
      case(3)
         LEVEL3 'interpolating from 3D field'
         call get_field(spm_file,spm_name,spm)			   	
      case default
         FATAL 'Not valid spm_init_method specified'
         stop 'init_spm'
   end select
   
   select case (spm_method)
       case(1) !erosion-sedimentation flux
          erosed_flux = .true.  
	  spm_pool = spm_pool_init
       case default
          FATAL 'Not valid spm_method specified'
         stop 'init_spm'
   end select 	  
#ifdef DEBUG
   write(debug,*) 'Leaving init_spm()'
   write(debug,*)
#endif
   return
   end subroutine init_spm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_spm()
!
! !INTERFACE:
   subroutine do_spm()
!
! !DESCRIPTION:
!  Description still missing
!
! !USES:
   use advection_3d, only: do_advection_3d
   use variables_3d, only: dt,cnpar,hun,hvn,ho,nuh,uu,vv,ww
   use domain,       only: au,av
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: dxu,dxv,dyu,dyv,arcd1
#else
   use domain, only: dx,dy,ard1
#endif
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
   integer	:: i,j,k,rc
   REALTYPE     :: spmtot
   REALTYPE	:: Res(0:kmax)
   REALTYPE	:: auxn(1:kmax-1),auxo(1:kmax-1)
   REALTYPE	:: a1(0:kmax),a2(0:kmax),a3(0:kmax),a4(0:kmax)
   REALTYPE , allocatable  ,dimension (:,:,:) :: ww_aux 
   REALTYPE	:: delxu(I2DFIELD),delxv(I2DFIELD)
   REALTYPE	:: delyu(I2DFIELD),delyv(I2DFIELD)
   REALTYPE	:: area_inv(I2DFIELD)
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_spm() # ',Ncall
#endif
   
#if defined(SPHERICAL) || defined(CURVILINEAR)
   delxu=dxu
   delxv=dxv
   delyu=dyu
   delyv=dyv
   area_inv=arcd1
#else
   delxu=dx
   delxv=dx
   delyu=dy
   delyv=dy
   area_inv=ard1
#endif
   allocate(ww_aux(I3DFIELD),stat=rc)    ! work array
   if (rc /= 0) stop 'init_spm: Error allocating memory (ww_aux)'
   
!  The vertical velocity to be used in the advection routine for spm is ww-ws
   ww_aux = ww-spm_ws
   call do_advection_3d(dt,spm,uu,vv,ww_aux,hun,hvn,ho,hn,   &
                        delxu,delxv,delyu,delyv,area_inv,az,au,av,	&
                        spm_hor_adv,spm_ver_adv,spm_strang,spm_AH)
!  Vertical diffusion of spm 
   do j=jjmin,jjmax
      do i=iimin,iimax
         if (az(i,j) .eq. 1) then
            if (kmax.gt.1) then
!           We impose a flux condition on bottom sediments
!           where flux is computed as the result of erosion and sedimentation
!           from a pool of nondynamic particulate matter   
               if(erosed_flux) then	           
	          if(spm_pool(i,j) > 1.e-12) then ! If there are sediments in the pool
		    Erosion_Flux = spm_erosion_const / rho_0 * (taub(i,j)*rho_0-spm_tauc_erosion )  
                    Erosion_Flux = max(Erosion_Flux,0.) 	                 
                  else
                     Erosion_Flux = _ZERO_
                  end if
                  Sedimentation_Flux = spm_ws(i,j,1) * spm(i,j,1) *(1.-taub(i,j)*rho_0 / spm_tauc_sedimentation)
                  Sedimentation_Flux = max(Sedimentation_Flux,0.)    
                  a4(1) = Erosion_Flux - Sedimentation_Flux
                  if (a4(1) > 1e-12) a4(1) = min(spm_pool(i,j)/dt, a4(1))
		  spm_pool(i,j) = spm_pool(i,j) -  dt * a4(1)
                  a4(1) = a4(1)* dt          
               else
                  a4(1) = _ZERO_  
               end if

!     Auxiliary terms, old and new time level,
               do k=1,kmax-1
	          auxo(k)=2.*(1-cnpar)*dt*nuh(i,j,k)/(hn(i,j,k+1)+hn(i,j,k))
	          auxn(k)=2.*   cnpar *dt*nuh(i,j,k)/(hn(i,j,k+1)+hn(i,j,k))
	       end do

!        Matrix elements for surface layer
	       k=kmax
	       a1(k)=-auxn(k-1)
	       a2(k)=hn(i,j,k)+auxn(k-1)
	       a4(k)=spm(i,j,k)*(hn(i,j,k)-auxo(k-1))+spm(i,j,k-1)*auxo(k-1) 

!        Matrix elements for inner layers
               do k=2,kmax-1
	          a3(k)=-auxn(k  )
	          a1(k)=-auxn(k-1)
       	          a2(k)=hn(i,j,k)+auxn(k)+auxn(k-1)
	          a4(k)=spm(i,j,k+1)*auxo(k)                          &
	               +spm(i,j,k  )*(hn(i,j,k)-auxo(k)-auxo(k-1))    &
	               +spm(i,j,k-1)*auxo(k-1)              
               end do 		    

!        Matrix elements for bottom layer
               k=1
	       a3(k)=-auxn(k  )
	       a2(k)=hn(i,j,k)+auxn(k)                            
               a4(k)=a4(k)+spm(i,j,k+1)*auxo(k)                           &
	            +spm(i,j,k  )*(hn(i,j,k)-auxo(k))                          

               call getm_tridiagonal(kmax,1,kmax,a1,a2,a3,a4,Res)
	       do k=1,kmax 
	          spm(i,j,k)=Res(k)
	       end do    
            end if 
         end if 
      end do
   end do    

#ifdef SALTWEDGE_TEST
!  Valid for saltwedge case, lateral zero-gradient BC for WE bound 
!  call do_bdy_3d(2,spm)
   spm(iimin,2,:)=spm(iimin+1,2,:) 
#endif
#ifdef AAAA
!  Checks conservation
   spmtot= _ZERO_  
   do k=1,kmax
      do j=jjmin,jjmax
         do i=iimin,iimax
            if (az(i,j) .eq. 1) then
               spmtot=spmtot+spm(i,j,k)*hn(i,j,k) 
            end if    
         end do
      end do
   end do
#endif

   call update_3d_halo(spm,spm,az,iimin,jjmin,iimax,jjmax,kmax,D_TAG)

   call wait_halo(D_TAG)

#ifdef FORTRAN90
   deallocate(ww_aux,stat=rc)
   if (rc /= 0) stop 'upstream_adv: Error de-allocating memory (ww_aux)'
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving do_spm()'
   write(debug,*)
#endif
   return
   end subroutine do_spm
!EOC

!-----------------------------------------------------------------------

   end module suspended_matter

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Manuel Ruiz, Hans Burchard and Karsten Bolding  !
!-----------------------------------------------------------------------
