#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: nomads_ascii_out() - produces NS_NOMADS specific output.
!
! !INTERFACE:
   subroutine nomads_ascii_out(loop,n,macro)
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: az,imin,imax,jmin,jmax,kmax,H,iimax,iimin,jjmax,jjmin
#if ( defined(SPHERICAL) || defined(CURVILINEAR) )
    use domain, only: xu,yu,yx,angle,xv,xx,dyu,az,arcd1
#else
    use domain, only: dx,dy
#endif
   use m2d,    only: z,D,U,DU,V,DV,ru,rv
   use variables_3d,    only: dt,hn,uu,hun,vv,hvn,S,T,rho
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)	:: loop,n
   logical, intent(in)	:: macro
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Adolf Stips
!
!  $Log: nomads_ascii_out.F90,v $
!  Revision 1.1  2002-05-02 14:01:53  gotm
!  Initial revision
!
!
! !LOCAL VARIABLES:
   REALTYPE	:: pi=3.141592654
   REALTYPE	:: area,V,Vtot,Vup,Vlow
   REALTYPE	:: Taver,Saver,Tupp,Supp,Tlow,Slow
   integer	:: i,j,k,k20
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'nomads_ascii_out() # ',Ncall
#endif

!  Calculate the mean salinity and temperature

#ifdef NS_NOMADS_TEST
   Tup = _ZERO_
   Tlow = _ZERO_
   Taver = _ZERO_
   Supp = _ZERO_
   Slow = _ZERO_
   Saver = _ZERO_
   Vup = _ZERO_
   Vtot = _ZERO_
   do j=jmin,jmax
      do i=imin,imax
         if(az(i,j) .gt. 0) then
	    area = _ONE_/arcd1(i,j)
	    depth = 0.5*hn(i,j,kmax)
	    if(depth .gt. 20.) then
	       k20 = kmax
	       Tupp = Tupp+20.*area*T(i,j,kmax)
	       Supp = Supp+20.*area*S(i,j,kmax)
	       Vup = Vup+20.*area
	    else
	       do k=kmax-1,1,-1
	          depth = depth + 0.5*(hn(i,j,k+1)+hn(i,j,k)) 
	          if(depth .ge. 20.0) EXIT
	       end do
	       k20 = k
            end if
            do k=kmax,k20+1,-1
	       V = area*hn(i,j,k)
	       Tupp = Tupp+V*T(i,j,k)
	       Supp = Supp+V*S(i,j,k)
	       Vup = Vup+V
	    end do
	    V = area*(20.0-20.0)
	    Tupp = Tupp+V*T(i,j,k)
	    Supp = Supp+V*S(i,j,k)
	    Vup = Vup+V
#if 0
            do k=kmax,1,-1
	       V = area*hn(i,j,k)
	       Vup = area*20.0
	       Vlow = V-Vup
	       Taver = vol*T(i,j,k)
	       Saver = vol*S(i,j,k)
            end do
#endif
         end if
      end do
   end do
   Tupp = Tupp/Vup
   Supp = Supp/Vup
   Taver = Taver/Vtot
   Saver = Saver/Vtot

   STDERR Tupp,Taver
   STDERR Supp,Saver

#if 0
! Mean kinetic energy:
    if (abs(loop/M-loop/float(M)).lt.1e-10) then
    MKE=0.
    APE=0.
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
	     MKE=MKE+0.5*1025.858*uu(i,j,k)**2/hun(i,j,k)*dx*dy+vv(i,j,k)**2/hvn(i,j,k)*dx*dy
	  end do
       end do
    end do
    do j=5,jmax-4
       do i=5,imax-4
          zs=z(i,j)
          do k=kmax,1,-1
             zzz=-(kmax-k+0.5)*(20./float(kmax))
	     zs=zs-0.5*hn(i,j,k)
	     densi=1025.-rho(i,j,k)*1025./9.82
!	     write(90,*) i,j,k,S(i,j,k),densi
	     APE=APE+dx*dy*9.82*(hn(i,j,k)*densi*zs-20./float(kmax)*1025.858*zzz)
	     zs=zs-0.5*hn(i,j,k)
	  end do
       end do
    end do
#endif

!    call ape_calc(buoyvec,zvec,volvec,count,APE,POT,POTX)

   return
   end subroutine nomads_ascii_out

!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2002 - Karsten Bolding & Adolf Stips                   !
!-----------------------------------------------------------------------
