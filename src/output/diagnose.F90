#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagnose() - produces a lot of diagnostic output.
!
! !INTERFACE:
   subroutine Diagnose(loop,n,macro)
!
! !DESCRIPTION:
!
! !USES:
   use domain, only: az,imin,imax,jmin,jmax,kmax,H,iimax,iimin,jjmax,jjmin
#if ! ( defined(SPHERICAL) || defined(CURVILINEAR) )
    use domain, only: dx,dy
#else
    use domain, only: dyu
#if defined(CURVILINEAR)
    use domain, only: xu,yu,yx,angle,xv,xx,dyu
#endif
#endif
   use m2d,    only: z,D,U,DU,V,DV,ru,rv
!HB   use m3d,    only: hn,uu,hun,vv,hvn,S,M,rho
#ifndef NO_3D
   use m3d,    only: M
   use variables_3d,    only: dt,hn,uu,hun,vv,hvn
#ifndef NO_BAROCLINIC
   use variables_3d, only: S,T,rho
#endif
#endif
   use meteo,      only: tausx
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: loop,n
   logical, intent(in)                 :: macro
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: diagnose.F90,v $
!  Revision 1.5  2006-01-27 19:55:56  hb
!  tab removed
!
!  Revision 1.4  2004-01-07 07:37:37  kbk
!  to compile with IFORT - TABS, etc.
!
!  Revision 1.3  2003/04/23 12:07:12  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.2  2003/04/07 12:32:58  kbk
!  parallel support + NO_3D, NO_BAROCLINIC
!
!  Revision 1.1.1.1  2002/05/02 14:01:52  gotm
!  recovering after CVS crash
!
!  Revision 1.1  2001/08/27 11:54:33  bbh
!  TVD-advection for momentum added, some bugs removed
!
!  Revision 1.1.1.1  2001/04/17 08:43:07  bbh
!  initial import into CVS
!
!
! !LOCAL VARIABLES:
   integer, parameter        :: MaxGauges=10
   integer                   :: gauges=0
   integer                   :: i,j,k,l
   integer:: ii(MaxGauges),jj(MaxGauges)
   REALTYPE:: zz,u1,v1,min_depth,uup,udn,vup,vdn,velup,veldn
   REALTYPE:: DMIN,DMAX
   integer :: DMINI,DMINJ,DMAXI,DMAXJ
   REALTYPE:: UMIN,UMAX
   integer :: UMINI,UMINJ,UMAXI,UMAXJ
   REALTYPE:: DUMIN,DUMAX
   integer :: DUMINI,DUMINJ,DUMAXI,DUMAXJ
   REALTYPE:: VMIN,VMAX
   integer :: VMINI,VMINJ,VMAXI,VMAXJ
   REALTYPE:: DVMIN,DVMAX
   integer :: DVMINI,DVMINJ,DVMAXI,DVMAXJ
   REALTYPE:: pi=3.141592654
   REALTYPE:: Flux1,Flux2,Flux3
#ifdef NOMADS_TEST
   REALTYPE:: MKE,APE,densi,zzz,zs,salmin,salmax,area
   INTEGER     :: salimax,saljmax
#endif
#ifdef CONSTANCE_TEST
   REALTYPE:: ZET,MKE,APE,POT,POTX,APE0,densi,zs,Tint,ccc,N0,b0,WIND=0.
   LOGICAL,save:: FIRST=.true.
   INTEGER     :: out,count
   integer      :: rc
   REALTYPE,save,dimension(:,:,:), allocatable:: dens0
   REALTYPE buoyvec(1:1000000)
   REALTYPE zvec(1:1000000)
   REALTYPE volvec(1:1000000)
#endif
   REALTYPE aa(1:500,0:25)
#ifdef HAIDVOGEL_TEST
   REALTYPE        :: MKE,APE,APE0,densi,zs,zzz
   LOGICAL,save    :: FIRST=.true.
#endif
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'Diagnose() # ',Ncall
#endif

#ifdef SYLT_TEST
   gauges = 8
!  (* Lister Hafen     i= 80; j= 52; h=  6.8   *)
!  (* Havneby          i=118; j= 94; h=  2.8   *)
!  (* Messpf. Westerl. i= 20; j=  1; h=  4.7   *)
!  (* Ballum Schleuse  i=141; j=127; h=- 0.6   *)
!  (* Munkmarsch       i= 30; j= 27; h=  0.6   *)
!  (* Rickelsbuell     i= 22; j=114; h=  0.0   *)
!  (* Pandertief       i= 36; j= 31; h=  1.8   *)
!  (* Buttersand       i= 28; j= 67; h=  1.9   *)

   ii(1)= 52
   jj(1)= 80
   ii(2)= 94
   jj(2)=118
   ii(3)=  1
   jj(3)= 20
   ii(4)=127
   jj(4)=141
   ii(5)= 27
   jj(5)= 30
   ii(6)=114
   jj(6)= 22
   ii(7)= 31
   jj(7)= 36
   ii(8)= 67
   jj(8)= 28
#endif

! Calculate min and max at every micro time step

   min_depth= 1000.
   Dmax=-1000.
   do j=jmin,jmax
      do i=imin,imax
         if (D(i,j).gt.Dmax) then
            Dmax=D(i,j)
            Dmaxi=i
            Dmaxj=j
         end if
         if (D(i,j).lt.min_depth) then
            min_depth=D(i,j)
            Dmini=i
            Dminj=j
         end if
      end do
   end do

   DUmin= 1000.
   DUmax=-1000.
   Umin= 1000.
   Umax=-1000.
   do j=jmin,jmax
      do i=imin,imax
         if (DU(i,j).gt.DUmax) then
            DUmax=DU(i,j)
            DUmaxi=i
            DUmaxj=j
         end if
         if (DU(i,j).lt.DUmin) then
            DUmin=DU(i,j)
            DUmini=i
            DUminj=j
         end if
         if (U(i,j)/DU(i,j).gt.Umax) then
            Umax=U(i,j)/DU(i,j)
            Umaxi=i
            Umaxj=j
         end if
         if (U(i,j)/DU(i,j).lt.Umin) then
            Umin=U(i,j)/DU(i,j)
            Umini=i
            Uminj=j
         end if
      end do
   end do

   DVmin= 1000.
   DVmax=-1000.
   Vmin= 1000.
   Vmax=-1000.
   do j=jmin,jmax
      do i=imin,imax
         if (DV(i,j).gt.DVmax) then
            DVmax=DV(i,j)
            DVmaxi=i
            DVmaxj=j
         end if
         if (DV(i,j).lt.DVmin) then
            DVmin=DV(i,j)
            DVmini=i
            DVminj=j
         end if
         if (V(i,j)/DV(i,j).gt.Vmax) then
            Vmax=V(i,j)/DV(i,j)
            Vmaxi=i
            Vmaxj=j
         end if
         if (V(i,j)/DV(i,j).lt.Vmin) then
            Vmin=V(i,j)/DV(i,j)
            Vmini=i
            Vminj=j
         end if
      end do
   end do

#ifdef CURVITEST
     i=5
     j=5
     Flux1=0.
     Flux2=0.
     Flux3=0.
     write(92,*) loop,U(i,j)/DU(i,j),z(i,j)
     do j=2,jmax-1
        i=1
        Flux1=Flux1+U(i,j)*DYU
        i=imax/3
        Flux2=Flux2+U(i,j)*DYU
        i=imax/2
        Flux3=Flux3+U(i,j)*DYU
     end do
     write(97,*) loop,Flux1,Flux2,Flux3

#ifdef CURVILINEAR
#ifdef PECS
     if (abs(loop/100-loop/100.).lt.1.e-10) then
        write(80,*) 'loop #',loop
       do i=imin,imax
          do j=jmin,jmax
              if (S(i,j,1).gt.13) then
                write(80,*) xx(i,j),yx(i,j)
                write(80,*) xx(i-1,j),yx(i-1,j)
                write(80,*) xx(i-1,j-1),yx(i-1,j-1)
                write(80,*) xx(i,j-1),yx(i,j-1)
                write(80,*) '>'
              end if        	
          end do
       end do
     end if
#endif
     if (loop.eq.n) then
        i=imax/2
        do j=2,jmax-1
            write(90,*) U(i,j)/DU(i,j),0.5*(z(i,j)+z(i+1,j)),yu(i,j)
!            write(90,*) U(i,j)/DU(i,j),0.5*(z(i,j)+z(i+1,j)),xu(i,j)
        end do
        do j=2,jmax-2
            write(91,*) V(i,j)/DV(i,j)+V(i+1,j)/DV(i+1,j),yx(i,j)
!            write(91,*) V(i,j)/DV(i,j)+V(i+1,j)/DV(i+1,j),xx(i,j)
        end do
       j=jmax/2
        do i=1,imax
           uup=0.5*(U(i-1,j+1)/DU(i-1,j+1)+U(i,j+1)/DU(i,j+1))
           udn=0.5*(U(i-1,j)/DU(i-1,j)+U(i,j)/DU(i,j))
           vup=0.5*(V(i,j+1)/DV(i,j+1)+V(i,j)/DV(i,j))
           vdn=0.5*(V(i,j-1)/DV(i,j-1)+V(i,j)/DV(i,j))
           velup=sqrt(uup**2+vup**2)
           veldn=sqrt(udn**2+vdn**2)
           uup=velup*cos(2.*pi*angle(i,j)/360.)
           vup=velup*sin(2.*pi*angle(i,j)/360.)
            write(93,*) 0.5*(vup+vdn),0.5*(z(i,j)+z(i,j+1)),xv(i,j)
            write(94,*) 0.5*(uup+udn),xx(i,j)
        end do
     end if
#else
     if (loop.eq.n) then
        j=jmax/2
        do i=1,imax-1
           write(95,*) V(i,j)/DV(i,j),0.5*(z(i,j)+z(i,j+1)),(i-0.5)*dx
        end do        
        do i=1,imax-2
           write(96,*) 0.5*(U(i,j)/DU(i,j)+U(i,j+1)/DU(i,j+1)),i*dx
        end do        
        i=imax/2
        do j=2,jmax-1
            write(90,*) U(i,j)/DU(i,j),0.5*(z(i,j)+z(i+1,j)),(j-1-0.5)*dy
        end do
        do j=2,jmax-2
            write(91,*) V(i,j)/DV(i,j)+V(i+1,j)/DV(i+1,j),(j-1)*dy
        end do
       j=jmax/2
        do i=1,imax
            write(93,*) V(i,j)/DV(i,j),0.5*(z(i,j)+z(i,j+1)),(i-0.5)*dx
        end do
        do i=1,imax
            write(94,*) 0.5*(U(i,j)/DU(i,j)+U(i,j+1)/DU(i,j+1)),i*dx
        end do
     end if
!#endif
#endif
#endif
#ifdef HAIDVOGEL_TEST

    if (abs(loop/M-loop/float(M)).lt.1e-10) then
    MKE=0.
    APE=0.
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
            MKE=MKE+0.5*1025.*(uu(i,j,k)**2/hun(i,j,k)*dx*dy+vv(i,j,k)**2/hvn(i,j,k)*dx*dy)
         end do
       end do
    end do
    do j=1,jmax
       do i=1,imax
          zs=z(i,j)
          do k=kmax,1,-1
            zs=zs-0.5*hn(i,j,k)
            densi=1025.-rho(i,j,k)*1025./9.82
            APE=APE+dx*dy*9.82*hn(i,j,k)*densi*zs
            zs=zs-0.5*hn(i,j,k)
         end do
       end do
    end do
    if (FIRST) then
       APE0=APE
       FIRST=.false.
    end if
    zzz=0
    do i=2,129
       do j=2,4
          zzz=zzz+z(i,j)
       end do
    end do        
    write(80,996) loop*dt/3600./24./float(M),MKE,APE-APE0,MKE+APE-APE0
    end if
#endif

#ifdef NOMADS_TEST

! Mean kinetic energy:
    if (abs(loop/M-loop/float(M)).lt.1e-10) then
!do i=1,imax
!   do j=1,jmax
!      write(85,*) i,j,T(i,j,kmax),S(i,j,kmax),1025.-rho(i,j,kmax)*1025./9.82
!   end do
!end do
!stop
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
!            write(90,*) i,j,k,S(i,j,k),densi
            APE=APE+dx*dy*9.82*(hn(i,j,k)*densi*zs-20./float(kmax)*1025.858*zzz)
            zs=zs-0.5*hn(i,j,k)
         end do
       end do
    end do
    salmin=10000.
    salmax=0.
    area=0.
    do j=1,jmax
       do i=1,imax
          do k=kmax,1,-1
            if (S(i,j,k).gt.salmax) then
               salmax=S(i,j,k)
               salimax=i
               saljmax=j
            end if
            if (S(i,j,k).lt.salmin) salmin=S(i,j,k)
            if ((S(i,j,k).lt.34.839).and.(k.eq.kmax)) area=area+dx*dy
         end do
       end do
    end do
    write(99,994) loop,MKE,APE,MKE+APE,salmin,salmax,area
!    write(*,994) loop,MKE,APE,MKE+APE,salmin,salmax,area
!    stop
    write(98,*) salimax,saljmax,salmax
    end if
    if (loop.eq.n) then
       do i=5,imax-4
          do j=5,jmax-4
            write(96,995) (i-4.5)*30./float(imax-8),            &
                          (j-4.5)*30./float(jmax-8),S(i,j,kmax),     &
           0.5*(uu(i,j,kmax)/hun(i,j,kmax)+uu(i-1,j,kmax)/hun(i-1,j,kmax)), &
           0.5*(vv(i,j,kmax)/hvn(i,j,kmax)+vv(i,j-1,kmax)/hvn(i,j-1,kmax)), &
           0.5*(uu(i,j,1)/hun(i,j,1)+uu(i-1,j,1)/hun(i-1,j,1)), &
           0.5*(vv(i,j,1)/hvn(i,j,1)+vv(i,j-1,1)/hvn(i,j-1,1))
        end do
       end do
    end if
#endif

#if 0
!    j=68
    i=109
    do j=32,51
!    do i=112,138
       if ((H(i,j).lt.-5).and.(H(i,j+1).gt.-5)) then
          aa(j,0)=-H(i,j+1)
          do k=1,kmax
             aa(j,k)=aa(j,k-1)+hn(i,j+1,k) 
          end do 
       end if
       if ((H(i,j).gt.-5).and.(H(i,j+1).lt.-5)) then
          aa(j,0)=-H(i,j)
          do k=1,kmax
             aa(j,k)=aa(j,k-1)+hn(i,j,k) 
          end do 
       end if
       if ((H(i,j).gt.-5).and.(H(i,j+1).gt.-5)) then
          aa(j,0)=-0.5*(H(i,j)+H(i,j+1))
          do k=1,kmax
             aa(j,k)=aa(j,k-1)+0.5*(hn(i,j,k)+hn(i,j+1,k)) 
          end do 
       end if 
    end do
    do k=1,kmax
!       do i=113,138
        do j=33,51
          if (H(i,j).gt.-5) then 
          write(96,*) j*6.,aa(j,k)
          write(96,*) (j-1.)*6.,aa(j-1,k)
          write(96,*) (j-1.)*6.,aa(j-1,k-1)
          write(96,*) j*6.,aa(j,k-1)
          write(96,*) j*6.,aa(j,k)
          write(96,*) 
         end if 
       end do 
    end do
    stop
#endif

#ifdef CONSTANCE_TEST
! Mean kinetic energy:
    if (abs(loop/M-loop/float(M)).lt.1e-10) then
    MKE=0.
    APE=0.
    ZET=0.
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
            MKE=MKE+uu(i,j,k)**2/hun(i,j,k)+vv(i,j,k)**2/hvn(i,j,k)
         end do
       end do
    end do
    MKE=MKE*0.5*dx*dy*1025.
    ccc=9.82/1025.*0.17
    if (FIRST) then
       allocate(dens0(I3DFIELD),stat=rc)    ! work array
       if (rc /= 0) stop 'diagnose: Error allocating memory (dens0)'
       dens0=rho
    end if
    do j=2,jmax-1 ! Calculate available potential energy
       do i=2,imax-1
          ZET=ZET+dx*dy*9.82*1025.*0.5*z(i,j)**2
          WIND=WIND+1025.*dt*dx*dy*tausx(i,j)*uu(i,j,kmax)/hun(i,j,kmax)
       end do
    end do
    if (FIRST) then
       FIRST=.false.
    end if
    if ((abs(loop/14400-loop/14400.).lt.1e-10).or.(loop.eq.10)) then
       do i=imin,imax
          do j=jmin,jmax
           out=85+loop/14400
            write(out,995) (i-0.5)*0.5,(j-0.5)*0.5,T(i,j,kmax),     &
            0.5*(uu(i,j,kmax)/hun(i,j,kmax)+uu(i-1,j,kmax)/hun(i-1,j,kmax)), &
            0.5*(vv(i,j,kmax)/hvn(i,j,kmax)+vv(i,j-1,kmax)/hvn(i,j-1,kmax)), &
           z(i,j)
        end do
       end do
    end if
    i=12
    j=27
    k=kmax
    write(99,*) loop*dt/3600./24./float(M),vv(i,j,k)/hvn(i,j,k)
    write(*,*) loop*dt/3600./24./float(M),vv(i,j,k)/hvn(i,j,k)


!    count=0
!    do j=1,jmax
!       do i=1,imax
!          zz=-H(i,j)
!          do k=1,kmax
!             if (az(i,j).eq.1) then
!                count=count+1
!                buoyvec(count)=rho(i,j,k)
!                zz=zz+0.5*hn(i,j,k)
!                zvec(count)=zz
!                zz=zz+0.5*hn(i,j,k)
!                volvec(count)=dx*dy*hn(i,j,k)
!             end if
!          end do
!       end do
!    end do


!    call ape_calc(buoyvec,zvec,volvec,count,APE,POT,POTX)


!    write(80,996) loop*dt/3600./24./float(M),MKE,WIND,ZET,APE,POT,POTX

    end if

#endif

 994  format(I5,1x,8(E13.6,1x))
 996  format(F10.5,1x,8(E13.6,1x))
 995  format(8(F10.5,1x))
 993  format(6(F12.8,1x))
 999  format (5F10.5)

      return
      end subroutine Diagnose


      subroutine ape_calc(buoyvec,zvec,volvec,count,APE,POT,POTX)

      implicit none

! INPUT VARIABLES
      REALTYPE buoyvec(1:1000000) ! Vector with buoyancy in m/s**2
      REALTYPE zvec(1:1000000)    ! Vector with z-levels of tracer points in m
      REALTYPE volvec(1:1000000)  ! Vector with volumes in m**3
      integer count               ! Length of vectors

! OUTPUT VARIABLES
      REALTYPE APE      ! Available potential energy
      REALTYPE POT      ! Potential energy, old definition
      REALTYPE POTX     ! Potential energy, new definition

! LOCAL VARIABLES
      REALTYPE,save :: POT0   ! Initial potential energy, old definition
      integer i,j,k
      REALTYPE vol,volidxi    ! dummy variables for volumes
      REALTYPE b0(1:1000000)  ! Buoyancy profile with minimum potential energy
      REALTYPE NN0(1:1000000) ! N**2 of b0
      REALTYPE,save :: b00(1:1000000)  ! Initial b0
      REALTYPE,save :: NN00(1:1000000) ! Initial N**2
      REALTYPE densi                   ! dummy variable for density

      integer idx(1:count)             ! index for increasing buoyancy
      integer zidx(1:count)            ! index for increasing z-levels

      REALTYPE :: grav=9.82
      REALTYPE :: rho0=1025.

      logical :: FIRST=.TRUE.


      call sort(buoyvec,count,idx)     ! Sort the buoyancy
      call sort(zvec,count,zidx)       ! Sort the depths

      b0=0.
      j=1
      vol=volvec(zidx(1))
      do i=1,count                     ! Redistribute buoyancy
         volidxi=volvec(idx(i))
111      if (vol.ge.volidxi) then
            b0(j)=b0(j)+buoyvec(idx(i))*volidxi
            vol=vol-volidxi
         else
            b0(j)=b0(j)+buoyvec(idx(i))*vol
            volidxi=volidxi-vol
            if (j.lt.count) then
               j=j+1
               vol=volvec(zidx(j))
               goto 111
            else
            end if
         end if
      end do
      do i=1,count
         b0(i)=b0(i)/volvec(zidx(i))
      end do

      APE=0.
      POT=0.
      do i=1,count                  ! Calculate potential energies
         j=i  ! lower index for NN0 calculation
         if (j.gt.1) then
222         j=j-1
            if ((j.gt.1).and.(abs(zvec(zidx(j))-zvec(zidx(i))).lt.1.e-1.or. &
                abs(b0(j)-b0(i)).lt.1.e-8 )) goto 222
         end if
         k=i   ! upper index for NN0 calculation
         if (k.lt.count) then
333        k=k+1
           if ((k.lt.count).and.(abs(zvec(zidx(k))-zvec(zidx(i))).lt.1.e-1.or. &
                abs(b0(k)-b0(i)).lt.1.e-8 )) goto 333
         end if
         NN0(i)=(b0(k)-b0(j))/(zvec(zidx(k))-zvec(zidx(j)))
         APE=APE+rho0*volvec(zidx(i))*(buoyvec(zidx(i))-b0(i))**2/NN0(i)
         densi=rho0-buoyvec(zidx(i))*rho0/grav
         POT=POT+grav*volvec(zidx(i))*zvec(zidx(i))*densi
      end do

      if (FIRST) then
         POT0=POT
         b00=b0
         NN00=NN0
         FIRST=.FALSE.
      end if
      POT=POT-POT0

      POTX=0.    ! Calculate potential energy relative to initial profile
      do i=1,count
         POTX=POTX+rho0*volvec(zidx(i))*(buoyvec(zidx(i))-b00(i))**2/NN00(i)
      end do


      return
      end subroutine ape_calc


      subroutine sort(vec,count,idx)

      implicit none

      integer count,vec_length,i,j,minj
      integer idx(1:count)
      REALTYPE vec(1:1000000),minvec
      logical sorted(1:1000000)

      sorted=.FALSE.
      do i=1,count
         minvec=1.e10
         do j=1,count
            if ((vec(j).le.minvec).and.(.not.sorted(j))) then
               minvec=vec(j)
               minj=j
            end if
         end do
         sorted(minj)=.true.
         idx(i)=minj
      end do

      return
      end subroutine sort
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 1999 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
