c###
	subroutine otide(model_name,tide,lat,long,jnasa,heure,istat)
c------------------------------------------------------------------------------
c  PROGRAMME : otide.f
c
c  DESCRIPTION : tide prediction routine (end user routine):
c	   tide= predicted tide in cm (output)
c          lat = latitude (between -90 et +90 deg ) (input)
c          long = longitude (either between -180 180 or 0 360).
c          jnasa= Day for prediction (as in TP GDR )
c          heure= Hour during the day jnasa
c          istat= Number of valid points in the interpolation
c
c  PROGRAMMEURS : J.M. MOLINES 
c 
c  DATE derniere retouche : 13/13/95
c------------------------------------------------------------------------------
        include 'common.h' 
c
	real*8 t1,dlat,dlon,heure
	character*(*) model_name
	real*4 lat,long,h1(extended_max)
	data ifirst /0/
c read the data files at first call. This is quite long
	if (ifirst.eq.0) then
	model=model_name
	 call init_data
	 ifirst=1
	endif
c
c       conversion  TOPEX date to adequate date for prediction
c
	call convert_date(jnasa,heure,t1)
c
	if (long.ge.180.) long=long-360.
	 dlon=dble(long)
	 dlat=dble(lat)

        call grenoble(t1,long,lat,h1,istat)
	tide=h1(13)
	return
	end
c##
                subroutine astronomics(tj)                                          
c------------------------------------------------------------------------------
c  PROGRAMME : astronomics.f
c
c  DESCRIPTION : This program initialize some astronomic data usefull for
c                nodal corrections.
c
c  PROGRAMMEURS : J.M. MOLINES et F. LYARD
c 
c  DATE derniere retouche : 19/01/94
c------------------------------------------------------------------------------
       include 'common.h'          
                      
      real*8 tj,rad,tgn2,at1,at2,u,tgi2
      real*8 ct0,ct1
      real*8 cn0,cn1
      real*8 cs0,cs1
      real*8 ch0,ch1
      real*8 cps0,cps1
      real*8 cp0,cp1  

c------------- tt mean solar angle relative to Greenwich ---------------
                    
	ct0=180.D+00
	ct1=360.D+00*3.6525D+04
      tt=dmod(ct0+ct1*tj,360.D+00)

c-------------- hp longitude of ascending lunar node ----------------
                                                           
	cn0=  259.1560563D+00
	cn1= -1934.1423972D+00   
      n=dmod(cn0+cn1*tj,360.D+00)

c------------------- hp mean solar longitude ----------------------

	ch0= 280.1895015D+00
	ch1= 36000.76892D+00 
	hp=dmod(ch0+ch1*tj,360.D+00)
                           
c-------------------- s mean lunar longitude ----------------------

	cs0= 277.0256206D+00
	cs1= 481267.892D+00
	s=dmod(cs0+cs1*tj,360.D+00)

c------------------- p1 longitude of solar perigee --------------------
                     
	cps0=281.2208568D+00
	cps1=1.719175D+00   
	p1=dmod(cps0+cps1*tj,360.D+00)

c------------------- p longitude of lunar perigee --------------------

	cp0=334.3837214D+00
	cp1=4069.0322056D+00
	p=dmod(cp0+cp1*tj,360.D+00)

      n=rad(n)
      hp=rad(hp)
      s=rad(s)
      p1=rad(p1)
      p=rad(p)

      u=9.13694997D-01-3.5692561D-02*dcos(n)
      iang=dacos(u)
      tgn2=dtan(n/2.)
      at1=datan(1.01883D+00*tgn2)
      at2=datan(6.4412D-01*tgn2)
      xi=-at1-at2+n
      if(n.gt.pi) then
        xi=xi-2.D+00*pi
      endif
      nu=at1-at2
   
c------------------------for constituents l2,k1,k2------------------------


      tgi2=dtan(iang/2.D+00)     
      pp=p-xi
      x1ra=dsqrt(1.-12.D+00*tgi2**2*cos(2.D+00*pp)+36.D+00*tgi2**4)
      r=datan(dsin(2.*pp)/(1./(6.*tgi2**2)-dcos(2.*pp)))
      nuprim=datan(dsin(2.D+00*iang)*dsin(nu)/
     &            (dsin(2.D+00*iang)*dcos(nu)+3.347D-01)) 
      nusec=0.5*datan(((dsin(iang)**2.D+00)*dsin(2.D+00*nu))
     &   /(dsin(iang)**2.D+00*dcos(2.*nu)+7.27D-02))                 

c------------------------------------------------------------------------
      return
      end
c###
      subroutine calend1(njd,nd,nm,na)
c------------------------------------------------------------------------------
c  PROGRAMME : calend1.f
c
c  DESCRIPTION : Conversion from Jcnes to  dd/mm/aaaa
c
c
c  PROGRAMMEURS : Origine GRGS Toulouse, C. BROSSIER
c 
c  DATE derniere retouche : 19/01/94
c------------------------------------------------------------------------------
      dimension n(12)
      data n /31,28,31,30,31,30,31,31,30,31,30,31/
      njul=njd+1
      na=njul/365
      nj=njul-na*365
      nb=(na+1)/ 4
      nj=nj-nb
      if(nj.gt.0) go to 1
      na=na+1949
      nm=12
      nd=nj+31
      return
    1 j=na-2-nb*4
      na=na+1950
      if(j.lt.0) go to 5000
 4000 if(60-nj)4500,7000,5000
 4500 nm1=60
      m=3
      go to 6000
 5000 nm1=0
      m=1
 6000 ndj=nm1+n(m)
      nj3=nj-ndj
      if(nj3.le.0) go to 8000
 6500 m=m+1
      nm1=ndj
      go to 6000
 7000 nm=2
      nd=29
      return
 8000 nm=m
      nd=nj-nm1
 9000 return
      end
c###
	subroutine convert_date(jnasa,heure,t1)
c------------------------------------------------------------------------------
c  PROGRAMME : convert_date.f
c
c  DESCRIPTION : Convert from Jnasa to Julian century as used in Schureman
c
c
c  PROGRAMMEURS : J.M. MOLINES 
c 
c  DATE derniere retouche : 13/03/1995
c------------------------------------------------------------------------------
c	
	integer tojul
	real*8 heure,t1,julian
	real*4 sec
c
c CAUTION 2922 = Number of days between NASA and CNES ...
c
	call calend1(jnasa+2922,ida,imonth,ian)
	iday=tojul(ian,imonth,ida)
        j=iday-1000*(ian-1900)
	sec=heure*3600.d0
	t1=julian(ian,j,sec)
	return
	end
c###
      real*8  function   deg(z)   
c------------------------------------------------------------------------------
c  PROGRAMME : deg.f
c
c  DESCRIPTION : conversion from radians to deg
c
c
c  PROGRAMMEURS : F. LYARD
c 
c  DATE derniere retouche : 19/01/94
c------------------------------------------------------------------------------
        include 'common.h'          
        real*8 z    
                          
      deg=(z*1.8D+02)/pi
      return
      end
c###
      subroutine grenoble(t,rlon,rlat,tide,istat)
c------------------------------------------------------------------------------
c  PROGRAMME : grenoble.f
c
c  DESCRIPTION : driver for tidal computation. The result is returned in tide.
c               tide is an array of 
c
c  PROGRAMMEURS : J.M. MOLINES et F. LYARD
c 
c  DATE derniere retouche : 13/13/95
c------------------------------------------------------------------------------
       include 'common.h' 
       dimension tide(extended_max)        

c----------------------------------------------------------------------- 
c istat = Number of grid points for data interpolation
c istat = 0  point (rlon,rlat) is out of the gridded area
c istat = 1,2,3  point (rlon,rlat) is on a frontier of the area
c istat = 4  point (rlon,rlat) is fully in the domain
c-----------------------------------------------------------------------

	real*8 delta,t
                  
      delta=(t-t_nodal)*3.6525d+04*24.D+00

	if(abs(delta).gt.delta_max) then
	  call init_corrections(t)
	  delta=0.D+00
	endif

c-----------------------------------------------------------------------
c      sreal(nb) = real part of the nb th tide wave
c      simag(nb) = imaginary  part of the nb th tide wave
c-----------------------------------------------------------------------

	call interpolation(rlon,rlat,istat)
                     
	if (istat.ne.0) then
	  h=0.0
	 
         do nb=1,nbwave
          phi=real(freq(nb)*delta)+v0_u(nb)  
          h=h+f(nb)*(sreal(nb)*cos(phi)+simag(nb)*sin(phi))
	  tide(nb)=h
         enddo 
	else
	 do nb=1,nbwave
	  tide(nb)=spec
	 enddo
	endif
    
      return
      end
c###
      subroutine init_corrections(t0)
c------------------------------------------------------------------------------
c  PROGRAMME : init_corrections.f
c
c  DESCRIPTION : compute nodal corrections
c
c
c  PROGRAMMEURS : J.M. MOLINES et F. LYARD
c 
c  DATE derniere retouche : 19/01/94
c------------------------------------------------------------------------------
       include 'common.h'          
                           
      real*8 t0
                  
      call nodalc(t0)     
                   
      do nb=1,nbwave   
c-->    convert v0_u from degrees to radians
	  v0_u(nb)=v0_u(nb)*pi/180.
      enddo   

	t_nodal=t0

      return
      end
c###
      subroutine interpolation(xx,yy,istat)
c------------------------------------------------------------------------------
c  PROGRAMME : interpolation.f
c
c  DESCRIPTION : Perform bi-linear interpolation at point xx,yy from the gridded
c                files.
c          istat returns the number of points used for the interpolation
c
c
c  PROGRAMMEURS : J.M. MOLINES et F. LYARD
c 
c  DATE derniere retouche : 13/10/1995
c------------------------------------------------------------------------------
       include 'common.h'  
                       
	istat=0     

	if (nbwave.eq.0) return

      do nb=1,nbwave 
        sreal(nb)=0.0
        simag(nb)=0.0  
      enddo           
                       
	x=xx      
        y=yy   

	xmax=xmin+(ni-1)*dx
	if (x.lt.xmin) x=x+360.0
	if (x.gt.xmax) x=x-360.0

	i0=int((x-xmin)/dx)+1
	j0=int((y-ymin)/dy)+1 
                                              
	ptot=0.0

      do i=i0,i0+1 
	  if(i.lt.0.or.i.gt.ni) goto 100
        do j=j0,j0+1 
	    if(j.lt.1.or.j.gt.nj) goto 200
c Check if the 8 Major constituents are available:
	    do nb=1,8
	    if(wrp(i,j,nb).gt.(spec-0.1).
     .       or.wip(i,j,nb).gt.(spec-0.1)) goto 200
	    end do
c If Minor constituents are not available set them to 0.0
	    do nb=9,nbwave
	    if(wrp(i,j,nb).gt.(spec-0.1).
     .        or.wip(i,j,nb).gt.(spec-0.1)) then
	      wrp(i,j,nb) = 0.0
	      wip(i,j,nb) = 0.0
            end if
	    end do
c
	    istat=istat+1
	    xij=xmin+(i-1)*dx  
	    yij=ymin+(j-1)*dy
c this line modified for taking into account interpolation near the origin
            pds=(dx-amod(abs(x-xij),dx))*(dy-abs(y-yij))
	    ptot=ptot+pds   
          do nb=1,nbwave 
            sreal(nb)=sreal(nb)+wrp(i,j,nb)*pds    
            simag(nb)=simag(nb)+wip(i,j,nb)*pds
	    enddo
200	    continue
        enddo 
100	  continue 
	enddo

	do nb=1,nbwave
	  if(ptot.ne.0.0) then
	    sreal(nb)=sreal(nb)/ptot  
	    simag(nb)=simag(nb)/ptot  
	  else    
	    istat=0
	    sreal(nb)=spec  
	    simag(nb)=spec  
	  endif
      enddo  

      return
      end
c###
      real*8 function julian(annee,jour,seconde)
c------------------------------------------------------------------------------
c  PROGRAMME :julian.f
c
c  DESCRIPTION : This function return the elapsed time since Jan 1 1900 0:0TU
c                in julian centuries .
c
c  PROGRAMMEURS : J.M. MOLINES et F. LYARD
c 
c  DATE derniere retouche : 19/01/94
c------------------------------------------------------------------------------
c
      real*8 xj,bissex
      integer annee,jour
	real*4 seconde

        xj=dble((annee-1900.)*365)
	xj=xj+dble(jour-1)
	xj=xj+dble(seconde)/8.64d+04
        bissex=dble(int((annee-1901)/4))
        julian=(xj+bissex)/3.6525d+04  
	
	return
      end
c###
	function lnblnk(string)
c------------------------------------------------------------------------------
c  PROGRAMME : lnblnk.f
c
c  DESCRIPTION : This function is usually part of libU77 on unix system. 
c                 provided for portability
c
c  PROGRAMMEURS : J.M. MOLINES
c 
c  DATE derniere retouche : 19/01/94
c------------------------------------------------------------------------------
	character string*(*)
	ll= len(string)
	ii= index(string,' ')
	if (ll.eq.0.or.ii.eq.1) then
	  lnblnk= 1
	  return
	else if (ii.eq.0) then
	  lnblnk= ll
	  return
	end if
	do 10 k= ll, ii-1, -1
	  if (string(k:k).ne.' ') go to 20
10	  continue
20	lnblnk= k
	return
	end
c###
         subroutine nodal_a      
c------------------------------------------------------------------------------
c  PROGRAMME : nodal_a.f
c
c  DESCRIPTION :Compute nodal corrections from SCHUREMAN (1958)
c     
c  CAUTION     : indexes used in this routine are internal to the code
c                and corresponds to the !original! ondes.dat file.
c
c  PROGRAMMEURS : J.M. MOLINES
c 
c  DATE derniere retouche : 19/01/94
c------------------------------------------------------------------------------
       include 'common.h'     
     
      do 100 i=1,nbwave 
        if(num(i).eq.1) goto 11 
        if(num(i).eq.2) goto 12 
        if(num(i).eq.3) goto 13  
        if(num(i).eq.5) goto 14
        if(num(i).eq.6) goto 14
        if(num(i).eq.7) goto 14
        if(num(i).eq.8) goto 14
        if(num(i).eq.9) goto 14
        if(num(i).eq.11) goto 15 
        if(num(i).eq.12) goto 12
        if(num(i).eq.13) goto 12 
        if(num(i).eq.14) goto 16 
        if(num(i).eq.27) goto 11
11        f(i)=dsin(iang)*dcos(iang/2)**2/0.38  
          goto 120
12        f(i)=1 
          goto 120 
13        f(i)=sqrt(0.8965*dsin(2*iang)**2+0.6001*dsin(2*iang)
     &              *dcos(nu)+0.1006)          
          goto 120
14        f(i)=real(dcos(iang/2.D+00)**4/9.154D-01)
          goto 120
15        f(i)=dcos(iang/2)**4/0.9154*x1ra        
          goto 120 
16        f(i)=sqrt(19.0444*dsin(iang)**4
     &               +2.7702*dsin(iang)**2*dcos(2*nu)+0.0981) 
          goto 120
      
120     continue   
100   continue
      return
      end
c###
                  subroutine  nodal_G  
c------------------------------------------------------------------------------
c  PROGRAMME : nodal_g.f
c
c  DESCRIPTION : compute V0+u from Schureman (1958)
c
c
c  PROGRAMMEURS : J.M. MOLINES et F. LYARD
c 
c  DATE derniere retouche : 19/01/94
c------------------------------------------------------------------------------
       include 'common.h'          

      real*8 deg
                        
      n=deg(n)
      hp=deg(hp)
      s=deg(s)    
      p1=deg(p1)
      p=deg(p)
      xi=deg(xi)
      nu=deg(nu)
      nuprim=deg(nuprim)
      nusec=deg(nusec)
      r=deg(r) 

      do 140 i=1,nbwave
c-------------------o1-------------------------------------------------- 
        if(num(i).eq.1) then
          v0=tt-2*s+hp+90.
          u=2*xi-nu
          v0_u(i)=tt-2*s+hp+90.+2*xi-nu
          goto 130
        endif
c-------------------p1--------------------------------------------------
        if(num(i).eq.2) then  
          v0=tt-hp+90.
          u=0
          v0_u(i)=tt-hp+90.
          goto 130
        endif      
c-------------------k1--------------------------------------------------
        if(num(i).eq.3) then           
          v0=tt+hp-90.                                             
          u=-nuprim
          v0_u(i)=tt+hp-90.-nuprim
          goto 130 
        endif
                                   
c-------------------2n2------------------------------------------------
        if(num(i).eq.5) then
          v0=2*tt-4*s+2*hp+2*p
          u=2*xi-2*nu
          v0_u(i)=2*tt-4*s+2*hp+2*p+2*xi-2*nu
          goto 130
        endif                                                           
c-------------------mu2-------------------------------------------------
        if(num(i).eq.6) then
          v0=2*tt-4*s+4*hp
          u=2*xi-2*nu
          v0_u(i)=2*tt-4*s+4*hp+2*xi-2*nu
          goto 130 
        endif
c-------------------n2--------------------------------------------------
       if(num(i).eq.7) then
          v0=2*tt-3*s+2*hp+p
          u=2*xi-2*nu
          v0_u(i)=2*tt-3*s+2*hp+p+2*xi-2*nu
          goto 130
        endif
c-------------------nu2-------------------------------------------------
        if(num(i).eq.8) then
          v0=2*tt-3*s+4*hp-p
          u=2*xi-2*nu
          v0_u(i)=2*tt-3*s+4*hp-p+2*xi-2*nu
          goto 130
        endif     
c-------------------m2--------------------------------------------------
        if(num(i).eq.9) then
          v0=2.*tt-2*s+2*hp
          u=2*xi-2*nu
          v0_u(i)=2.*tt+real(-2.D+00*s+2.D+00*hp+2.D+00*xi-2.D+00*nu)
          goto 130
        endif
c------------------l2---------------------------------------------------
        if(num(i).eq.11) then
          v0=2*tt-s+2*hp-p+180
          u=2*xi-2*nu-r
          v0_u(i)=2*tt-s+2*hp-p+180+2*xi-2*nu-r
          goto 130
        endif
c------------------t2---------------------------------------------------
        if(num(i).eq.12) then
          v0=2*tt-hp+p1
          u=0
          v0_u(i)=2*tt-hp+p1
          goto 130
        endif
c------------------s2---------------------------------------------------
        if(num(i).eq.13) then
          v0=2*tt
          u=0
          v0_u(i)=2*tt
          goto 130
        endif
c------------------k2---------------------------------------------------
        if(num(i).eq.14) then
          v0=2*tt+2*hp
          u=-2*nusec
          v0_u(i)=2*tt+2*hp-2*nusec
          goto 130
        endif
c------------------q1---------------------------------------------------
        if(num(i).eq.27) then                 
          v0=tt-3*s+hp+p+90.
          u=2*xi-nu
          v0_u(i)=tt-3*s+hp+p+90.+2*xi-nu
          goto 130
        endif
c
130     v0_u(i)=amod(v0_u(i),360.00)    
c	  write(*,*) mod(v0,360.),u
	  

140   continue

      return
      end
c###
               subroutine nodalc(tj)
c------------------------------------------------------------------------------
c  PROGRAMME : nodalc.f
c
c  DESCRIPTION : Nodal correction calls
c
c
c  PROGRAMMEURS : J.M. MOLINES et F. LYARD
c 
c  DATE derniere retouche : 19/01/94
c------------------------------------------------------------------------------
                    
      real*8 tj
                                                                
      call astronomics(tj) 
      call nodal_a
      call nodal_G  

	return
      end
      real*8 function rad(w)      
c------------------------------------------------------------------------------
c  PROGRAMME : rad.f
c
c  DESCRIPTION : conversion de degres en radians
c
c
c
c  PROGRAMMEURS : J.M. MOLINES et F. LYARD
c 
c  DATE derniere retouche : 19/01/94
c------------------------------------------------------------------------------
       include 'common.h'          
c------------------------------------------------------------------------
c
c     cette fonction permet de convertir la mesure d'un angle exprimee
c     en degres --> radians.
c                         
c------------------------------------------------------------------------
      real*8 w    
                          
      rad=(w*pi)/1.8D+02
      return
      end
	integer function tojul (ian,imois,ijour)
c------------------------------------------------------------------------------
c  PROGRAMME : tojul.f
c
c  DESCRIPTION : renvoi une date au format aajjj quand elle est donnee au
c                format aaaa mm jj 
c
c
c
c  PROGRAMMEURS : J.M. MOLINES
c 
c  DATE derniere retouche : 19/01/94
c------------------------------------------------------------------------------
	dimension jmois(12)
c	integer tojul
	data jmois /31,28,31,30,31,30,31,31,30,31,30,31/
	if(mod(ian,4).eq.0) then
	jmois(2)=29
	else
	jmois(2)=28
	endif
	if(ijour.gt.jmois(imois)) then
	print *,' ERROR in tojul : There are only ',jmois(imois),
     |  'days in month number ',imois
	stop 'FATAL ERROR'
	endif
	ia=ian-1900
	s=0
	do m=1,imois-1
	s=s+jmois(m)
	enddo
	ida=s+ijour
	tojul=ia*1000+ida
	return
	end
c###
      subroutine init_data
c------------------------------------------------------------------------------
c  PROGRAMME : init_data
c
c  DESCRIPTION : Read the data files for each tidal constituents. The path
c               to the directory where the data files are located is given
c               by an environment variable WAVE_PATH, that must be setenv
c               before running the program.
c              (e.g :setenv WAVE_PATH /users/toto/Grenoble_DATA ).
c               nbext is the number of available constituents. Actually 13.
c              All files are named {wave_name}.model (ASCII files) or
c                                  {wave_name}.model.bimg (BINARY files).
c                       and model is a valid extension for the files
c
c  PROGRAMMEURS : J.M. MOLINES et F. LYARD
c 
c  DATE derniere retouche : 13/03/1995
c------------------------------------------------------------------------------
       include 'common.h'  

      character*255    pathname    

      character*6      name(13),extended(extended_max)
      character*3 rep
      character*255 WAVE_PATH  
      real*8 frequency(13)
      integer*4        code(13),unit2 ,lex,lmod 
c
c these values are internal to the program and to the nodal correction
c routines. Do not change under any circumstances.
	data (code(i),name(i),frequency(i),i=1,13) /
     | 27,'q1',13.39866087990d0,
     | 1,'o1' ,13.94303558000d0,
     | 2,'p1' ,14.95893136000d0,
     | 3,'k1' ,15.04106864000d0,
     | 5,'2n2',27.89535481990d0,
     | 6,'mu2',27.96820844000d0,
     | 7, 'n2',28.43972952010d0,
     | 8,'nu2',28.51258314000d0,
     | 9, 'm2',28.98410422000d0,
     | 11,'l2',29.52847892000d0,
     | 12,'t2',29.95893332010d0,
     | 13,'s2',30.00000000000d0,
     | 14,'k2',30.08213728000d0/
c
	call getenv('WAVE_PATH',WAVE_PATH)
	lwa=lnblnk(WAVE_PATH)
	WAVE_PATH(lwa+1:lwa+1)='/'
	lwa=lwa+1
	print *, 'Initialisation ...please be patient ...'
	print *,'Used environment variable :'
	print *,'WAVE_PATH=',WAVE_PATH(1:lwa)
#ifndef BIMG_FILE 
	print *,'Reading ASCII files'
#else
	print *,'Reading BIMG files'
#endif
c
	nbext=13
c
	extended(1)='m2'               
	extended(2)='s2'               
	extended(3)='k1'               
	extended(4)='o1'               
	extended(5)='n2' 
	extended(6)='p1'               
	extended(7)='k2'               
	extended(8)='q1'               
c   
	extended(9)='2n2'               
	extended(10)='mu2'               
	extended(11)='nu2'               
	extended(12)='l2'               
	extended(13)='t2' 

	unit2=20
             
	nb=0      
742	continue
      do nn=1,nbext 

	  nb=nb+1 
	  lex=lnblnk(extended(nn))
	  lmod=lnblnk(model)
#ifndef BIMG_FILE
	  pathname=WAVE_PATH(1:lwa)//
     |     extended(nn)(1:lex)//'.'//model(1:lmod)
#else
	  pathname=WAVE_PATH(1:lwa)//
     |     extended(nn)(1:lex)//'.'//model(1:lmod)//'.bimg'
#endif
c
c
120	  wave(nb)=extended(nn)
	if (rep.eq.'N') goto 101 
	print *,'trying to read ',pathname(1:lnblnk(pathname))
        call read_tide(nb,unit2,pathname,istat) 

	  if (istat.eq.0) then 
	    nbwave=nb
          do i=1,13         
            if(wave(nb).eq.name(i)) then
              freq(nb)=frequency(i)
              num(nb)=code(i)   
c              write(*,*) wave(nb),real(freq(nb)),' degree/hour'
	        freq(nb)=freq(nb)*pi/180.d+00
              goto 100
            endif
	    enddo 
	    stop 'unknown wave...'
100	    continue
	  else
	    nb=nb-1
	  endif
101	continue
	enddo  

10	format(i3,3x,a6,3x,d14.11) 
             
                 
	call init_corrections(0.0d0)
	print *, 'Now ready !!!'
      return
      end
c###
      subroutine read_tide(nb,iunit,pathname,istat)
c------------------------------------------------------------------------------
c  PROGRAMME :read_tide.f
c
c  DESCRIPTION : Read data files either  ASCII formatted or BIMG files
c                (Check the CPPFLAG in Makefile)
c                to convert ASCII to BIMG use ascii2bimg.f program
c
c
c
c  PROGRAMMEURS :  F. LYARD
c 
c  DATE derniere retouche : 13/03/1995
c------------------------------------------------------------------------------
       include 'common.h'   

	integer iunit
	character*255 pathname
	character*80 line1,line2,line3,line4
                    
	real a(30),G(30)  

	conv=pi/180. 

	spec=32767.

	istat=0
c
c default files are ASCII files, (BIMG_FILE not defined)
c
#ifndef BIMG_FILE    
	 open(iunit,file=pathname,status='old',err=100)
	read(iunit,*) xmin,xma,line1
	read(iunit,*) ymin,yma
	read(iunit,*) dx,dy
	read(iunit,*)   ni,nj
	read(iunit,*) amask,Gmask
	print *,line1

	do j=1,nj
	  do i=1,ni,30
	    read(iunit,10) (a(k),k=1,30)
	    read(iunit,10) (G(k),k=1,30)
	    do k=1,30
	      ii=k+(i-1)
	      if (a(k).ne.amask.and.G(k).ne.Gmask) then
	        wrp(ii,j,nb)=a(k)*cos(G(k)*conv)
	        wip(ii,j,nb)=a(k)*sin(G(k)*conv)
	      else
	        wrp(ii,j,nb)=spec
	        wip(ii,j,nb)=spec
	      endif
	    enddo
	  enddo

c
c tide(-0.5,phi)=tide(359.5,phi)
c convenient for interpolation
c over a global solution
c
          wrp(0,j,nb)=wrp(ni,j,nb)
	  wip(0,j,nb)=wip(ni,j,nb)
          wrp(ni+1,j,nb)=wrp(1,j,nb)
          wip(ni+1,j,nb)=wip(1,j,nb)

	enddo	    
10	format(30f7.2)
20	format(30f7.1)

	close(iunit) 
#else   
	open(iunit,file=pathname,status='old',form='unformatted',
     |       err=100)
	read(iunit)line1
	read(iunit)line2
	read(iunit)line3
	read(iunit)line4
	read(iunit)nimax,njmax,nz,nt,ndim,icod
	print *, line1
Cjw970415	print *, line2
Cjw970415	print *, line3
Cjw970415	print *, line4
	ni=nimax-1
	nj=njmax
	read(iunit)xmin,ymin,dx,dy,amask
	read(iunit) z0
	read(iunit) tim_tag
	read(iunit) ((wrp(i,j,nb),i=0,ni),j=1,nj)
	read(iunit) ((wip(i,j,nb),i=0,ni),j=1,nj)
c row ni+1 is filled with value at i=1 for interpolation
        do j=1,nj
          wrp(ni+1,j,nb)=wrp(1,j,nb)
          wip(ni+1,j,nb)=wip(1,j,nb)
        enddo
	 do j=1,nj
	  do i=0,ni+1
	   a(1)=wrp(i,j,nb)
	   G(1)=wip(i,j,nb)
	      if (a(1).ne.amask.and.G(1).ne.amask) then
	        wrp(i,j,nb)=a(1)*cos(G(1)*conv)
	        wip(i,j,nb)=a(1)*sin(G(1)*conv)
	      else
	        wrp(i,j,nb)=spec
	        wip(i,j,nb)=spec
	      endif
	  enddo
	 enddo	
	close(iunit)
#endif
    
	ni=ni+1

      print *,'  Wave ',wave(nb),' ok...' 
	return
            
100	continue
       print *,' ERROR opening data file', pathname
	istat=-1
	stop 'FATAL ERROR'

      return
      end
