C
C comturb.h created Mon Jun  8 12:18:54 MET DST 1998 
C
C Created using /home/kbk/scr/newdoc by kbk
C


C
C
C $Id: comturb.h,v 1.1 2002-05-02 14:01:56 gotm Exp $
C
C

C
C
C $Log: comturb.h,v $
C Revision 1.1  2002-05-02 14:01:56  gotm
C Initial revision
C
C Revision 1.1.1.1  2001/04/17 08:43:08  bbh
C initial import into CVS
C
C
C
      double precision cde,sige        ! Parameters for k-epsilon model 


      common /comturb/ 
     &		cde,sige


      double precision cmue,tkemin,epsmin,sigk,ce1,ce2,ce3,z0s 

      parameter(cmue=0.5625d0) 
      parameter(tkemin=1.0d-7) 
      parameter(epsmin=5.0d-10) 
      parameter(sigk=1.0) 
      parameter(ce1=1.44) 
      parameter(ce2=1.92) 
      parameter(ce3=-0.9) 
      parameter(z0s=0.1)

