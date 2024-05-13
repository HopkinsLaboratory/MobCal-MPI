C     PROGRAM MOBCAL_MPI v2.0.1
C NEW ENTRY
C     Mobility Calculation program published in
C       Haack, A.; Ieritano, C.; Hopkins, W. S. MobCal-MPI 2.0:
C       An Accurate and Parallelized Package for Calculating Field-
C       Dependent Collision Cross Sections and Ion Mobilities.
C       Analyst 2023, 148 (14), 3257–3273. DOI: 10.1039/D3AN00545C.
C
C**************************************************************************
C**************************************************************************
C     Adapted from original code with reference:
c     PROGRAM MOBCAL
c
c     Program to Calculate Mobilities
c
c     See: M. F. Mesleh, J. M. Hunter, A. A. Shvartsburg, G. C. Schatz,
c     and M. F. Jarrold, "Structural Information from Ion Mobility
c     Measurements: Effects of the Long Range Potential" J. Phys. Chem.
c     1996, 100, 16082-16086; A. A. Shvartsburg and M. F. Jarrold,
c     "An Exact Hard Spheres Scattering Model for the Mobilities of
c     Polyatomic Ions" Chem. Phys. Lett. 1996, 261, 86.
c
c     Version of 06/15/00 junko.f
c     Modified to include Na+ parameters and how masses are read
c     Version of 09/14/98
c     Corrected center of mass calculation in fcoord and ncoord and
c     changed how structural asymmetry parameter is calculated
c     Version of 08/15/98
c     Modified to calculate structural asymmetry parameter
c     Version of 02/20/98
c     Modified GSANG: removed bug that caused trajectory start position
c     to go into endless loop for long structures. Also fixed TESTXRAND
c     and summary print out.
c     Version of 01/29/98
c     Modified MOBIL2: removed a bug that messed up b2max calculation
c     for large molecules
c     Version of 12/24/97
c     Modified to permit averaging over multiple conformations 
c     Version of 12/08/97
c     Extensively modified to allow for a molecule consisting of
c     an unlimited variety of different atoms. Exact hard sphere
c     scattering incorporated. RANLUX random number generator 
c     incorporated. Several obsolete subroutines removed.
c     Version of 10/16/97
c     Modified to allow uniform and non-uniform charge distributions
c     Version of 08/23/97
c     Modified to allow calculations with a non-uniform charge distribution
c     Version of 21/10/19
C
C     Version 17/01/21
C     Adding multiple temp runs in same input file
C     Version 14/04/21
C     Added calculation of omega*(1,4) and omega*(2,3)
C     Version 21/06/21
C     velocity integration is now over a linear grid. this makes the
C     higher order integrals accurate as well.
C     Also added omega*(2,4), omega*(3,3) and omega*(3,4)
C     modified the output printed: velocity grid, omega*(l,s) summary
C
C     Version 01/05/23
C     implemented 2nd order two-temperature theory
C     updated the output significantly
C     updated the error analysis
C     for details, see Analyst 2023, 148 (14), 3257–3273.
C     
C***  Explicit seed and correction to calculation of ion mobility (K); does not affect previous CCS calculations
c
c     MOBIL2 calculates the mobility using a Lennard-Jones plus
c     ion-induced dipole potential. MOBIL2 uses Monte Carlo integrations
c     over orientation and impact parameter and a numerical integration
c     over g* (the reduced velocity). MOBIL2 includes second order
c     corrections, though these don't appear to be important.
c
c     GSANG/DERIV/DIFFEQ are subroutines derived from code provided by
c     George Schatz. They use Runge-Kutta-Gill and Adams-Moulton predictor
c     -corrector integration methods. These subroutines are now set-up to
c     work in three dimensions.
c
c     DLJPOT calculates the potential and the derivatives of the potential.
c     The potential is given by a sum of 6-12 two body interactions and
c     ion-induced dipole interactions. The charge can be distributed over 
c     all atoms equally, a specific charge distribution can be employed, or
c     the charge can be set to zero (only the 6-12 part of the potential
c     is used). Each atom in the cluster or molecule can have different 
c     Lennard-Jones parameters.
c     
c     MOBIL4 calculates the exact hard sphere scattering mobility and
c     the projection approximation mobility. Adapted from code written
c     by Alexandre Shvartsburg.
c
c     XRAND/RANLUX are random number generators. RANLUX is the standard
c     ranlux number generator of M. Luscher (code from F. James). XRAND
c     is a function that calls either RANLUX or RAND (the standard F77
c     random number generator).
c
c     FCOORD reads in the coordinates, integer masses, and partial charges.  
c
c     RANTATE/ROTATE rotates the coordinates of the cluster.  
c
c     TRAJ calculates a series of trajectories for closer examination.
c
c     TRAJONE calculates one trajectory for a given velocity, impact
c     parameter, and angles (theta, phi, and gamma in degrees).
c
c     POTENT determines the average potential. (Only for near spherical
c     molecules or clusters).
c     
c
c     ***************************************************************
C******************************************************************************
C******************************************************************************
C     Major changes include inclusion of "softer" exp-6 Lennard-Jones
C     potential and more exhaustive paramater set which needs to be
C     included as input
c
c
      use mpi
      implicit double precision (a-h,m-z)
      parameter(inpmax=201)
      dimension tmc(100),
     ?asympp(100)
      character*50 filen1,filen2,unit,dchar,xlabel,infile
      parameter (ixlen=1000)
      character*2 flend(32)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xeo,xk,xn,xmv,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(ixlen),k_const(ixlen),tcharge
      common/coordinates/fx(ixlen),fy(ixlen),fz(ixlen),
     ?ox(ixlen),oy(ixlen),oz(ixlen)
      common/ljparameters/eolj(ixlen),rolj(ixlen),eox4(ixlen),
     ?ro6lj(ixlen),ro12lj(ixlen),dro6(ixlen),dro12(ixlen)
      common/ljparameters2/alphai(ixlen),Ni(ixlen),Ai(ixlen),Gi(ixlen),
     1 RijStar(ixlen),eij(ixlen)
      common/mmffN2/alphaN2,NiN2,AiN2,GiN2,RN2Star,MMFF_B,MMFF_beta
      common/mmffHe/alphaHe,NiHe,AiHe,GiHe,RHeStar
      common/hsparameters/rhs(ixlen),rhs2(ixlen)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/temperature/tt1,tt2,ittgrid
      common/xrandom/i1,i2,i3,i4,i5,i6
      common/scaling/enersca,distsca
      common/empcorr/Aec,Bec
C****NEED COMMON BLOCK FOR MPI ELEMENTS
      common/mobilcal/b2max(inpmax),parab2max(inpmax)
      common/mobilgrid/pgst(inpmax),wgst(inpmax),q1st(inpmax),
     1  q2st(inpmax),q3st(inpmax),q1STD(inpmax),q2STD(inpmax),
     1  q3STD(inpmax)
      common/mpicon/imyrank,inprocs,imp_per_node,inp_per_node,s_time
      common/runparams/itn,inp,imp,igas
      dimension tmmarray(100),csarray(100),ENarray(100),sdevarray(100)
C******
      data flend/'_1','_2','_3','_4','_5','_6','_7','_8',
     1 '_9','10','11','12','13','14','15','16','17','18',
     1 '19','20','21','22','23','24','25','26','27','28',
     1 '29','30','31','32'/
C****Start MPI Stuff
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,inprocs,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,imyrank,ierr)
      s_time=MPI_WTIME()
c
C***Read input file and output file name from command prompt
      if(imyrank.eq.0)then
       call get_command_argument(1,filen1)
	   call get_command_argument(2,filen2)
      endif
c
C************
c     Define parameters for MOBIL2
c
c     Number of complete cycles for average mobility calculation in 
c     MOBIL2. Default value is 10.
c
      itn=10
c
c     Number of points in velocity integration in MOBIL2. Default 
c     value is 20.
c
      inp=40
c
c     Number of points in Monte Carlo integrations of impact parameter
c     and orientation in MOBIL2. Default value is 500.
c
      imp=500
      igas=2
C****Only node 0 opens output and input file
      if(imyrank.eq.0)then
       open(unit=8,file=filen2)
      else
       open(unit=1000+imyrank,file='output'//flend(imyrank)//'.out',
     1  status='unknown')
      endif
      
C******
c
c     print switches ip=1  print scattering angles
c                    it=1  print trajectory
c                    iu1=1 print initial coordinates from FCOORD
c                    iu2=1 print angles and coordinates from ROTATE 
c                    iu3=1 print angles from ROTATE
c                    iv=1  print all potentials from POTENT
c                    im2=1 print brief information from MOBIL2
c                    im4=1 print brief information from MOBIL4
c                    igs=1 print out angles, v, and b from MOBIL4 into 
c                          a temporary file called hold
c
      ip=0
      it=0
      iu1=0
      iu2=0
      iu3=0
      iv=0
      im2=0
      im4=0
      igs=0
c
c     Constants from Handbook of Chemistry and Physics, 70th Edition
c     xmv yields reduced mobility K0
c
      pi=3.14159265358979323846d0
      cang=180.d0/pi
      xe=1.60217733d-19
      xk=1.380658d-23
      xn=6.0221367d23
      xeo=8.854187817d-12
      xmv=8.20573660809596d-5*273.15
c
C*** Define MMFF parameters for Diatomic Nitrogen
      alphaN2=1.739688487d0
      NiN2=5.918d0
      AiN2=3.1781544d0
      GiN2=1.16175013d0
      RN2Star=AiN2*(alphaN2**(0.25d0))
C*** Define MMFF parameters for Helium
      alphaHe=0.205236d0
      NiHe=1.42d0
      AiHe=4.40d0
      GiHe=1.209d0
      RHeStar=AiHe*(alphaHe**(0.25d0))
C****************
      MMFF_B=0.2d0
      MMFF_beta=12.0d0
C*****
      if(imyrank.eq.0)then
       call fcoord(filen1,unit,dchar,xlabel,asympp(1))
       close(9)
      endif
C*******************
C****Broadcast******
      call MPI_BCAST(itn,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(i2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(inp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(imp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(igas,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(inatom,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(icoord,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(m2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tt1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tt2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ittgrid,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(romax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(fx,ixlen,MPI_DOUBLE_PRECISION,0,
     1 MPI_COMM_WORLD,ierr)
      call MPI_BCAST(fy,ixlen,MPI_DOUBLE_PRECISION,0,
     1 MPI_COMM_WORLD,ierr)
      call MPI_BCAST(fz,ixlen,MPI_DOUBLE_PRECISION,0,
     1 MPI_COMM_WORLD,ierr)
      call MPI_BCAST(pcharge,ixlen,MPI_DOUBLE_PRECISION,0,
     1 MPI_COMM_WORLD,ierr)
      call MPI_BCAST(eij,ixlen,MPI_DOUBLE_PRECISION,0,
     1 MPI_COMM_WORLD,ierr)
      call MPI_BCAST(RijStar,ixlen,MPI_DOUBLE_PRECISION,0,
     1 MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ox,ixlen,MPI_DOUBLE_PRECISION,0,
     1 MPI_COMM_WORLD,ierr)
      call MPI_BCAST(oy,ixlen,MPI_DOUBLE_PRECISION,0,
     1 MPI_COMM_WORLD,ierr)
      call MPI_BCAST(oz,ixlen,MPI_DOUBLE_PRECISION,0,
     1 MPI_COMM_WORLD,ierr)
C******k*************
C***Added loop here over temperature array
C*******************
c     Temperature
      if(imyrank.eq.0)then
       write(8,*)'Temperature range from input read:'
       write(8,'(A,f7.2,A)') '  T_bath    = ', tt1, 'K'
       write(8,'(A,f7.2,A)') '  T_eff max = ', tt2, 'K'
       write(8,'(A,i4)')     '  # of T_eff grid points =', ittgrid
      endif
      eo=1.34d-03*xe
      ro=3.043d0*1.0d-10
      ro2=ro*ro
      if(imyrank.eq.0)write(8,600) eo/xe,ro/1.0d-10
      if(imyrank.eq.0)then
       write(8,'(A,f6.3,A,f6.3)')'  rho_e=',enersca,' rho_s=',distsca
      endif
c     empiciral correction to 2TT mobility
      if(correct.gt.0)then
       Aec=6.11d-2
       Bec=143.0d0
       if(imyrank.eq.0)then
        write(8,*)'Empirical correction: 1+A*exp(-B/(E/N))'
        write(8,'(A,f8.4)')  '  A = ', Aec
        write(8,'(A,f6.2,A)')'  B = ', Bec, ' Td'
       endif
      else
       Aec=0.0d0
       Bec=1.0d0
       if(imyrank.eq.0)then
        write(8,*)'Empirical correction turned off.'
       endif
      endif
c
c     Constant for ion-induced dipole potential
c
c     Polarizability of helium = 0.204956d-30 m3
c     xeo is permitivity of vacuum, 8.854187817d-12 F.m-1
c
      if(igas.eq.2)then
       dipol=1.740d-30/(2.d0*4.d0*pi*xeo)
       dipol=dipol*xe*xe 
       if(imyrank.eq.0)write(8,601) dipol
c
c     Mass constants
c
       m1=28.0134d0
       mu=((m1*m2)/(m1+m2))/(xn*1.0d3)
      elseif(igas.eq.1)then
       dipol=0.204956d-30/(2.0d0*4.0d0*pi*xeo)
       dipol=dipol*xe*xe
       if(imyrank.eq.0)write(8,601) dipol
c
c     Mass constants
c
       m1=4.0026d0
       mu=((m1*m2)/(m1+m2))/(xn*1.0d3)
      endif
c
C****PATCH PATCH
C**** Some depreciated print statements. Remove?
C      if(imyrank.eq.0)write(8,602) mconst
c
c
C      if(imyrank.eq.0)write(8,603) t
C****PATCH PATCH
c
c     Define parameters for random number generator
c
c     If i5=1 RANLUX is used otherwise RAND is used. If RAND is used 
c     i2 is the seed integer. If RANLUX is used i2, i3, and i4 are seed
c     integers. i3 and i4 (which are used to start RANLUX in a particular
c     place) are usually set to zero. i1 contain the "luxury level" - how
c     good the random number generator is. Values of 0 to 4 can be used
c     and the default value in RANLUX is 3. At this level "any 
c     theoretically possible correlations have a very small chance of
c     being observed".
c
      i1=3
      i3=0
      i4=0
      i5=1
      if(i5.ne.1) then 
      if(imyrank.eq.0)write(8,604) i2
      else
      if(imyrank.eq.0)write(8,619) 
      endif
c     initialize the random number generators
C****Broadcast i2 and initialize several times to mix things up...
      call MPI_BCAST(i2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      do i=1,imyrank+1
       call rluxgo(i1,i2,i3,i4)
       call srand(i2)
       i2=int(xrand()*1d8+xrand()*1d7+
     1  xrand()*1d6+xrand()*1d5+
     1  xrand()*1d4+xrand()*1d3+
     1  xrand()*1d2)
      enddo
C************
C***Calc inp per node
      if(mod(inp,inprocs).ne.0)then
       iadd=inprocs-mod(inp,inprocs)
       inp=inp+iadd
      endif
      inp_per_node=inp/inprocs
C*********************    
C***Test imp to make sure it fits evenly over inprocs
      if(mod(imp,inprocs).ne.0)then
       iadd=inprocs-mod(imp,inprocs)
       imp=imp+iadd
      endif
      imp_per_node=imp/inprocs
 
c      imp=500
c
c     Minimum value of (1-cosX). This quantity determines the maximum
c     impact parameter at each velocity. Default value is 0.0005.
c
      cmin=0.0005d0
c
c     Define some parameters for trajectories: sw1 defines the potential
c     energy where the trajectory starts and dtsf1 is related to the 
c     time step at the start of the trajectory. As the trajectory comes
c     close to a collision the time step is reduced. sw2 defines the 
c     potential energy where the time step is reduced and dtsf2 defines
c     the reduced time step. Default values are: 
c     sw1 = 0.00005   dtsf1=0.5
c     sw2 = 0.0025    dtsf2=0.05
c     inwr is the number of integration steps before the program tests
c     to see if the trajectory is done or lost. ifail is the number of
c     failed trajectories that are permitted (a failed trajectory is one
c     that does not conserve energy to within 1%. Default values are: 
c     inwr = 1        ifail = 100
c
      sw1=0.00005d0
      sw2=0.005d0 
      dtsf1=0.5d0
      dtsf2=0.1d0
      inwr=1
      ifail=100
      ifailc=0
C******
      tmm=0.d0
      cs=0.0d0
C***********************************************************
C***********************************************************
      call pre_mobil2()
      deltt=0.0d0
      if(ittgrid.gt.1)then
       deltt=(tt2-tt1)/dble(ittgrid)
      endif
C*** Only need to call mobil2 for irank=0
      if(imyrank.eq.0)then
C*** Print out Header
       write(8,777)ittgrid
       write(8,751)
       do itloop=1,ittgrid+1
        ttin=tt1+dble(itloop-1)*deltt
C*** Moved these from above inside the T loop
c     Mobility constant
c
        mconst=dsqrt(18.d0*pi)/16.d0
        mconst=mconst*dsqrt(xn*1.0d3)*dsqrt((1.d0/m1)+(1.d0/m2))
        mconst=mconst*xe*dabs(tcharge)/dsqrt(xk)
        dens=xn/xmv
        mconst=mconst/dens
C****
        call mobil2(ttin,tmm,cs,sdevpc,EN)
        tmmarray(itloop)=tmm*1.d4
        ENarray(itloop)=EN*1.d21
        csarray(itloop)=cs*1.d20
        sdevarray(itloop)=sdevpc
       enddo
      endif
C***********************************************************
C***********************************************************
c     print out summary
C********************************
C***Only print out for 0 node
C********************************
      if(imyrank.eq.0)then
c
       write(8,605) filen1,xlabel
       write(8,615)itn,inp,imp,itn*inp*imp,ifailc
       if(igas.eq.2)then
        write(8,*)'Mobility Calculated under N2 gas'
       elseif(igas.eq.1)then
        write(8,*)'Mobility Calculated under He gas'
       endif
C*       write(8,*) 'enersca=',enersca,' distsca=',distsca
       e_time=MPI_WTIME()
       write(8,*)'Job Completed in',e_time-s_time,'s'
       call flush(8)
      endif
C****End temperature loop here
C      enddo
C****Collective output:
      write(8,*)' '
      write(8,'(22("*"),"Mobility Summary",22("*"))')
      write(8,'(t1,21("*"),"at T_bath=",f7.2,"K",21("*"))')tt1
      write(8,'(t1,a)')' Teff [K]   E/N [Td]  K0 [cm**2/Vs]  CCS [A**2]
     1  uncertainty'
      write(8,'(t1,60("-"))')
      deltt=0.0d0
      if(ittgrid.gt.1)then
       deltt=(tt2-tt1)/dble(ittgrid)
      endif
      do itloop=1,ittgrid+1
       ttin=tt1+dble(itloop-1)*deltt
       if(imyrank.eq.0)then
        write(8,'(t1,1x,f7.2,4x,f7.2,6x,f7.4,5x,f8.2,5x,f6.2,"%")')
     1   ttin,ENarray(itloop),tmmarray(itloop),csarray(itloop),
     1   sdevarray(itloop)
       endif
      enddo
      write(8,'(t1,60("*"))')
      write(8,*)'   '
C********************************
C***Only print out for 0 node
C********************************
      if(imyrank.eq.0)then
       close(8)
      else
       close(1000+imyrank)
      endif
C********************************
C********************************
      call MPI_FINALIZE(ierr)
C****
601   format(1x,'dipole constant =',1pe11.4)
602   format(1x,'mobility constant =',1pe11.4)
603   format(1x,'temperature =',1pe11.4)
619   format(/1x,'using RANLUX',/)
604   format(1x,'using RAND with seed integer =',i8)
635   format(/)
621   format(/1x,'coordinate set =',i5)
637   format(/1x,'structural asymmetry parameter =',f8.4)
608   format(1x,'using no charge - only LJ interactions')
612   format(/1x,'inverse average EHS mobility =',1pe11.4,
     1 /1x,'average EHS cross section =',e11.4)
611   format(/1x,'inverse average PA mobility =',1pe11.4,
     1 /1x,'average PA cross section =',e11.4)
609   format(/1x,'mobility calculation by MOBIL4 (HS scattering)')
610   format(/1x,'number of Monte Carlo trajectories =',i7,
     1 /1x,'maximum number of reflections encountered =',i3)
616   format(1x,'temperature =',1pe11.4)
617   format(1x,'using RAND with seed integer =',i8)
620   format(1x,'using RANLUX with seed integer =',i8)
636   format(1x,'structural asymmetry parameter =',1pe11.4) 
607   format(1x,'using a calculated (non-uniform) charge distribution')
606   format(1x,'using a uniform charge distribution')
605   format(///1x,'SUMMARY',//1x,'program version =
     1 MobCal-MPI_2.f',
     1 /1x,'input file name = ',a50,/1x,'input file label = ',a50)
615   format(//1x,'number of complete cycles (itn) = ',i7,/1x,
     1 'number of velocity points (inp) = ',i7,/1x,
     1 'number of random points   (imp) = ',i7,/1x,
     1 'total number of points          = ',i7,/1x,
     1 'number of failed trajectories   = ',i7)
614   format(/1x,'trajectory parameters',/1x,'sw1 =',1pe11.4,7x,
     1 'sw2 =',e11.4,/1x,'dtsf1 =',e11.4,5x,'dtsf2 =',e11.4)
613   format(/1x,'calculating trajectories to obtain momentum transfer
     1   cross sections')
634   format(/1x,'minimum and maximum number of reflections =',
     1 i3,2x,i3)
622   format(//4x,'set',5x,'PA CS',7x,'PA MOB^-1',6x,'EHSS CS',
     1 6x,'EHSS MOB^-1',4x,'ASYMP')  
626   format(/1x,'number of Monte Carlo trajectories =',i7)
625   format(/1x,'mobility calculation by MOBIL4 (HS scattering)')
623   format(1x,i5,3x,1pe11.4,3x,e11.4,3x,e11.4,3x,e11.4,3x,f8.4) 
624   format(/3x,'AVGE',2x,1pe11.4,3x,e11.4,3x,e11.4,3x,e11.4,3x,f8.4)
628   format(/1x,'trajectory parameters',/1x,'sw1 =',1pe11.4,7x,
     1 'sw2 =',e11.4,/1x,'dtsf1 =',e11.4,5x,'dtsf2 =',e11.4)
631   format(1x,i5,3x,1pe11.4,3x,e11.4)
629   format(//1x,'number of complete cycles (itn) = ',i7,/1x,
     1 'number of velocity points (inp) = ',i7,/1x,
     1 'number of random points   (imp) = ',i7,/1x,
     1 'total number of points          = ',i7,/)
630   format(//4x,'set',5x,'TM CS',7x,'TM MOB^-1')
633   format(/1x,'total number of failed trajectories =',i4)
632   format(/3x,'AVGE',2x,1pe11.4,3x,e11.4)
627   format(/1x,'mobility calculation by MOBIL2 (trajectory method)')
600   format(1x,'van der Waals scaling parameters: eo=',
     1  1pe11.4,'eV, ro=',e11.4,'A')
777   format(//1x,61('-'),/1x,'Starting CCS calculations on ittgrid = '
     1  ,i3,' temperature points',/1x,61('-'))
750   format(1x,a7,1x,f5.2)
751   format(/1x,'Weighting functions w_s(gst)',/1x,
     1 'w1(gst) = 1/( 1*Tst**3) * gst**05 * exp(-gst**2/Tst)',/1x,
     1 'w2(gst) = 1/( 3*Tst**4) * gst**07 * exp(-gst**2/Tst)',/1x,
     1 'w3(gst) = 1/(12*Tst**5) * gst**09 * exp(-gst**2/Tst)',/1x,
     1 'w4(gst) = 1/(60*Tst**6) * gst**11 * exp(-gst**2/Tst)')
      stop
      end
c
c     ***************************************************************
c
      subroutine fcoord(filen1,unit,dchar,xlabel,asymp)
c
c     Reads in coordinates and other parameters.
c
      use mpi
      implicit double precision (a-h,m-z)
      character*50 filen1,unit,dchar,xlabel
      parameter (ixlen=1000,inpmax=201)
      dimension imass(ixlen),xmass(ixlen) 
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xeo,xk,xn,xmv,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(ixlen),k_const(ixlen),tcharge
      common/coordinates/fx(ixlen),fy(ixlen),fz(ixlen),
     ?ox(ixlen),oy(ixlen),oz(ixlen)
      common/ljparameters/eolj(ixlen),rolj(ixlen),eox4(ixlen),
     ?ro6lj(ixlen),ro12lj(ixlen),dro6(ixlen),dro12(ixlen)
      common/ljparameters2/alphai(ixlen),Ni(ixlen),Ai(ixlen),Gi(ixlen),
     1 RijStar(ixlen),eij(ixlen)
      common/mmffN2/alphaN2,NiN2,AiN2,GiN2,RN2Star,MMFF_B,MMFF_beta
      common/mmffHe/alphaHe,NiHe,AiHe,GiHe,RHeStar
      common/hsparameters/rhs(ixlen),rhs2(ixlen)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/temperature/tt1,tt2,ittgrid
      common/xrandom/i1,i2,i3,i4,i5,i6
      common/scaling/enersca,distsca
      common/mobilcal/b2max(inpmax),parab2max(inpmax)
      common/mpicon/imyrank,inprocs,imp_per_node,inp_per_node,s_time
      common/runparams/itn,inp,imp,igas
      character*1000 readline
C*** Initialize temperatures in case full swath not provided
      tt1=-999.0d0
      tt2=-999.0d0
      ittgrid=-999
c
      write(8,603) filen1
      open (9,file=filen1)
      read(9,'(a30)') xlabel
      write(8,601) xlabel
      read(9,*) icoord
      write(8,650) icoord
c 
      read(9,*) inatom
      write(8,612) inatom
      read(9,'(a30)') unit
      read(9,'(a30)') dchar
c
      if(unit.eq.'au') write(8,611)
      if(unit.eq.'ang') write(8,614)
      if(unit.ne.'au'.and.unit.ne.'ang') then
       write(8,610)
       close(8)
       call MPI_FINALIZE(ierr)
       stop
      endif
c
      read(9,*) correct
c      write(8,613) correct
c
      read(9,'(t1,a1000)')readline
      ilen=len(trim(readline))+1
      isw=0
      k=0
      do i=1,ilen
       if(isw.eq.1)then
        if(readline(i:i).eq.' '.or.i.eq.ilen)then
         ifin=i-1
         if(i.eq.ilen)ifin=i
         isw=0
         k=k+1
         if(k.eq.1)then
          read(readline(istart:ifin),*)itn
         elseif(k.eq.2)then
          read(readline(istart:ifin),*)inp
         elseif(k.eq.3)then
          read(readline(istart:ifin),*)imp
         elseif(k.eq.4)then
          read(readline(istart:ifin),*)igas
         elseif(k.eq.5)then
          read(readline(istart:ifin),*)i2
         elseif(k.eq.6)then
          read(readline(istart:ifin),*)tt1
         elseif(k.eq.7)then
          read(readline(istart:ifin),*)tt2
         elseif(k.eq.8)then
          read(readline(istart:ifin),*)ittgrid
         endif
        endif
       endif
       if(isw.eq.0)then
        if(readline(i:i).ne.' ')then
         istart=i
         isw=1
        endif
       endif
      enddo
      if(tt1.lt.0.0d0)then
       write(8,616)
       call MPI_FINALIZE(ierr)
       close(8)
       stop
      elseif(tt2.lt.0.0d0)then
       tt2=tt1
       ittgrid=0
      endif
      if(ittgrid.lt.0)then
       write(8,617)
       ittgrid=10
      endif
C****node 0 can read in coord. then broadcast
      if(igas.eq.2)then
       enersca=1.275d0
       distsca=0.825d0
      else
C** updated Oct 2023
       enersca=0.93d0
       distsca=0.97d0
c       enersca=0.81d0
c       distsca=0.98d0
      endif
c
      if(dchar.eq.'equal') write(8,630) 
      if(dchar.eq.'calc') write(8,631)
      if(dchar.eq.'none') write(8,633)
      if(dchar.ne.'equal'.and.dchar.ne.'calc'.and.dchar
     1 .ne.'none') then
      write (8,632)
      close(8)
      call MPI_FINALIZE(ierr)
      stop
      endif
c
      tcharge=0.0d0
      acharge=0.0d0
      do 2000 iatom=1,inatom
C****Read in vdw parameters from input
       read(9,*) fx(iatom),fy(iatom),fz(iatom),
     1  xmass(iatom),pcharge(iatom),alphai(iatom),
     1  Ni(iatom),Ai(iatom),Gi(iatom)
C*****Generate pair paramters (RijStar and eij) following combination
C*****rules
       RiiStar=Ai(iatom)*(alphai(iatom)**0.25d0)
       if(igas.eq.2)then
        Rsum=RiiStar+RN2Star
        gammaij=(RiiStar-RN2Star)/Rsum
        coeff1=-MMFF_beta*gammaij*gammaij
        RijStar(iatom)=0.5d0*Rsum*(1.0d0+MMFF_B*(1.0d0-dexp(coeff1)))

        coeff2=181.16d0*Gi(iatom)*GiN2*alphai(iatom)*alphaN2
        coeff3=dsqrt(alphai(iatom)/Ni(iatom))+dsqrt(alphaN2/NiN2)
       elseif(igas.eq.1)then
        Rsum=RiiStar+RHeStar
        gammaij=(RiiStar-RHeStar)/Rsum
        coeff1=-MMFF_beta*gammaij*gammaij
        RijStar(iatom)=0.5d0*Rsum*(1.0d0+MMFF_B*(1.0d0-dexp(coeff1)))

        coeff2=181.16d0*Gi(iatom)*GiHe*alphai(iatom)*alphaHe
        coeff3=dsqrt(alphai(iatom)/Ni(iatom))+dsqrt(alphaHe/NiHe)
       endif
       eij(iatom)=coeff2/(coeff3*(RijStar(iatom)**6.0d0))
C***Unit conversion from ang to meters and kcal/mol to J
       RijStar(iatom)=RijStar(iatom)*1.0d-10*distsca
       eij(iatom)=eij(iatom)*(4.184d3/xn)*enersca
C****well depth in MM3 is 1.1195eps
       eij(iatom)=eij(iatom)/1.1195d0
 
       tcharge=tcharge+pcharge(iatom)
       acharge=acharge+dabs(pcharge(iatom))
       if(unit.eq.'au') then
        fx(iatom)=fx(iatom)*0.52917706d0
        fy(iatom)=fy(iatom)*0.52917706d0
        fz(iatom)=fz(iatom)*0.52917706d0
       endif
2000  continue
      write(8,615) tcharge,acharge
c
      if(dchar.eq.'equal') then
       do 2011 iatom=1,inatom
c 2011 pcharge(iatom)=1.d0/dble(inatom)
        pcharge(iatom)=tcharge/dble(inatom)
2011   continue
      endif
c
      if(dchar.eq.'none') then
       do 2012 iatom=1,inatom
        pcharge(iatom)=0.d0
2012   continue
      endif
c
C      if(dchar.ne.'none') write(8,615) tcharge,acharge
C**********
      m2=0.d0
      do 2021 iatom=1,inatom
       m2=m2+xmass(iatom)
2021  continue
      write(8,604) m2
c
c
C      if(iu1.eq.1) write(8,620)
c
      fxo=0.d0
      fyo=0.d0
      fzo=0.d0
      do 2009 iatom=1,inatom
       fxo=fxo+(fx(iatom)*xmass(iatom))
       fyo=fyo+(fy(iatom)*xmass(iatom))
       fzo=fzo+(fz(iatom)*xmass(iatom))
2009  continue
      fxo=fxo/m2
      fyo=fyo/m2
      fzo=fzo/m2
C      write(8,623) fxo,fyo,fzo
      do 2010 iatom=1,inatom
      fx(iatom)=(fx(iatom)-fxo)*1.d-10
      fy(iatom)=(fy(iatom)-fyo)*1.d-10
      fz(iatom)=(fz(iatom)-fzo)*1.d-10
C      if(iu1.eq.1) write(8,600) fx(iatom),fy(iatom),fz(iatom),
C     ?imass(iatom),pcharge(iatom),eij(iatom)/xe,RijStar(iatom)*1.0d10
 2010 continue
C      if(iu1.eq.1) write(8,621)
c
      if(icoord.eq.1) close(9)
c
      do 3000 iatom=1,inatom
       ox(iatom)=fx(iatom)
       oy(iatom)=fy(iatom)
       oz(iatom)=fz(iatom)
3000  continue
c
      romax=0.d0
      do 3001 iatom=1,inatom
       if(RijStar(iatom).gt.romax) romax=RijStar(iatom)
3001  continue
c
      romax=romax+(1.1055d-10/2.0)
c
c     determine structural asymmetry parameter
c
      theta=0.d0
      asymp=0.d0
      do 5010 igamma=0,360,2
       do 5000 iphi=0,180,2
        gamma=dble(igamma)/cang
        phi=dble(iphi)/cang
        call rotate
        xyzsum=0.d0
        yzsum=0.d0
        do 5005 iatom=1,inatom
         xyz=dsqrt(fx(iatom)**2+fy(iatom)**2+fz(iatom)**2)
         yz=dsqrt(fy(iatom)**2+fz(iatom)**2)
         xyzsum=xyzsum+xyz
         yzsum=yzsum+yz
5005    continue
        hold=((pi/4.d0)*xyzsum)/yzsum
        if(hold.gt.asymp) asymp=hold
5000   continue
5010  continue
c
603   format(1x,'input file name = ',a50)
601   format(1x,'input file label = ',a50)
650   format(1x,'number of coordinate sets =',i5)
612   format(1x,'number of atoms =',i4)
611   format(1x,'coordinates in atomic units')
614   format(1x,'coordinates in angstroms')
610   format(1x,'units not specified')
616   format(1x,'No temperatures specified, exiting program')
617   format(1x,'Temperature grid not specified, 
     1  setting to default of 10')
613   format(1x,'correction factor for coordinates =',1pe11.4)
630   format(1x,'using a uniform charge distribution')
631   format(1x,'using a calculated (non-uniform) charge distribution')
633   format(1x,'using no charge - only LJ interactions')
632   format(1x,'charge distribution not specified')
615   format(1x,'total charge =',1pe11.4,/1x,
     1 'total absolute charge =',e11.4)
620   format(/9x,'initial coordinates',9x,'mass',3x,'charge',
     1 9x,'LJ parameters',/)
602   format(1x,'type not defined for atom number',i3)
604   format(1x,'mass of ion =',1pd11.4,' amu')
623   format(1x,'center of mass coordinates = ',1pe11.4,',',e11.4,
     1 ',',e11.4)
600   format(1x,1pe11.4,1x,e11.4,1x,e11.4,1x,i3,1x,
     1 e11.4,1x,e11.4,1x,e11.4)
621   format(/)
      return
      end
c
c     ***************************************************************
c
      subroutine rantate
c
c     Rotates the cluster/molecule to a random orientation.  
c
      use mpi
      implicit double precision (a-h,m-z)
      parameter (ixlen=1000,inpmax=201)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xeo,xk,xn,xmv,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(ixlen),k_const(ixlen),tcharge
      common/coordinates/fx(ixlen),fy(ixlen),fz(ixlen),
     ?ox(ixlen),oy(ixlen),oz(ixlen)
      common/ljparameters/eolj(ixlen),rolj(ixlen),eox4(ixlen),
     ?ro6lj(ixlen),ro12lj(ixlen),dro6(ixlen),dro12(ixlen)
      common/ljparameters2/alphai(ixlen),Ni(ixlen),Ai(ixlen),Gi(ixlen),
     1 RijStar(ixlen),eij(ixlen)
      common/mmffN2/alphaN2,NiN2,AiN2,GiN2,RN2Star,MMFF_B,MMFF_beta
      common/mmffHe/alphaHe,NiHe,AiHe,GiHe,RHeStar
      common/hsparameters/rhs(ixlen),rhs2(ixlen)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
      common/scaling/enersca,distsca
      common/mobilcal/b2max(inpmax),parab2max(inpmax)
      common/mpicon/imyrank,inprocs,imp_per_node,inp_per_node,s_time
c
      rnt=xrand()
      rnp=xrand()
      rng=xrand()
      theta=rnt*2.d0*pi
      phi=dasin((rnp*2.d0)-1.d0)+(pi/2.d0)
      gamma=rng*2.d0*pi
      call rotate
c
      return
      end
c
c     ***************************************************************
c
      subroutine rotate
c
c     Rotates the cluster/molecule.  
c
      use mpi
      implicit double precision (a-h,m-z)
      parameter (ixlen=1000,inpmax=201)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xeo,xk,xn,xmv,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(ixlen),k_const(ixlen),tcharge
      common/coordinates/fx(ixlen),fy(ixlen),fz(ixlen),
     ?ox(ixlen),oy(ixlen),oz(ixlen)
      common/ljparameters/eolj(ixlen),rolj(ixlen),eox4(ixlen),
     ?ro6lj(ixlen),ro12lj(ixlen),dro6(ixlen),dro12(ixlen)
      common/ljparameters2/alphai(ixlen),Ni(ixlen),Ai(ixlen),Gi(ixlen),
     1 RijStar(ixlen),eij(ixlen)
      common/mmffN2/alphaN2,NiN2,AiN2,GiN2,RN2Star,MMFF_B,MMFF_beta
      common/mmffHe/alphaHe,NiHe,AiHe,GiHe,RHeStar
      common/hsparameters/rhs(ixlen),rhs2(ixlen)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
      common/scaling/enersca,distsca
      common/mobilcal/b2max(inpmax),parab2max(inpmax)
      common/mpicon/imyrank,inprocs,imp_per_node,inp_per_node,s_time
      common/runparams/itn,inp,imp,igas
c
C      if(imyrank.eq.0)then
C       if(iu2.eq.1.or.iu3.eq.1) write(8,610) theta*cang,phi*cang,
C     1  gamma*cang
C      else
C       if(iu2.eq.1.or.iu3.eq.1) write(1000+imyrank,610)theta*cang,
C     1  phi*cang,
C     1  gamma*cang
C      endif
c
      do 1000 iatom=1,inatom
       rxy=dsqrt(ox(iatom)*ox(iatom)+(oy(iatom)*oy(iatom)))
       if(rxy.gt.0.0d0)then
        otheta=dacos(ox(iatom)/rxy)
        if(oy(iatom).lt.0.d0) otheta=(2.d0*pi)-otheta
        ntheta=otheta+theta
       endif
       fx(iatom)=dcos(ntheta)*rxy
       fy(iatom)=dsin(ntheta)*rxy
1000  continue
c
      do 2000 iatom=1,inatom
       rzy=dsqrt(oz(iatom)*oz(iatom)+(fy(iatom)*fy(iatom)))
       if(rzy.gt.0.d0)then
        ophi=dacos(oz(iatom)/rzy)
        if(fy(iatom).lt.0.d0) ophi=(2.d0*pi)-ophi
        nphi=ophi+phi
       endif
       fz(iatom)=dcos(nphi)*rzy
       fy(iatom)=dsin(nphi)*rzy
2000  continue
c
      do 3000 iatom=1,inatom
       rxy=dsqrt(fx(iatom)*fx(iatom)+(fy(iatom)*fy(iatom)))
       if(rxy.gt.0.d0)then
        ogamma=dacos(fx(iatom)/rxy)
        if(fy(iatom).lt.0.d0) ogamma=(2.d0*pi)-ogamma
        ngamma=ogamma+gamma
       endif
       fx(iatom)=dcos(ngamma)*rxy
       fy(iatom)=dsin(ngamma)*rxy
3000  continue
c
C      if(iu2.ne.0)then
C       if(imyrank.eq.0)then
C        write(8,620)
C       else
C        write(1000+imyrank,620)
C       endif
C       do iatom=1,inatom
C        if(imyrank.eq.0)then
C         write(8,600) ox(iatom),oy(iatom),oz(iatom),fx(iatom),
C     1    fy(iatom),fz(iatom)
C        else
C         write(1000+imyrank,600) ox(iatom),oy(iatom),
C     1    oz(iatom),fx(iatom),
C     1    fy(iatom),fz(iatom)
C        endif
C       enddo
C      endif
c
600   format(1x,1pe11.4,2(1x,e11.4),5x,3(1x,e11.4))
610   format(//1x,'coordinates rotated by ROTATE',//1x,
     1 'theta=',1pe11.4,1x,'phi=',e11.4,1x,'gamma=',1pe11.4,/)
620   format(9x,'initial coordinates',24x,'new coordinates',/)
      return
      end
c
c     *****************************************************
c
      subroutine dljpotHe(x,y,z,pot,dpotx,dpoty,dpotz,dmax,dchar)
c
c     Subroutine to calculate L-J + ion-dipole potential.
c
      use mpi
      implicit double precision (a-h,m-z)
      parameter (ixlen=1000,inpmax=201)
      character*4 dchar
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xeo,xk,xn,xmv,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(ixlen),k_const(ixlen),tcharge
      common/coordinates/fx(ixlen),fy(ixlen),fz(ixlen),
     ?ox(ixlen),oy(ixlen),oz(ixlen)
      common/ljparameters/eolj(ixlen),rolj(ixlen),eox4(ixlen),
     ?ro6lj(ixlen),ro12lj(ixlen),dro6(ixlen),dro12(ixlen)
      common/ljparameters2/alphai(ixlen),Ni(ixlen),Ai(ixlen),Gi(ixlen),
     1 RijStar(ixlen),eij(ixlen)
      common/mmffN2/alphaN2,NiN2,AiN2,GiN2,RN2Star,MMFF_B,MMFF_beta
      common/mmffHe/alphaHe,NiHe,AiHe,GiHe,RHeStar
      common/hsparameters/rhs(ixlen),rhs2(ixlen)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
      common/scaling/enersca,distsca
      common/mobilcal/b2max(inpmax),parab2max(inpmax)
      common/mpicon/imyrank,inprocs,imp_per_node,inp_per_node,s_time
      common/runparams/itn,inp,imp,igas
      dimension pottry(2,3)
      dimension pot_mol(3)
      dimension dpotxtry(2,3)
      dimension dpotx_mol(3)
      dimension dpotytry(2,3)
      dimension dpoty_mol(3)
      dimension dpotztry(2,3)
      dimension dpotz_mol(3)
c
c     nitrogen : five charge model (Allen Tildesley, 14 page)
c     Every data from B3LYP//aug-cc-pVDZ
      dmax=2.d0*romax
      rx=0.d0
      ry=0.d0
      rz=0.d0
      e00=0.d0
      de00x=0.d0
      de00y=0.d0
      de00z=0.d0
      sum1=0.d0
      sum2=0.d0
      sum3=0.d0
      sum4=0.d0
      sum5=0.d0
      sum6=0.d0
      do 1100 iatom=1,inatom
c
       xx=x-fx(iatom)
       xx2=xx*xx
c
       yy=y-fy(iatom)
       yy2=yy*yy
c
       zz=z-fz(iatom)
       zz2=zz*zz
c
       rxyz2=xx2+yy2+zz2
c
       rxyz=dsqrt(rxyz2)
C*****
C***** Exp6 potential calculation
       R=rxyz/RijStar(iatom)
       R2=R*R
       R6=R2*R2*R2
       preexp=1.84d5
       expterm=preexp*dexp(-12.0d0*R)
       predisp=2.25d0
       dispterm=predisp/R6
C*****

       if(rxyz.lt.dmax) dmax=rxyz
       rxyz3=rxyz2*rxyz
       rxyz5=rxyz3*rxyz2
C**** Exp6 potential
       e00=e00+eij(iatom)*(expterm-dispterm)
C**** Exp6 derivative
       de00=eij(iatom)*(expterm*(-12.0d0)+dispterm*6.0d0/R)
     1   /rxyz/RijStar(iatom)
       de00x=de00x+(de00*xx)
       de00y=de00y+(de00*yy)
       de00z=de00z+(de00*zz)
C     ion-induced dipole potential
       if(dchar.ne.'none')then
        if(pcharge(iatom).ne.0.d0)then
         rxyz3i=pcharge(iatom)/rxyz3
         rxyz5i=-3.d0*pcharge(iatom)/rxyz5
         rx=rx+(xx*rxyz3i)
         ry=ry+(yy*rxyz3i)
         rz=rz+(zz*rxyz3i)
C     ion-induced dipole derivative
         sum1=sum1+(rxyz3i+(xx2*rxyz5i))
         sum2=sum2+(xx*yy*rxyz5i)
         sum3=sum3+(xx*zz*rxyz5i)
         sum4=sum4+(rxyz3i+(yy2*rxyz5i))
         sum5=sum5+(yy*zz*rxyz5i)
         sum6=sum6+(rxyz3i+(zz2*rxyz5i))
        endif
       endif
c
 1100 continue
c
      pot=e00-(dipol*((rx*rx)+(ry*ry)+(rz*rz)))
      dpotx=de00x-(dipol*((2.0d0*rx*sum1)+(2.0d0*ry*sum2)
     1   +(2.0d0*rz*sum3)))
      dpoty=de00y-(dipol*((2.0d0*rx*sum2)+(2.0d0*ry*sum4)
     1   +(2.0d0*rz*sum5)))
      dpotz=de00z-(dipol*((2.0d0*rx*sum3)+(2.0d0*ry*sum5)
     1   +(2.0d0*rz*sum6)))
c
      return
      end
c
c     *******************************************************
c     *****************************************************
c
      subroutine dljpotN2(x,y,z,pot,dpotx,dpoty,dpotz,dmax,dchar)
c
c     Subroutine to calculate L-J + ion-dipole potential.
c
      use mpi
      implicit double precision (a-h,m-z)
      parameter (ixlen=1000,inpmax=201)
      character*4 dchar
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xeo,xk,xn,xmv,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(ixlen),k_const(ixlen),tcharge
      common/coordinates/fx(ixlen),fy(ixlen),fz(ixlen),
     ?ox(ixlen),oy(ixlen),oz(ixlen)
      common/ljparameters/eolj(ixlen),rolj(ixlen),eox4(ixlen),
     ?ro6lj(ixlen),ro12lj(ixlen),dro6(ixlen),dro12(ixlen)
      common/ljparameters2/alphai(ixlen),Ni(ixlen),Ai(ixlen),Gi(ixlen),
     1 RijStar(ixlen),eij(ixlen)
      common/mmffN2/alphaN2,NiN2,AiN2,GiN2,RN2Star,MMFF_B,MMFF_beta
      common/mmffHe/alphaHe,NiHe,AiHe,GiHe,RHeStar
      common/hsparameters/rhs(ixlen),rhs2(ixlen)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/temperature/tt1,tt2,ittgrid
      common/xrandom/i1,i2,i3,i4,i5,i6
      common/scaling/enersca,distsca
      common/mobilcal/b2max(inpmax),parab2max(inpmax)
      common/mpicon/imyrank,inprocs,imp_per_node,inp_per_node,s_time
      common/runparams/itn,inp,imp,igas
      dimension pottry(2,3)
      dimension pot_mol(3)
      dimension dpotxtry(2,3)
      dimension dpotx_mol(3)
      dimension dpotytry(2,3)
      dimension dpoty_mol(3)
      dimension dpotztry(2,3)
      dimension dpotz_mol(3)
c
c     nitrogen : five charge model (Allen Tildesley, 14 page)
c     Every data from B3LYP//aug-cc-pVDZ
      dmax=2.d0*romax
      bond=1.0976d-10
      Ptfn=0.d0
      xkT=tt1*xk
      pc=-0.4825d0
c     -1.447 D A 
      pc_center=-(pc)
c     B3LYP/aug-cc-pVDZ
c      dipolzz=2.263d-30/(2.d0*4.d0*pi*xeo)*1.0
      dipolzz=1.710d-30/(2.d0*4.d0*pi*xeo)
c      dipolzz=1.7543d-30/(2.d0*4.d0*pi*xeo)
      dipolzz=dipolzz*xe*xe 
      dipolxx=1.710d-30/(2.d0*4.d0*pi*xeo)
c      dipolxx=1.500d-30/(2.d0*4.d0*pi*xeo)*1.0
c      dipolxx=1.7543d-30/(2.d0*4.d0*pi*xeo)
      dipolxx=dipolxx*xe*xe 
      pot_min=1.0d8
c
      do 1300 isamp=1,3
       do 1200 ibatom=1,2
        rx=0.d0
        ry=0.d0
        rz=0.d0
        e00=0.d0
        de00x=0.d0
        de00y=0.d0
        de00z=0.d0
        sum1=0.d0
        sum2=0.d0
        sum3=0.d0
        sum4=0.d0
        sum5=0.d0
        sum6=0.d0
        qpol=0.d0
        dqpolx=0.d0
        dqpoly=0.d0
        dqpolz=0.d0
        xc=0.0d0
        yc=0.0d0
        zc=0.0d0
        if(isamp.eq.1)then
         xc=(bond/2.d0)*(2.d0*dble(ibatom)-3.d0)
         dpolx=dipolzz
         dpoly=dipolxx
         dpolz=dipolxx
        elseif(isamp.eq.2)then
         yc=(bond/2.d0)*(2.d0*dble(ibatom)-3.d0)
         dpolx=dipolxx
         dpoly=dipolzz
         dpolz=dipolxx
        elseif(isamp.eq.3)then
         zc=(bond/2.d0)*(2.d0*dble(ibatom)-3.d0)
         dpolx=dipolxx
         dpoly=dipolxx
         dpolz=dipolzz
        endif
c
        do 1100 iatom=1,inatom
c
         xx_center=x-fx(iatom)
         xx=xx_center+xc
         xx_center2=xx_center*xx_center
         xx2=xx*xx
c
         yy_center=y-fy(iatom)
         yy=yy_center+yc
         yy_center2=yy_center*yy_center
         yy2=yy*yy
c
         zz_center=z-fz(iatom)
         zz=zz_center+zc
         zz_center2=zz_center*zz_center
         zz2=zz*zz
c
         rxyz_center2=xx_center2+yy_center2+zz_center2
         rxyz2=xx2+yy2+zz2
c
         rxyz_center=dsqrt(rxyz_center2)
         rxyz=dsqrt(rxyz2)
C*****
C***** vdW potential -> CoM
         R=rxyz/RijStar(iatom)
         R2=R*R
         R6=R2*R2*R2
         R7=R6*R
C***** Exp6 potential calculation
         preexp=1.84d5
         expterm=preexp*dexp(-12.0d0*R)
         predisp=2.25d0
         dispterm=predisp/R6
         e00pre=(expterm-dispterm)
C*****
c
         if(rxyz.lt.dmax) dmax=rxyz
         rxyz3=rxyz2*rxyz
         rxyz5=rxyz3*rxyz2
         rxyz6=rxyz5*rxyz
         rxyz8=rxyz5*rxyz3
         rxyz12=rxyz6*rxyz6
         rxyz14=rxyz12*rxyz2
         rxyz_center3=rxyz_center2*rxyz_center
         rxyz_center5=rxyz_center3*rxyz_center2
C**** vdW potential
         e00=e00+eij(iatom)*e00pre
C**** vdW derivative
c****   Exp-6 derivative
         de00=eij(iatom)*(expterm*(-12.0d0)+dispterm*6.0d0/R)
     1      /rxyz/RijStar(iatom)
c     cartesian components
         de00x=de00x+(de00*xx)
         de00y=de00y+(de00*yy)
         de00z=de00z+(de00*zz)
C     ion-induced dipole potential
         if(dchar.ne.'none')then
          if(pcharge(iatom).ne.0.d0)then
           rxyz3i=pcharge(iatom)/rxyz_center3
           rxyz5i=-3.d0*pcharge(iatom)/rxyz_center5
           rx=rx+(xx_center*rxyz3i)
           ry=ry+(yy_center*rxyz3i)
           rz=rz+(zz_center*rxyz3i)
C     ion-induced dipole derivative
           sum1=sum1+(rxyz3i+(xx_center2*rxyz5i))
           sum2=sum2+(xx_center*yy_center*rxyz5i)
           sum3=sum3+(xx_center*zz_center*rxyz5i)
           sum4=sum4+(rxyz3i+(yy_center2*rxyz5i))
           sum5=sum5+(yy_center*zz_center*rxyz5i)
           sum6=sum6+(rxyz3i+(zz_center2*rxyz5i))
c     ion-partial charge coulomb potential(quadrupole)       
           const_k=pcharge(iatom)*(xe*xe)/(4.d0*pi*xeo)
           qpol=qpol+(pc_center*const_k/rxyz_center)
           qpol=qpol+(pc*const_k/rxyz)
c     ion-partial charge coulomb derivative(quadrupole)       
           dqpolx=dqpolx-((pc_center*const_k/rxyz_center3)*(xx_center))
           dqpoly=dqpoly-((pc_center*const_k/rxyz_center3)*(yy_center))
           dqpolz=dqpolz-((pc_center*const_k/rxyz_center3)*(zz_center))
           dqpolx=dqpolx-((pc*const_k/rxyz3)*(xx))
           dqpoly=dqpoly-((pc*const_k/rxyz3)*(yy))
           dqpolz=dqpolz-((pc*const_k/rxyz3)*(zz))
          endif
         endif
c
 1100   continue
c
c*AH    half the vdW pot and its derivative for averaging both Ns
        e00=e00*0.5d0
        de00x=de00x*0.5d0
        de00y=de00y*0.5d0
        de00z=de00z*0.5d0
c*AH
        pottry(ibatom,isamp)=e00-0.5d0*(((dpolx*rx*rx)+(dpoly*ry*ry)
     1   +(dpolz*rz*rz)))+qpol
        dpotxtry(ibatom,isamp)=de00x-0.5d0*((dpolx*2.0d0*rx*sum1)
     1   +(dpoly*2.0d0*ry*sum2)+(dpolz*2.d0*rz*sum3))+dqpolx
        dpotytry(ibatom,isamp)=de00y-0.5d0*((dpolx*2.0d0*rx*sum2)
     1   +(dpoly*2.0d0*ry*sum4)+(dpolz*2.d0*rz*sum5))+dqpoly
        dpotztry(ibatom,isamp)=de00z-0.5d0*((dpolx*2.0d0*rx*sum3)
     1   +(dpoly*2.0d0*ry*sum5)+(dpolz*2.0d0*rz*sum6))+dqpolz
c 
c
 1200  continue
       pot_mol(isamp)=pottry(1,isamp)+pottry(2,isamp)
       tpot=pot_mol(isamp)
       if(pot_min.ge.pot_mol(isamp)) pot_min=pot_mol(isamp)
       dpotx_mol(isamp)=dpotxtry(1,isamp)+dpotxtry(2,isamp)
       dpoty_mol(isamp)=dpotytry(1,isamp)+dpotytry(2,isamp)
       dpotz_mol(isamp)=dpotztry(1,isamp)+dpotztry(2,isamp)
 1300 continue
c
      do 1410 isamp=1,3
       temp_pot=pot_mol(isamp)-pot_min
       Ptfn=Ptfn+dexp(-temp_pot/xkT)
 1410 continue
c
      pot=0.d0
      dpotx=0.d0
      dpoty=0.d0
      dpotz=0.d0
      do 1400 isamp=1,3
       temp_pot=pot_mol(isamp)-pot_min
       weight=dexp(-temp_pot/xkT)/Ptfn
       pot=pot+weight*pot_mol(isamp)
       dpotx=dpotx+weight*dpotx_mol(isamp)
       dpoty=dpoty+weight*dpoty_mol(isamp)
       dpotz=dpotz+weight*dpotz_mol(isamp)
 1400 continue
c
      return
      end
c
c     *******************************************************
c
      subroutine gsang(v,b,erat,ang,d1,istep)
c
c     Calculates trajectory. Adapted from code written by George Schatz
c     A fixed step integrator is used which combines a runge-kutta-gill
c     initiator with an adams-moulton predictor-corrector propagator.
c     
      use mpi
      implicit double precision (a-h,m-z)
      integer ns,nw,l
      dimension w(6),dw(6)
      parameter (ixlen=1000,inpmax=201)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xeo,xk,xn,xmv,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(ixlen),k_const(ixlen),tcharge
      common/coordinates/fx(ixlen),fy(ixlen),fz(ixlen),
     ?ox(ixlen),oy(ixlen),oz(ixlen)
      common/ljparameters/eolj(ixlen),rolj(ixlen),eox4(ixlen),
     ?ro6lj(ixlen),ro12lj(ixlen),dro6(ixlen),dro12(ixlen)
      common/ljparameters2/alphai(ixlen),Ni(ixlen),Ai(ixlen),Gi(ixlen),
     1 RijStar(ixlen),eij(ixlen)
      common/mmffN2/alphaN2,NiN2,AiN2,GiN2,RN2Star,MMFF_B,MMFF_beta
      common/mmffHe/alphaHe,NiHe,AiHe,GiHe,RHeStar
      common/hsparameters/rhs(ixlen),rhs2(ixlen)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
      common/scaling/enersca,distsca
      common/mobilcal/b2max(inpmax),parab2max(inpmax)
      common/mpicon/imyrank,inprocs,imp_per_node,inp_per_node,s_time
      common/runparams/itn,inp,imp,igas
c
      vy=-v
      vx=0.d0
      vz=0.d0
      vxyz=dabs(vy)
c
c     determine time step
c
      top=(v/95.2381d0)-0.5d0
      if(v.ge.1000.d0) top=10.d0
      if(v.ge.2000.d0) top=10.0d0-((v-2000.d0)*7.5d-3)
      if(v.ge.3000.d0) top=2.5d0
c
      dt1=top*dtsf1*1.0d-11/v
      dt2=dt1*dtsf2
      dt=dt1
c
c     determine trajectory start position
c
      e0=0.5d0*mu*v*v
      x=b
      z=0.d0     
c
      ymin=0.d0
      ymax=0.d0
C      do 200 i=1,inatom
C       if(fy(i).gt.ymax) ymax=fy(i)
C       if(fy(i).lt.ymin) ymin=fy(i)
C200   continue
      ymax=ymax/1.0d-10
      ymin=ymin/1.0d-10
      iymin=nint(ymin)-1
      iymax=nint(ymax)+1
C***********************************
C***********************************
      id2=iymax
      istop=0
      iswitch=0
      ireturn=0
      icheck=0
      pot=0.0d0
      y=dble(id2)*1.0d-10
      if(igas.eq.2)then
       call dljpotN2(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
      elseif(igas.eq.1)then
       call dljpotHe(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
      endif
      if(dabs(pot/e0).gt.sw1)then
C****
       do while(dabs(pot/e0).gt.sw1)
        id2=id2+10
        y=dble(id2)*1.0d-10
        if(igas.eq.2)then
         call dljpotN2(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
        elseif(igas.eq.1)then
         call dljpotHe(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
        endif
       enddo
       do while(dabs(pot/e0).lt.sw1)
        id2=id2-1
        y=dble(id2)*1.0d-10
        if(igas.eq.2)then
         call dljpotN2(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
        elseif(igas.eq.1)then
         call dljpotHe(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
        endif
       enddo
      else
       do while(dabs(pot/e0).lt.sw1)
        id2=id2-1
        y=dble(id2)*1.0d-10
        if(igas.eq.2)then
         call dljpotN2(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
        elseif(igas.eq.1)then
         call dljpotHe(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
        endif
        if((id2.lt.iymin).and.(icheck.eq.0)) then
C        if(it.eq.1) write(8,621)
         ang=0.d0
         erat=1.000d0
         ireturn=1
        endif
       enddo
C****
      endif
      if(ireturn.eq.1)return
      y=dble(id2)*1.0d-10
      etot=e0+pot
C      if(it.eq.1) write(8,622) y*1.0d10
      d1=y
C***********************************
c     initial coordinates and momenta
c
      w(1)=x
      w(2)=vx*mu
      w(3)=y
      w(4)=vy*mu
      w(5)=z
      w(6)=vz*mu
      tim=0.d0
    
C      if(it.eq.1) write(8,623)  
c
c     initialize the time derivatives of the coordinates and momenta
c
      call deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
      ns=0
      nw=0
      l=0
C**********************************
      istop=0
      ireturn=0
      do while(istop.eq.0)
       call diffeq(l,tim,dt,w,dw,pot,dmax)
       nw=nw+1
       if (nw.eq.inwr)then
        ns=ns+nw
        nw=0
c     check if trajectory has become "lost" (too many steps)
c     b is in m, v is in m/s
c
        if(ns.gt.30000) then
         if(imyrank.eq.0)then
          write(8,105) b,v
         else 
          write(1000+imyrank,105) b,v
         endif
         ang=pi/2.d0
         e=0.5d0*mu*(dw(1)**2+dw(3)**2+dw(5)**2)
         erat=(e+pot)/etot
         istep=ns
         ireturn=1
         istop=1
        else
c
c     check if the trajectory is finished
c
         istop=1
         if(dmax.lt.romax)istop=0
         if(dabs(pot/e0).gt.sw2.and.dt.eq.dt1.and.istop.eq.1) then
          dt=dt2
          l=0
         endif
         if(dabs(pot/e0).lt.sw2.and.dt.eq.dt2.and.istop.eq.1) then
          dt=dt1 
          l=0
         endif
         if(dabs(pot/e0).gt.sw1)istop=0
         if(ns.lt.50)istop=0
        endif
c
       endif
      enddo
      if(ireturn.eq.1)return
C**********************************
      istep=ns
c
c     determine scattering angle 
c
      if(dw(1).gt.0.d0) then
       num=dw(3)*(-v)
       den=v*dsqrt(dw(1)**2+dw(3)**2+dw(5)**2)
       ang=dacos(num/den)
      endif
      if(dw(1).lt.0.d0) then
       num=dw(3)*(-v)
       den=v*dsqrt(dw(1)**2+dw(3)**2+dw(5)**2)
       ang=(-dacos(num/den))
      endif
c
c     check for energy conservation
c
      e=0.5d0*mu*(dw(1)**2+dw(3)**2+dw(5)**2)
      erat=(e+pot)/etot
622   format(1x,'trajectory start position =',1pe11.4)
621   format(1x,'trajectory not started - potential too small')
623   format(//1x,'trajectory ns, x,  y,  z,  kin e, dt,    tot e',/1x,
     1 '               vx, vy, vz, pot e, pot/e0'/)    
105   format(1x,'trajectory lost: b =',1pe11.4,' v =',e11.4)
108   format(/1x,'specific trajectory parameters',//1x,
     1 'v =',1pe11.4,4x,'b =',e11.4)
107   format(1x,'time steps, dt1 =',1pe11.4,' dt2 =',e11.4)
      return
c       
c
      end
c
c     *******************************************************
c
      subroutine diffeq(l,tim,dt,w,dw,pot,dmax)
c
c     Integration subroutine - uses 5th order runge-kutta-gill to 
c     initiate and 5th order adams-moulton predictor-corrector to 
c     propagate. Parameter l is initially set to zero and then 
c     incremented to tell the subroutine when to switch between 
c     integration methods. DIFFEQ calls subroutine DERIV to define 
c     the equations of motion to be integrated.
c
      use mpi
      implicit double precision (a-h,m-z)
      dimension w(6),dw(6)
      parameter (ixlen=1000,inpmax=201) 
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xeo,xk,xn,xmv,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(ixlen),k_const(ixlen),tcharge
      common/coordinates/fx(ixlen),fy(ixlen),fz(ixlen),
     ?ox(ixlen),oy(ixlen),oz(ixlen)
      common/ljparameters/eolj(ixlen),rolj(ixlen),eox4(ixlen), 
     ?ro6lj(ixlen),ro12lj(ixlen),dro6(ixlen),dro12(ixlen)
      common/ljparameters2/alphai(ixlen),Ni(ixlen),Ai(ixlen),Gi(ixlen),
     1 RijStar(ixlen),eij(ixlen)
      common/mmffN2/alphaN2,NiN2,AiN2,GiN2,RN2Star,MMFF_B,MMFF_beta
      common/mmffHe/alphaHe,NiHe,AiHe,GiHe,RHeStar
      common/hsparameters/rhs(ixlen),rhs2(ixlen)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
      common/scaling/enersca,distsca
      common/mobilcal/b2max(inpmax),parab2max(inpmax)
      common/mpicon/imyrank,inprocs,imp_per_node,inp_per_node,s_time
      common/runparams/itn,inp,imp,igas
      dimension a(4),b(4),c(4),ampc(5),amcc(4),array(6,40),savw(40),
     1 savdw(40),q(40)
      data a/0.50d0,0.292893218814d0,1.70710678118d0,0.1666666666667d0/
      data b/2.0d0,1.0d0,1.0d0,2.0d0/
      data c/-0.5d0,-0.292893218814d0,-1.70710678118d0,-0.5d0/
      data ampc/-0.111059153612d0,0.672667757774d0,-1.70633621697d0,
     ?2.33387888707d0,-1.8524668225d0/
      data amcc/0.0189208128941d0,-0.121233356692d0,0.337771548703d0,
     ?-0.55921513665d0/
      data var,cvar,acst/2.97013888888d0,0.990972222222d0,
     ?0.332866152768d0/
      save hvar,hcvar
c
      if(l.ge.0)then
C      if (l) 4,1,3
       if(l.eq.0)then
        do 2 j=1,6
         q(j)=0.0
2       continue
        hvar=dt*var
        hcvar=dt*cvar
        dt=0.5*dt
       endif
       l=l+1
c
c     This is the runge-kutta-gill part...the steps are broken up into
c     half stekps to improve accuracy.
c
       k=0
       istop=0
       do while(istop.eq.0)
        do 17 j=1,4
         if (((-1)**j).gt.0) tim=tim+0.5*dt
         call deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
         do 7 i=1,6
          dw(i)=dt*dw(i)
          r=a(j)*(dw(i)-b(j)*q(i))
          w(i)=w(i)+r
          q(i)=q(i)+3.0*r+c(j)*dw(i)
7        continue
17      continue
        call deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
        if(k.gt.0)istop=1
        if(k.le.0)k=1
       enddo
C********
       if((l-6).ge.0)then
        l=-1
        dt=2.0*dt
       else
        do 8 j=1,6
         array(l,j)=dw(j)
8       continue
       endif
C********
      else
c     This is the adams-moulton predictor-corrector part.
c
       do 10 j=1,6
        savw(j)=w(j)
        savdw(j)=dw(j)
        array(6,j)=savdw(j)
        do 9 i=1,5
         array(6,j)=array(6,j)+ampc(i)*array(i,j)
9       continue
        w(j)=array(6,j)*hvar+w(j)
10     continue
       tim=tim+dt
       call deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
       do 12 j=1,6
        array(6,j)=acst*dw(j)
        do 11 i=1,4
         array(i,j)=array(i+1,j)
         array(6,j)=array(i,j)*amcc(i)+array(6,j)
11      continue
        array(5,j)=savdw(j)
        w(j)=savw(j)+hcvar*(array(5,j)+array(6,j))
12     continue
       call deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
      endif
C********************
      return
      end
c
c     ***************************************************************
c
      subroutine deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
c
c     Defines Hamilton's equations of motion as the time derivatives 
c     of the coordinates and momenta.
c
      use mpi
      implicit double precision (a-h,m-z)
      dimension w(6),dw(6)
      parameter (ixlen=1000,inpmax=201)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xeo,xk,xn,xmv,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(ixlen),k_const(ixlen),tcharge
      common/coordinates/fx(ixlen),fy(ixlen),fz(ixlen),
     ?ox(ixlen),oy(ixlen),oz(ixlen)
      common/ljparameters/eolj(ixlen),rolj(ixlen),eox4(ixlen),
     ?ro6lj(ixlen),ro12lj(ixlen),dro6(ixlen),dro12(ixlen)
      common/ljparameters2/alphai(ixlen),Ni(ixlen),Ai(ixlen),Gi(ixlen),
     1 RijStar(ixlen),eij(ixlen)
      common/mmffN2/alphaN2,NiN2,AiN2,GiN2,RN2Star,MMFF_B,MMFF_beta
      common/mmffHe/alphaHe,NiHe,AiHe,GiHe,RHeStar
      common/hsparameters/rhs(ixlen),rhs2(ixlen)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
      common/scaling/enersca,distsca
      common/mobilcal/b2max(inpmax),parab2max(inpmax)
      common/mpicon/imyrank,inprocs,imp_per_node,inp_per_node,s_time
      common/runparams/itn,inp,imp,igas
c
c     From Hamilton's equations, the time derivatives of the coordinates
c     are the conjugates divided by the mass.
c
      dw(1)=w(2)/mu
      dw(3)=w(4)/mu
      dw(5)=w(6)/mu
c
c     Hamilton's equations for the time derivatives of the momenta
c     evaluated by using the coordinate derivatives together with the
c     chain rule.
c
      x=w(1)
      y=w(3)
      z=w(5)
c
c    These are analytical derivatives.
C
c
      if(igas.eq.2)then
       call dljpotN2(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
      elseif(igas.eq.1)then
       call dljpotHe(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
      endif
      dw(2)=-(dpotx)
      dw(4)=-(dpoty)
      dw(6)=-(dpotz)
c
c
      return
      end
c
c     ***************************************************************
c
      subroutine pre_mobil2 ()
c
c     Subroutine to determine grids for mobil2 
c
      use mpi
      implicit double precision (a-h,m-z)
      parameter (ixlen=1000,inpmax=201)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xeo,xk,xn,xmv,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(ixlen),k_const(ixlen),tcharge
      common/coordinates/fx(ixlen),fy(ixlen),fz(ixlen),
     ?ox(ixlen),oy(ixlen),oz(ixlen)
      common/ljparameters/eolj(ixlen),rolj(ixlen),eox4(ixlen),
     ?ro6lj(ixlen),ro12lj(ixlen),dro6(ixlen),dro12(ixlen)
      common/ljparameters2/alphai(ixlen),Ni(ixlen),Ai(ixlen),Gi(ixlen),
     1 RijStar(ixlen),eij(ixlen)
      common/mmffN2/alphaN2,NiN2,AiN2,GiN2,RN2Star,MMFF_B,MMFF_beta
      common/mmffHe/alphaHe,NiHe,AiHe,GiHe,RHeStar
      common/hsparameters/rhs(ixlen),rhs2(ixlen)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/temperature/tt1,tt2,ittgrid
      common/xrandom/i1,i2,i3,i4,i5,i6
      common/scaling/enersca,distsca
C****
      common/mobilcal/b2max(inpmax),parab2max(inpmax)
      common/mobilgrid/pgst(inpmax),wgst(inpmax),q1st(inpmax),
     1  q2st(inpmax),q3st(inpmax),q1STD(inpmax),q2STD(inpmax),
     1  q3STD(inpmax)
      common/mpicon/imyrank,inprocs,imp_per_node,inp_per_node,s_time
      common/runparams/itn,inp,imp,igas
      dimension q12D(inpmax,50),q22D(inpmax,50),q32D(inpmax,50)
      dimension cosx(0:500)
c     
      if(im2.eq.0)then
       if(imyrank.eq.0)then
        write(8,631)
        write(8,603) sw1,sw2,dtsf1,dtsf2,inwr,ifail
       else
        write(1000+imyrank,631)
        write(1000+imyrank,603) sw1,sw2,dtsf1,dtsf2,inwr,ifail
       endif
      endif
      it=0
      iu2=0
c
c     determine maximum extent and orientate along x axis
c
      if(im2.eq.0)then
       if(imyrank.eq.0)then
        write(8,632)
       else
        write(1000+imyrank,632)
       endif
      endif
      rmax=0.d0
      do 1000 iatom=1,inatom
       r=dsqrt((ox(iatom)*ox(iatom))+(oy(iatom)*oy(iatom))+
     1  (oz(iatom)*oz(iatom)))
       if(r.gt.rmax) then
        rmax=r
        ihold=iatom
       endif
1000  continue
c
      rzy=dsqrt((oz(ihold)*oz(ihold))+(oy(ihold)*oy(ihold)))
      phi=dacos(oz(ihold)/rzy)
      phi=phi+(pi/2.d0)
      if(oy(ihold).lt.0.d0) phi=(2.d0*pi)-phi
      phi=(2.d0*pi)-phi
      theta=0.d0 
      gamma=0.d0
      call rotate
      rxy=dsqrt((fx(ihold)*fx(ihold))+(fy(ihold)*fy(ihold)))
      gamma=dacos(fx(ihold)/rxy)
      if(fy(ihold).lt.0.d0) gamma=(2.d0*pi)-gamma
      gamma=(2.d0*pi)-gamma
      if(im2.eq.0) iu3=1
      if(ip.eq.1) iu2=1
      call rotate
      iu3=0
      iu2=0
      hold=fx(ihold)/rmax
      if(hold.lt.0.9999999999d0.or.hold.gt.1.0000000001d0.or.
     1 fy(ihold).gt.1.0d-20.or.fz(ihold).gt.1.0d-20.or.
     1 fy(ihold).lt.-1.0d-20.or.fz(ihold).lt.-1.0d-20) then 
       if(imyrank.eq.0)write(8,601)
       do 1001 iatom=1,inatom
        hold=dsqrt(fx(iatom)**2+fy(iatom)**2+fz(iatom)**2)
        if(imyrank.eq.0)then
         write(8,602) iatom,fx(iatom),fy(iatom),fz(iatom),hold
        else
         write(1000+imyrank,602) iatom,fx(iatom),fy(iatom),fz(iatom),
     1    hold
        endif
1001   continue     
       if(imyrank.eq.0)close(8)
       call MPI_FINALIZE(ierr)
       stop
      endif
c
c     determine rmax, emax, and r00 along x, y, and z directions
c
C     if(ip.eq.1) write(8,689)
      irn=1000
      ddd=(rmax+romax)/dble(irn)
c
      y=0.d0
      z=0.d0
      emaxx=0.d0
      istop=0
       do 1101 ir=1,irn
        x=rmax+romax-(dble(ir)*ddd) 
        if(igas.eq.2)then
         call dljpotN2(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
        elseif(igas.eq.1)then
         call dljpotHe(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
        endif
C*** Changed this back to break loop. Didn't make sense before
        if(pot.gt.0.d0)then
         goto 1170
        else
         r00x=x
         if(pot.lt.emaxx) then
          rmaxx=x
          emaxx=pot
         endif
        endif
1101  continue
1170  continue
      if(im2.eq.0)then
       if(imyrank.eq.0)then
        write(8,614) emaxx/xe,rmaxx*1.0d10,r00x*1.0d10
       else
        write(1000+imyrank,614) emaxx/xe,rmaxx*1.0d10,r00x*1.0d10
       endif
      endif    
c
C      x=0.d0
C      z=0.d0
C      emaxy=0.d0
C      istop=0
C      do while(istop.eq.0)
C       do 1100 ir=1,irn
C        y=rmax+romax-(dble(ir)*ddd)
C        call dljpot(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
C        if(pot.gt.0.d0)then
C         istop=1
C        else
C         r00y=y
C         if(pot.lt.emaxy) then
C          rmaxy=y
C          emaxy=pot
C         endif
C        endif
C1100   continue
C      enddo
C      if(im2.eq.0)then
C       if(imyrank.eq.0)then
C        write(8,613) emaxy/xe,rmaxy*1.0d10,r00y*1.0d10
C       else
C        write(1000+imyrank,613) emaxy/xe,rmaxy*1.0d10,r00y*1.0d10
C       endif
C      endif
c
C      x=0.d0
C      y=0.d0
C      emaxz=0.d0
C      istop=0
C      do while(istop.eq.0)
C       do 1102 ir=1,irn
C        z=rmax+romax-(dble(ir)*ddd)
C        call dljpot(x,y,z,pot,dpotx,dpoty,dpotz,dmax,'calc')
C        if(pot.gt.0.d0)then
C         istop=1
C        else
C         r00z=z
C         if(pot.lt.emaxz) then
C          rmaxz=z
C          emaxz=pot
C         endif
C        endif
C1102   continue
C      enddo
C      if(im2.eq.0)then
C       if(imyrank.eq.0)then
C        write(8,615) emaxz/xe,rmaxz*1.0d10,r00z*1.0d10
C       else
C        write(1000+imyrank,615) emaxz/xe,rmaxz*1.0d10,r00z*1.0d10
C       endif
C      endif
c
c     set-up integration over gst
c     
      tst=xk*tt1/eo
C****PATCH PATCH: vestigial printing. remove?
C      if(im2.eq.0)then 
C       if(imyrank.eq.0)then
C        write(8,600) tst
C       else
C        write(1000+imyrank,600) tst
C       endif
C      endif
      tst3=tst*tst*tst
c
C*      dgst=5.0d-7*6.d0*dsqrt(tst)
C*      gst=dgst
C*      sum=0.d0
C*      sum1=0.d0
C*      sum2=0.d0
c
C*      do 2020 i=1,inp
C*       sum1=sum1+dsqrt(dble(i))
C*2020  continue
c
C*      if(im2.eq.0)then
C*       if(imyrank.eq.0)then
C*        write(8,611)
C*       else
C*        write(1000+imyrank,611)
C*       endif
C*      endif
C*      do 2000 i=1,inp
C*       hold1=dsqrt(dble(i))
C*       hold2=dsqrt(dble(i-1))
C*       sum2=sum2+hold2
C*       wgst(i)=hold1/sum1
C*       gstt=tst3*(sum2+(hold1/2.d0))/sum1
c
C*       istop=0
C*******
C*       do while(istop.eq.0)
C*        sum=sum+(dexp(-gst*gst/tst)*gst*gst*gst*gst*gst*dgst)
C*        gst=gst+dgst
C*        if(sum.gt.gstt) pgst(i)=gst-(dgst/2.d0)
C*        istop=1
C*        if(sum.lt.gstt)istop=0
C*       enddo
C*******
C*       hold1=dsqrt((pgst(i)*pgst(i)*eo)/(0.5d0*mu))
C*       hold2=0.5d0*mu*hold1*hold1/(xk*t)
C*       hold3=dexp(-pgst(i)*pgst(i)/tst)*pgst(i)**5.d0
C*       if(im2.eq.0)then
C*        if(imyrank.eq.0)then
C*         write(8,610) pgst(i),wgst(i),hold1,hold2,
C*     1    hold3,sum/tst3
C*        else
C*         write(1000+imyrank,610) pgst(i),wgst(i),hold1,hold2,
C*     1    hold3,sum/tst3
C*        endif
C*       endif
C*2000  continue
C*************
C*******Update: use linear spacing now
C******* chose limits such that integration is accurate for
C******* all OM_ls with s=1-4
C******* pg_max corresponds to GAMMA(6,x)/120 = 1.E-03
C******* pg_min corresponds to GAMMA(3,x)/2   = 1-1.E-04
C******* where x=pgst^2/tst, GAMMA(n,x) is incomplete gamma func
      pg_max=dsqrt(16.455*tst*tt2/tt1)
      pg_min=dsqrt(0.0862*tst)
C******* make the g grid
      dgst=(pg_max-pg_min)/(dble(inp-1))
      pgst(1)=pg_min
      do i=2,inp
       pgst(i)=pgst(i-1)+dgst
      enddo
      if(imyrank.eq.0)then
       write(8,611) inp
       write(8,711) pg_min,pg_max
      else
       write(1000+imyrank,611)
       write(1000+imyrank,711)pg_min,pg_max
      endif
C*** PATCH PATCH:Removing this print
C*** printing of pgst grid
C      if(im2.eq.0)then
C       if(imyrank.eq.0)then
C        write(8,611)
C       else
C        write(1000+imyrank,611)
C       endif
C      endif
C***PATCH PATCH
      do i=1,inp
       hold1=dsqrt((pgst(i)*pgst(i)*eo)/(0.5d0*mu))
       hold2=0.5d0*mu*hold1*hold1/(xk*tt1)
       hold3=dexp(-pgst(i)*pgst(i)/tst)*pgst(i)**5.d0
       hold3=hold3/tst3
       hold4=hold3*pgst(i)**6.d0*(1.d0/(60.d0*tst3))
C***PATCH PATCH: Moving this out of the loop and reducing prints
C       if(im2.eq.0)then
C        if(imyrank.eq.0)then
C         write(8,610) pgst(i),hold1,hold2,hold3,hold4
C        else
C         write(1000+imyrank,610) pgst(i),hold1,hold2,
C     1    hold3,hold4
C        endif
C       endif
C***PATCH PATCH
      enddo
c
c
C*******re-declare all necessary variables
C******
c     determine b2max
C*****
      do ig=1,inpmax
       parab2max(ig)=0.0d0
       b2max(ig)=0.0d0
      enddo
C*****
c
      dbst2=1.d0
      dbst22=dbst2/10.d0
      cmin=0.0005d0
      if(im2.eq.0)then
       if(imyrank.eq.0)then
        write(8,652) imp,cmin
       else
        write(1000+imyrank,652) imp,cmin
       endif
      endif
c********************************************
      istart=inp-inp_per_node*imyrank
      ifinish=inp-inp_per_node*(imyrank+1)+1
C********************************************
      do 3030 ig=istart,ifinish,-1
       gst2=pgst(ig)*pgst(ig)
       v=dsqrt((gst2*eo)/(0.5d0*mu))
       ibst=nint(rmaxx/ro)-6
       if(ig.lt.istart)ibst=nint(b2max(ig+1)/dbst2)-6
       if(ibst.lt.0) ibst=0
C*****************
       istop=0
       do while(istop.eq.0)
        bst2=dbst2*dble(ibst)
        b=ro*dsqrt(bst2)
        call gsang(v,b,erat,ang,d1,istep)
        cosx(ibst)=1.d0-dcos(ang)
C      if(ip.eq.1) write(8,651) b,bst2,ang,cosx(ibst),erat
        if(ibst.ge.5)then
         if(cosx(ibst).lt.cmin.and.cosx(ibst-1).lt.cmin.and.
     1    cosx(ibst-2).lt.cmin.and.cosx(ibst-3).lt.cmin.and.
     1    cosx(ibst-4).lt.cmin) istop=1
        endif
        if(istop.eq.0)then
        ibst=ibst+1
C***Error checking
         if(ibst.gt.750) then
          if(imyrank.eq.0)then
           write(8,653)
           call flush(8)
          else
           write(1000+imyrank,653)
           call flush(1000+imyrank)
          endif
          if(imyrank.eq.0)close(8)
          call MPI_FINALIZE(ierr)
          stop     
         endif
        endif
C***Printing
       enddo
       b2max(ig)=dble(ibst-5)*dbst2
       istop=0
       do while(istop.eq.0)
        b2max(ig)=b2max(ig)+dbst22
        b=ro*dsqrt(b2max(ig))
        call gsang(v,b,erat,ang,d1,istep)
        if((1.d0-dcos(ang)).le.cmin)istop=1 
       enddo
C*****************
 3030 continue
C********************************************
C********Collect b2max values and broadcast them
C******
      call MPI_REDUCE(b2max,parab2max,inpmax,
     1 MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      do ig=1,inpmax
       b2max(ig)=parab2max(ig)
      enddo
      call MPI_BCAST(b2max,inpmax,MPI_DOUBLE_PRECISION,0,
     1 MPI_COMM_WORLD,ierr)
C*******
      if(im2.eq.0) then
       if(imyrank.eq.0)then
        write(8,637) 
        do ig=1,inp
         write(8,630) pgst(ig),b2max(ig),ro*dsqrt(b2max(ig))*1.0d10
        enddo
       else
        write(1000+imyrank,637) 
        do ig=1,inp
         write(1000+imyrank,630) 
     1    pgst(ig),b2max(ig),ro*dsqrt(b2max(ig))*1.0d10
        enddo
       endif
      endif
C**** Calculate q1st terms here
      do 4011 ig=1,inp
       q1st(ig)=0.d0
       q2st(ig)=0.d0
       q3st(ig)=0.d0
       q1STD(ig)=0.d0
       q2STD(ig)=0.d0
       q3STD(ig)=0.d0
 4011 continue
      if(imyrank.eq.0)then
       write(8,781)itn
      else
       write(1000+imyrank,781)itn
      endif
      do 4040 ic=1,itn
       if(imyrank.eq.0)then
        if(ip.eq.1) write(8,681) ic
        e_time=MPI_WTIME()
        write(8,681) ic
        write(8,*)'in',e_time-s_time,'s'
        call flush(8)
       else
        if(ip.eq.1) write(1000+imyrank,681) ic
        e_time=MPI_WTIME()
        write(1000+imyrank,681) ic
        write(1000+imyrank,*)'in',e_time-s_time,'s'
        call flush(1000+imyrank)
       endif
c
c* print Q*(l) at every iteration - header
c
c       if(imyrank.eq.0)then
c         write(8,679)
c       else
c         write(1000+imyrank,679)
c       endif
c
      do 4010 ig=1,inp
        gst2=pgst(ig)*pgst(ig)
        v=dsqrt((gst2*eo)/(0.5d0*mu))
C        if(imyrank.eq.0)then
C         if(ip.eq.1) write(8,682) ic,ig,gst2,v
C        endif
        temp1=0.0d0
        temp2=0.0d0
        temp3=0.0d0
        mpitemp1=0.0d0
        mpitemp2=0.0d0
        mpitemp3=0.0d0
        do 4000 im=1,imp_per_node
         rnb=xrand()
         call rantate
         bst2=rnb*b2max(ig)
         b=ro*dsqrt(bst2)
         call gsang(v,b,erat,ang,d1,istep)
         hold1=1.d0-dcos(ang)
         hold2=dsin(ang)*dsin(ang)
         hold3=1.d0-(dcos(ang))**3
         temp1=temp1
     1    +(hold1*b2max(ig)/dble(imp))
         temp2=temp2
     1    +(1.5d0*hold2*b2max(ig)/dble(imp))
         temp3=temp3
     1    +(hold3*b2max(ig)/dble(imp))
4000    continue
c*AH
        call MPI_REDUCE(temp1,mpitemp1,1,
     1   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(temp2,mpitemp2,1,
     1   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(temp3,mpitemp3,1,
     1   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
c* 231004
        call MPI_BCAST(mpitemp1,1,MPI_DOUBLE_PRECISION,0,
     1 MPI_COMM_WORLD,ierr)
        call MPI_BCAST(mpitemp2,1,MPI_DOUBLE_PRECISION,0,
     1 MPI_COMM_WORLD,ierr)
        call MPI_BCAST(mpitemp3,1,MPI_DOUBLE_PRECISION,0,
     1 MPI_COMM_WORLD,ierr)
c* 231004 end
c
C* printing of Q*(l) for each loop of itn
c
c       if(imyrank.eq.0)then
c         write(8,676) pgst(ig),mpitemp1,mpitemp2,mpitemp3
c       else
c         write(1000+imyrank,676) pgst(ig),mpitemp1,mpitemp2,mpitemp3
c       endif
c
c collect values of all cores into 2D array only for node 0
        if(imyrank.eq.0)then
         q12D(ig,ic)=mpitemp1
         q22D(ig,ic)=mpitemp2
         q32D(ig,ic)=mpitemp3
        endif
4010  continue
4040  continue
C*****
      if(imyrank.eq.0)then
       do ig=1,inp
c
c     MEAN and STD of qst
        do ic=1,itn
         q1st(ig)=q1st(ig)+q12D(ig,ic)/dble(itn)
         q2st(ig)=q2st(ig)+q22D(ig,ic)/dble(itn)
         q3st(ig)=q3st(ig)+q32D(ig,ic)/dble(itn)
        enddo
        STD1sum=0.0d0
        STD2sum=0.0d0
        STD3sum=0.0d0
        do ic=1,itn
         STD1sum=STD1sum+(q1st(ig)-q12D(ig,ic))**2
         STD2sum=STD2sum+(q2st(ig)-q22D(ig,ic))**2
         STD3sum=STD3sum+(q3st(ig)-q32D(ig,ic))**2
        enddo
C* q1STD is confidence-invervall, hence sig/sqrt(itn)*(99% factor)
C* for the 99% factor, we choose a Gaussian distribution
        q1STD(ig)=dsqrt(STD1sum/dble(itn*(itn-1)))*2.576d0
        q2STD(ig)=dsqrt(STD2sum/dble(itn*(itn-1)))*2.576d0
        q3STD(ig)=dsqrt(STD3sum/dble(itn*(itn-1)))*2.576d0
       enddo
      endif
      call MPI_BCAST(q1st,inpmax,MPI_DOUBLE_PRECISION,0,
     1 MPI_COMM_WORLD,ierr)
      call MPI_BCAST(q2st,inpmax,MPI_DOUBLE_PRECISION,0,
     1 MPI_COMM_WORLD,ierr)
      call MPI_BCAST(q3st,inpmax,MPI_DOUBLE_PRECISION,0,
     1 MPI_COMM_WORLD,ierr)
      call MPI_BCAST(q1STD,inpmax,MPI_DOUBLE_PRECISION,0,
     1 MPI_COMM_WORLD,ierr)
      call MPI_BCAST(q2STD,inpmax,MPI_DOUBLE_PRECISION,0,
     1 MPI_COMM_WORLD,ierr)
      call MPI_BCAST(q3STD,inpmax,MPI_DOUBLE_PRECISION,0,
     1 MPI_COMM_WORLD,ierr)
c
      if(im2.eq.0)then
       if(imyrank.eq.0)then
        write(8,672) itn,inp,imp,itn*inp*imp
        write(8,675)
        do ig=1,inp
         write(8,676) pgst(ig),q1st(ig),q1STD(ig),
     1    q2st(ig),q2STD(ig),q3st(ig),q3STD(ig)
        enddo
       else
        write(1000+imyrank,672) itn,inp,imp,itn*inp*imp
        write(1000+imyrank,675)
        do ig=1,inp
         write(1000+imyrank,676) pgst(ig),q1st(ig),q1STD(ig),
     1    q2st(ig),q2STD(ig),q3st(ig),q3STD(ig)
        enddo
       endif
      endif
c
C***********************************
673   format(//1x,'f value for second order correction=',1pe11.4,/1x
     1 '(see Hirschfelder, Curtis, Bird, 1954, p. 606)')
674   format(//1x,'mean OMEGA*(1,1) =',1pe11.4,/1x,
     1 'standard deviation =',
     1 e11.4,/1x,'standard error of mean =',e11.4)
C*676   format(1x,1pe11.4,1x,e11.4,1x,e11.4)
C*C676   format(1x,1pe11.4,3(1x,e11.4))
676   format(1x,1pe11.4,3(1x,e11.4,' +/-',e9.2))
675   format(//1x,'Final averaged values of Q*(l): ',//5x,
     1 'gst ',15x,'q1st',21x,'q2st',21x,'q3st')
679   format(/1x,'Current values of Q*(l): ',//5x,
     1 'gst ',8x,'q1st',8x,'q2st',8x,'q3st')
622   format(1x,i3,4x,1pe11.4,4x,e11.4,4x,e11.4,4x,e11.4)
685   format(/1x,'summary of mobility calculations',//1x,'cycle',
     1 5x,'cs/A^2',6x,'avge cs/A^2',8x,'Ko^-1',7x,'avge Ko^-1')
620   format(/1x,'OMEGA(1,1)*=',1pe11.4,/)
670   format(/1x,'v =',1pe11.4,5x,'q1st =',e11.4,/)
684   format(1x,1pe11.4,7(e11.4))
683   format(/5x,'b/A',8x,'ang',6x,'(1-cosX)',4x,'e ratio',4x,'theta',
     1 7x,'phi',7x,'gamma')
682   format(/1x,'ic =',i3,1x,'ig =',i4,1x,'gst2 =',1pe11.4,
     1 1x,'v =',e11.4)
680   format(1x,'start mobility calculation')
651   format(1x,1pe11.4,6(1x,e11.4))
653   format(1x,'ibst greater than 750')
637   format(/5x,'gst',13x,'b2max/ro2',9x,'b/A')
630   format(1x,1pe13.6,5x,e13.6,5x,e13.6)
672   format(//1x,'number of complete cycles (itn) = ',i7,/1x,
     1 'number of velocity points (inp) = ',i7,/1x,
     1 'number of random points   (imp) = ',i7,/1x,
     1 'total number of points          = ',i7,/)
781   format(//1x,'Starting iterations over itn = ',i3,' cycles',
     1 /1x,41('-'))
681   format(//1x,'Entering cycle number, ic =',i3)
650   format(/1x,'gst2 =',1pe11.4,1x,'v =',e11.4,/6x,'b',
     1 10x,'bst2',7x,'X ang',7x,'cos(X)',6x,'e ratio')
652   format(/1x,'Set up Monte-Carlo integration:',/1x,'imp=',
     1 i4,' random points over impact parameter b and orientation',
     1 /1x,'determine b2max value for each gst grid point using',
     1 /1x,'minimum value of (1-cosX) =',1pe11.4)
C*610   format(1x,1pe11.4,5(1x,e11.4))
C*611   format(//1x,'set-up gst integration - integration over velocity',
C*     1 //5x,'pgst',8x,'wgst',9x,'v',9x,'ke/kt',7x,'gst^5*',5x,'frac of',
C*     1 /48x,'exp(gst^2/tst)',3x,'sum',/)
610   format(1x,1pe11.4,2(2x,e11.4),3x,e11.4,5x,e11.4)
611   format(/1x,'Set up velocity integration:',/1x,
     1 'inp=',i3,' gst points on regular grid between')
711   format(1x,'gst_min =',e11.4,' and gst_max =',e11.4)
600   format(/1x,'t*=',1pe11.4)
615   format(1x,'along z axis emax =',1pe11.4,'eV rmax =',
     1 e11.4,'A r00 =',e11.4,'A',/)
613   format(1x,'along y axis emax =',1pe11.4,'eV rmax =',
     1 e11.4,'A r00 =',e11.4,'A')
614   format(1x,'along x axis emax =',1pe11.4,'eV rmax =',
     1 e11.4,'A r00 =',e11.4,'A')
601   format(/1x,'Problem orientating along x axis',/)
632   format(/1x,'maximum extent orientated along x axis')
631   format(/1x,82('-'),/1x,'Calculating trajectories to obtain
     1 momentum transfer cross sections Q*(l), l=1,2,3',/1x,82('-'),
     1 //1x,'Set-ups',/1x,7('-'),/)
603   format(1x,'global trajectory parameters',/1x,'sw1 =',1pe11.4,7x,
     1 'sw2 =',e11.4,/1x,'dtsf1 =',e11.4,5x,'dtsf2 =',e11.4,/1x,
     1 'inwr =',i3,14x,'ifail =',i5)
602   format(1x,i4,5x,1pe11.4,3(5x,e11.4))
689   format(/)
      return
      end
c
c     ***************************************************************
c
      subroutine mat_el (tion,OM11,OM12,OM13,OM14,OM22,OM23,OM24,a_2T)
c
c     Subroutine to get the matrix elements a_rs(l) within 
c     2TT at current Teff. see Table 6-2-1 of Mason&McDaniel 1988
c
      use mpi
      implicit double precision (a-h,m-z)
      parameter (ixlen=1000,inpmax=201)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xeo,xk,xn,xmv,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(ixlen),k_const(ixlen),tcharge
      common/coordinates/fx(ixlen),fy(ixlen),fz(ixlen),
     ?ox(ixlen),oy(ixlen),oz(ixlen)
      common/ljparameters/eolj(ixlen),rolj(ixlen),eox4(ixlen),
     ?ro6lj(ixlen),ro12lj(ixlen),dro6(ixlen),dro12(ixlen)
      common/ljparameters2/alphai(ixlen),Ni(ixlen),Ai(ixlen),Gi(ixlen),
     1 RijStar(ixlen),eij(ixlen)
      common/mmffN2/alphaN2,NiN2,AiN2,GiN2,RN2Star,MMFF_B,MMFF_beta
      common/mmffHe/alphaHe,NiHe,AiHe,GiHe,RHeStar
      common/hsparameters/rhs(ixlen),rhs2(ixlen)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/temperature/tt1,tt2,ittgrid
      common/xrandom/i1,i2,i3,i4,i5,i6
      common/scaling/enersca,distsca
C****
      common/mobilcal/b2max(inpmax),parab2max(inpmax)
      common/mobilgrid/pgst(inpmax),wgst(inpmax),q1st(inpmax),
     1  q2st(inpmax),q3st(inpmax),q1STD(inpmax),q2STD(inpmax),
     1  q3STD(inpmax)
      common/mpicon/imyrank,inprocs,imp_per_node,inp_per_node,s_time
      common/runparams/itn,inp,imp,igas
      dimension a_2T(0:2,0:2,0:2)
c
       do i1=0,2
        do i2=0,2
         do i3=0,2
          a_2T(i1,i2,i3)=0.0d0
         enddo
        enddo
       enddo
c
       e_2T=m1/(m1+m2)
       d_2T=m1*tion/(m2*tt1+m1*tion)
       f_2T=d_2T*(1.0d0-e_2T)*(tion-tt1)/tion
       Ast=OM22/OM11
       Bst=(5.0d0*OM12-4.0d0*OM13)/OM11
       Cst=OM12/OM11
       Dst=OM14/OM11
       Est=OM23/OM22
       Gst=OM24/OM22
c
       Ahat=e_2T*Ast
       Bhat=4.0d0*Bst-5.0d0
       Chat=6.0d0*Cst-5.0d0
       Dhat=32.0d0*Dst-21.0d0
       Ehat=8.0d0*Est-7.0d0
       Ghat=80.0d0*Gst-63.0d0
       dsq=d_2T**2
       esq=e_2T**2
       fsq=f_2T**2
       e2_3f2=esq+3.0d0*fsq
       e2_f2=esq+fsq
c
       a_2T(0,0,0)=0.0d0
       a_2T(0,1,0)=0.0d0
       a_2T(0,2,0)=0.0d0
       a_2T(1,0,0)=-8.0d0*f_2T/d_2T
       a_2T(1,1,0)=(8.0/3.0)*(2.0*(1.0-e_2T)+f_2T*Chat)
       a_2T(1,2,0)=-(8.0*d_2T/15.0)*(4.0*(1.0-d_2T)*Chat-3.0*f_2T*Bhat)
       a_2T(2,0,0)=-(4.0/dsq)*(10.0*e_2T*fsq-f_2T*e2_f2*Chat
     1             -4.0*fsq*Ahat)
       a_2T(2,1,0)=-(4.0/(3.0*d_2T))*(10.0*f_2T*(1.0-4.0*e_2T+4.0*esq)
     1          +2.0*(dsq*(1.0-d_2T)-f_2T*d_2T*(2.0-3.0*d_2T+11.0*f_2T)
     1             +fsq*(4.0+7.0*f_2T))*Chat-3.0*f_2T*e2_f2*Bhat
     1             +16.0*f_2T*(1.0-e_2T)*Ahat+4.0*fsq*d_2T*Ehat*Ahat)
       a_2T(2,2,0)=(4.0/15.0)*(40.0*(1.0-e_2T)*(1.0-2.0*e_2T+2.0*esq)
     1             +f_2T*(44.0-128.0*d_2T+80.0*f_2T-dsq+90.0*d_2T*f_2T
     1             -170.0*fsq)*Chat+3.0*(5.0*f_2T*(dsq+2.0*esq+2.0*fsq)
     1             -4.0*(1.0-d_2T)*e2_3f2)*Bhat+15.0*f_2T*e2_f2*Dhat
     1             +32.0*((1.0-e_2T)**2)*Ahat+8.0*f_2T*(4.0*(1.0-d_2T)
     1             -5.0*f_2T)*Ehat*Ahat+4.0*fsq*Ghat*Ahat)
       a_2T(0,0,1)=8.0/3.0
       a_2T(0,1,1)=-(8.0*d_2T/15.0)*Chat
       a_2T(0,2,1)=-(8.0*dsq/105.0)*(4.0*Chat+3.0*Bhat)
       a_2T(1,0,1)=-(4.0/(3.0*d_2T))*(10.0*f_2T*(1.0-2.0*e_2T)
     1             +e2_3f2*Chat+8.0*f_2T*Ahat)
       a_2T(1,1,1)=(4.0/15.0)*(10.0*(3.0-6.0*e_2T+4.0*esq)
     1             +(5.0*e2_3f2-5.0*dsq+22.0*f_2T*(1.0-d_2T))*Chat
     1             -3.0*e2_3f2*Bhat+16.0*(1.0-e_2T)*Ahat
     1             +8.0*f_2T*Ehat*Ahat)
       a_2T(0,0,2)=(16.0/15.0)*(5.0*(1.0-e_2T)+f_2T*Chat+3.0*Ahat)
       a_2T(0,1,2)=(16.0*d_2T/105.0)*((7.0*(d_2T-1.0)-3.0*f_2T)*Chat
     1             +3.0*f_2T*Bhat-3.0*Ehat*Ahat)
c
c
      return
      end
c
c     ***************************************************************
c
      subroutine trunc_it (g_2T,vDD,EEps)
c
c     Subroutine to calculate drift velocity and field strength within
c     2TT using matrix elements gamma(Teff). see Viehland&Mason1978
c
      use mpi
      implicit double precision (a-h,m-z)
      parameter (ixlen=1000,inpmax=201)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xeo,xk,xn,xmv,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(ixlen),k_const(ixlen),tcharge
      common/coordinates/fx(ixlen),fy(ixlen),fz(ixlen),
     ?ox(ixlen),oy(ixlen),oz(ixlen)
      common/ljparameters/eolj(ixlen),rolj(ixlen),eox4(ixlen),
     ?ro6lj(ixlen),ro12lj(ixlen),dro6(ixlen),dro12(ixlen)
      common/ljparameters2/alphai(ixlen),Ni(ixlen),Ai(ixlen),Gi(ixlen),
     1 RijStar(ixlen),eij(ixlen)
      common/mmffN2/alphaN2,NiN2,AiN2,GiN2,RN2Star,MMFF_B,MMFF_beta
      common/mmffHe/alphaHe,NiHe,AiHe,GiHe,RHeStar
      common/hsparameters/rhs(ixlen),rhs2(ixlen)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/temperature/tt1,tt2,ittgrid
      common/xrandom/i1,i2,i3,i4,i5,i6
      common/scaling/enersca,distsca
C****
      common/mobilcal/b2max(inpmax),parab2max(inpmax)
      common/mobilgrid/pgst(inpmax),wgst(inpmax),q1st(inpmax),
     1  q2st(inpmax),q3st(inpmax),q1STD(inpmax),q2STD(inpmax),
     1  q3STD(inpmax)
      common/mpicon/imyrank,inprocs,imp_per_node,inp_per_node,s_time
      common/runparams/itn,inp,imp,igas
      dimension psi(0:2,0:2,1:2),Eps(1:2),g_2T(0:2,0:2,0:2)
c
      iNapx=2
      vDD=0.0d0
      EEps=0.0d0
      do 4001 inx=1,iNapx
c fix moments
       psi(0,0,inx)=1.0d0
       psi(0,1,inx)=0.0d0
c step 2: Eps(inx), psi(1,0)_inx and vD_inx
       ssum1=0.0d0
       do is=1,(inx-1)
        ssum1=ssum1+g_2T(0,is,1)*psi(1,is,inx-is)
       enddo
       ssum2=0.0d0
       do is=2,inx
        ssum2=ssum2+g_2T(1,is,0)*psi(0,is,inx+1-is)
       enddo
       Epsi=ssum1+dsqrt((ssum1**2)-2.0d0*(g_2T(1,0,0)+ssum2))
       Eps(inx)=Epsi/2.0d0
       ssum3=0.0d0
       do is=1,inx
        ssum3=ssum3+g_2T(0,is,1)*psi(1,is,inx-is)
       enddo
       psi(1,0,inx)=(Eps(inx)-ssum3)/g_2T(0,0,1)
c check if done
       if(inx.eq.iNapx)then
        vDD=psi(1,0,inx)
        EEps=Eps(inx)
        goto 4002
       endif
c step 3: more moments in current order of approx
c    l,r=2,0
       term1=Eps(inx)*(2.0*psi(2-1,0,inx))
       term2=0.0d0
       term0=0.0d0
       psi(2,0,inx)=(term1-term2-term0)/g_2T(0,0,2)
c    l,r=1,1
       term1=Eps(inx)*(5.0*psi(1-1,1,inx)-4.0*psi(1+1,1-1,inx))/3.0
       term2=0.0d0
       term0=g_2T(1,0,1)*psi(1,0,inx)
       psi(1,1,inx)=(term1-term2-term0)/g_2T(1,1,1)
c    l,r=0,2
       term1=Eps(inx)*(2.0*psi(0+1,2-1,inx))
       term2=0.0d0
       term0=g_2T(2,0,0)*psi(0,0,inx)+g_2T(2,1,0)*psi(0,1,inx)
       psi(0,2,inx)=(term1-term2-term0)/g_2T(2,2,0)
4001   continue
4002   continue
c
c
      return
      end
c
c     ***************************************************************
c
      subroutine mobil2 (t,mob,cs,sdevpc,EN)
c
c     Subroutine to determine average mobility by trajectory method. 
c     All integrations Monte Carlo (except over velocity).
c
      use mpi
      implicit double precision (a-h,m-z)
      parameter (ixlen=1000,inpmax=201)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xeo,xk,xn,xmv,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(ixlen),k_const(ixlen),tcharge
      common/coordinates/fx(ixlen),fy(ixlen),fz(ixlen),
     ?ox(ixlen),oy(ixlen),oz(ixlen)
      common/ljparameters/eolj(ixlen),rolj(ixlen),eox4(ixlen),
     ?ro6lj(ixlen),ro12lj(ixlen),dro6(ixlen),dro12(ixlen)
      common/ljparameters2/alphai(ixlen),Ni(ixlen),Ai(ixlen),Gi(ixlen),
     1 RijStar(ixlen),eij(ixlen)
      common/mmffN2/alphaN2,NiN2,AiN2,GiN2,RN2Star,MMFF_B,MMFF_beta
      common/mmffHe/alphaHe,NiHe,AiHe,GiHe,RHeStar
      common/hsparameters/rhs(ixlen),rhs2(ixlen)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/temperature/tt1,tt2,ittgrid
      common/xrandom/i1,i2,i3,i4,i5,i6
      common/scaling/enersca,distsca
      common/empcorr/Aec,Bec
C****
      common/mobilcal/b2max(inpmax),parab2max(inpmax)
      common/mobilgrid/pgst(inpmax),wgst(inpmax),q1st(inpmax),
     1  q2st(inpmax),q3st(inpmax),q1STD(inpmax),q2STD(inpmax),
     1  q3STD(inpmax)
      common/mpicon/imyrank,inprocs,imp_per_node,inp_per_node,s_time
      common/runparams/itn,inp,imp,igas
      dimension a_2T(0:2,0:2,0:2),g_2T(0:2,0:2,0:2)
c     
       do i1=0,2
        do i2=0,2
         do i3=0,2
          a_2T(i1,i2,i3)=0.0d0
          g_2T(i1,i2,i3)=0.0d0
         enddo
        enddo
       enddo
C****
c     Calculate Omega(1,1)*, Omega(1,2)*, Omega(1,3)*, and Omega(2,2)*
c     by integrating Q(1)* or Q(2)* over all orientations, and initial 
c     relative velocities.
c
C** Define reduced temperature based upon input T for calculation
      tst=xk*t/eo
      if(im2.eq.0)then 
       if(imyrank.eq.0)then
        write(8,642) t,tst
       else
        write(1000+imyrank,642) t,tst
       endif
      endif
      tst3=tst*tst*tst
C**** Going to remove (ic) from these:
       mom11st=0.d0
       mom12st=0.d0
       mom13st=0.d0
       mom14st=0.d0
       mom22st=0.d0
       mom23st=0.d0
       mom24st=0.d0
       mom33st=0.d0
       mom34st=0.d0
C** their variances
       VARmom11st=0.d0
       VARmom12st=0.d0
       VARmom13st=0.d0
       VARmom14st=0.d0
       VARmom22st=0.d0
       VARmom23st=0.d0
       VARmom24st=0.d0
       VARmom33st=0.d0
       VARmom34st=0.d0
       pg_max=dsqrt(16.455*xk*tt2/eo)
       pg_min=dsqrt(0.0862*xk*tt1/eo)
       dgst=(pg_max-pg_min)/(dble(inp-1))
C******
       do 4010 ig=1,inp
        gst2=pgst(ig)*pgst(ig)
        v=dsqrt((gst2*eo)/(0.5d0*mu))
C        if(imyrank.eq.0)then
C         if(ip.eq.1) write(8,682) ic,ig,gst2,v
C        endif
C***** Need to implement new weighting functions for linear
C***** included additional omegas (up to l=3 and s=4)
        exp_fac=dexp(-pgst(ig)**2/tst)*(pgst(ig)**5)*dgst/(tst**3)
        wgst(ig)=exp_fac/dgst
        s1_fac=1.d0*exp_fac
        s2_fac=(pgst(ig)**2)*exp_fac*(1.d0/(3.d0*tst))
        s3_fac=(pgst(ig)**4)*exp_fac*(1.d0/(12.d0*tst*tst))
        s4_fac=(pgst(ig)**6)*exp_fac*(1.d0/(60.d0*tst*tst*tst))
C***Removing (ic) from these:
        mom11st=mom11st+(q1st(ig)*s1_fac)
        mom12st=mom12st+(q1st(ig)*s2_fac)
        mom13st=mom13st+(q1st(ig)*s3_fac)
        mom14st=mom14st+(q1st(ig)*s4_fac)
        mom22st=mom22st+(q2st(ig)*s2_fac)
        mom23st=mom23st+(q2st(ig)*s3_fac)
        mom24st=mom24st+(q2st(ig)*s4_fac)
        mom33st=mom33st+(q3st(ig)*s3_fac)
        mom34st=mom34st+(q3st(ig)*s4_fac)
C** their VARs
        VARmom11st=VARmom11st+((q1STD(ig)*s1_fac)**2)
        VARmom12st=VARmom12st+((q1STD(ig)*s2_fac)**2)
        VARmom13st=VARmom13st+((q1STD(ig)*s3_fac)**2)
        VARmom14st=VARmom14st+((q1STD(ig)*s4_fac)**2)
        VARmom22st=VARmom22st+((q2STD(ig)*s2_fac)**2)
        VARmom23st=VARmom23st+((q2STD(ig)*s3_fac)**2)
        VARmom24st=VARmom24st+((q2STD(ig)*s4_fac)**2)
        VARmom33st=VARmom33st+((q3STD(ig)*s3_fac)**2)
        VARmom34st=VARmom34st+((q3STD(ig)*s4_fac)**2)
4010   continue
C*
C*  printing Omega values
       if(im2.eq.0) write(8,677) mom11st,mom12st,mom13st,mom14st,
     1  mom22st,mom23st,mom24st,mom33st,mom34st,pi*ro*ro*1.d20
C*** Error of CCS
       sterr=dsqrt(VARmom11st)
c* new way is already std err.
       cs=mom11st*pi*ro*ro
       sterr2=sterr*pi*ro*ro
       sdevpc=100.d0*sterr/mom11st
C*
C*   Use 2TT to obtain mobility at set field strength
C*   tt1 is T_bath and t=T_eff -> hence, we have OM(l,s) at Teff!
C*   m1 is m_N2 and m2=m_ion
C*
       T_ion=((m1+m2)*t-m2*tt1)/m1
c
       call mat_el(T_ion,mom11st,mom12st,mom13st,mom14st,mom22st,
     1             mom23st,mom24st,a_2T)
c
c normalize the matrix elements to a_00(1) = a_2T(0,0,1)
       do i1=0,2
        do i2=0,2
         do i3=0,2
          g_2T(i1,i2,i3)=a_2T(i1,i2,i3)/a_2T(0,0,1)
         enddo
        enddo
       enddo
c
c perform 2TT truncation iteration procedure to 2nd order
       call trunc_it(g_2T,vDD,EEps)
c
c transform drift velocity and field into usable units
       EN=EEps*16.0d0/(3.0d0*dsqrt(pi)*xe*dabs(tcharge))
       EN=EN*dsqrt(m1/(m1+m2))*xk*dsqrt(T_ion*t)*cs
       vD=vDD*dsqrt(2.0d0*xk*T_ion/(m2/(xn*1.0d3)))
c
c calcualte mobility
        mob_1st=mconst/(dsqrt(t)*cs)
       if(t.eq.tt1)then
        EN=0.0d0
        vD=0.0d0
        mob_2nd=mob_1st
       else
        dens=xn/xmv
        mob_2nd=vD/(EN*dens)
       endif
       mob=mob_2nd
       sig_mob1=mob_1st/cs*sterr2
       sig_mob2=mob_2nd/cs*sterr2
c empirical correction
       empcorr=1.0d0+Aec*dexp(-Bec*1.d-21/EN)
       mob=mob*empcorr
       sig_mob=sig_mob2*empcorr
c printout
       if(correct.gt.0)then
        write(8,671) EN*1.d21,vD,mob_1st,sig_mob1,mob_2nd,sig_mob2,
     1               mob,sig_mob,cs*1.d20,sterr2*1.d20
       else
        write(8,672) EN*1.d21,vD,mob_1st,sig_mob1,mob_2nd,sig_mob2,
     1               cs*1.d20,sterr2*1.d20
       endif
c
c
C******************* old 2nd order correction *******
c
c     Use omegas to obtain higher order correction factor to mobility
c
C*C       ayst=mom22st/mom11st
C*C       best=((5.d0*mom12st)-(4.d0*mom13st))/mom11st
C*C       cest=mom12st/mom11st
C*C       term=((4.d0*ayst)/(15.d0))+(.5d0*((m2-m1)**2.d0)/(m1*m2))
C*C       u2=term-(.08333d0*(2.4d0*best+1.d0)*(m1/m2))
C*C       w=(m1/m2)
C*C       delta=((((6.d0*cest)-5.d0)**2.d0)*w)/(60.d0*(1.d0+u2))
C*C       f=1.d0/(1.d0-delta)
C*C       if(im2.eq.0) write(8,677) mom11st,mom12st,mom13st,mom14st,
C*C     1  mom22st,mom23st,mom24st,mom33st,mom34st,pi*ro*ro*1.d20
C*C       if(im2.eq.0) write(8,673) f,u2,w,delta
C*C       mob=(mconst*f)/(dsqrt(t)*cs)
C*C       write(8,671) mob,1.d0/mob,cs*1.d20,sterr2*1.d20
C*C       if(imyrank.eq.0)then
C*C        e_time=MPI_WTIME()
C*C        write(8,*)'Job Completed in',e_time-s_time,'s'
C*C        call flush(8)
C*C       endif
C***********************************
       call flush(8)
C***********************************
C      endif
C***********************************
C***********************************
671   format(/1x,'red. field stregth E/N   =  ',f8.3,' Td'
     1 /1x,     'drift velocity vD        =',f8.1,'   m/s'
     1 /1x,     'red. mobility (1st ord.) =',e11.4,' +/-',e9.2,' m2/Vs'
     1 /1x,     'red. mobility (2nd ord.) =',e11.4,' +/-',e9.2,' m2/Vs'
     1 /1x,     'red. mobility (emp corr) =',e11.4,' +/-',e9.2,' m2/Vs'
     1 /1x      'av. coll cross section   =',e11.4,' +/-',e9.2,' A**2')
672   format(/1x,'red. field stregth E/N   =  ',f8.3,' Td'
     1 /1x,     'drift velocity vD        =',f8.1,'   m/s'
     1 /1x,     'red. mobility (1st ord.) =',e11.4,' +/-',e9.2,' m2/Vs'
     1 /1x,     'red. mobility (2nd ord.) =',e11.4,' +/-',e9.2,' m2/Vs'
     1 /1x      'av. coll cross section   =',e11.4,' +/-',e9.2,' A**2')
677   format(1x,'OM*(l,s)',5x,'s=1',9x,'s=2',9x,'s=3',9x,'s=4',
     1 /1x,'   l=1  ',1x,1pe11.4,1x,1pe11.4,1x,1pe11.4,1x,1pe11.4,
     1 /1x,'   l=2  ',12x       ,1x,1pe11.4,1x,1pe11.4,1x,1pe11.4,
     1 /1x,'   l=3  ',24x                  ,1x,1pe11.4,1x,1pe11.4,
     1 /1x,'(multiply by ',e11.4,' for values in Angstrom**2)')
c
673   format(//1x,'f value for second order correction =',1pe11.4,/,
     1 ' u2 =',e11.4,2x,'w =',e11.4,2x,'delta =',e11.4,/,
     1 ' (see Hirschfelder, Curtis, Bird, 1954, p. 606)')
642   format(//1x,'Collision integrals OM*(l,s) at Teff=',f7.2,'K',
     1 ' (Tst=',1pe11.4,')',/1x,63('-'))
674   format(/1x,'mean OMEGA*(1,1) =',1pe11.4,/1x,
     1 'standard deviation =',
     1 e11.4,/1x,'standard error of mean =',e11.4)
C*676   format(1x,1pe11.4,1x,e11.4,1x,e11.4)
676   format(2x,1pe11.4,4(1x,e11.4))
675   format(//1x,'Weighting functions w(s) for Omega(l,s)',//5x,
C*     1 'gst2',8x,'wgst',8x,'q1st')
     1 'gst',8x,'w(s=1)',6x,'w(s=2)',6x,'w(s=3)',6x,'w(s=4)')
622   format(1x,i3,4x,1pe11.4,4x,e11.4,4x,e11.4,4x,e11.4)
685   format(/1x,'summary of mobility calculations',//1x,'cycle',
     1 5x,'cs/A^2',6x,'avge cs/A^2',8x,'Ko^-1',7x,'avge Ko^-1')
620   format(/1x,'OMEGA(1,1)*=',1pe11.4,/)
670   format(/1x,'v =',1pe11.4,5x,'q1st =',e11.4,/)
684   format(1x,1pe11.4,7(e11.4))
682   format(/1x,'ic =',i3,1x,'ig =',i4,1x,'gst2 =',1pe11.4,
     1 1x,'v =',e11.4)
680   format(1x,'start mobility calculation')
651   format(1x,1pe11.4,6(1x,e11.4))
653   format(1x,'ibst greater than 750')
637   format(/5x,'gst',13x,'b2max/ro2',9x,'b/A')
630   format(1x,1pe11.4,5x,e11.4,5x,e11.4)
681   format(//1x,'Entering cycle number, ic =',i3)
C*610   format(1x,1pe11.4,5(1x,e11.4))
C*611   format(//1x,'set-up gst integration - integration over velocity',
C*     1 //5x,'pgst',8x,'wgst',9x,'v',9x,'ke/kt',7x,'gst^5*',5x,'frac of',
C*     1 /48x,'exp(gst^2/tst)',3x,'sum',/)
610   format(1x,1pe11.4,2(2x,e11.4),3x,e11.4,5x,e11.4)
611   format(//1x,'set-up gst integration - integration over velocity',
     1 //5x,'pgst',9x,'v',12x,'ke/kt',6x,'1/N gst^05*',5x,'1/N gst^11*',
     1 /41x,'exp(-gst^2/tst)',1x,'exp(-gst^2/tst)',/)
615   format(1x,'along z axis emax =',1pe11.4,'eV rmax =',
     1 e11.4,'A r00 =',e11.4,'A',/)
613   format(1x,'along y axis emax =',1pe11.4,'eV rmax =',
     1 e11.4,'A r00 =',e11.4,'A')
614   format(1x,'along x axis emax =',1pe11.4,'eV rmax =',
     1 e11.4,'A r00 =',e11.4,'A')
601   format(/1x,'Problem orientating along x axis',/)
632   format(/1x,'maximum extent orientated along x axis')
631   format(/1x,'mobility calculation by MOBIL2 (trajectory method)',/)
603   format(1x,'global trajectory parameters',/1x,'sw1 =',1pe11.4,7x,
     1 'sw2 =',e11.4,/1x,'dtsf1 =',e11.4,5x,'dtsf2 =',e11.4,/1x,
     1 'inwr =',i3,14x,'ifail =',i5)
602   format(1x,i4,5x,1pe11.4,3(5x,e11.4))
689   format(/)
      return
      end
c
c     ***************************************************************
c
      double precision function xrand()
      implicit double precision (a-h,m-z)
      dimension rvec(10)
      common/xrandom/i1,i2,i3,i4,i5,i6
c
c     XRAND is a random number generator that uses RANLUX if i5=1
c     otherwise it uses the standard RAND subroutine available in 
c     FORTRAN 77. If RANLUX is used i1 contains the luxury level
c     (1-4, 4 is the highest, 3 is default in RANLUX). i2, i3, and
c     i4 are seed integers. i3 and i4 will normally be zero. If the
c     standard RAND subroutine is to be employed, i2 contains the 
c     seed integer. RANLUX was downloaded from http://kumo.swcp.
c     com/fortran/random2.f90. 
c 
      i6=i6+1
c
      if (i5.eq.1) then
      call ranlux(rvec,1)
      xrand=rvec(1)
      return
      else
      xrand=rand()
      return
      endif
c
      end
c
c     ***************************************************************
c

      SUBROUTINE RANLUX(RVEC,LENV)
C         Subtract-and-borrow random number generator proposed by
C         Marsaglia and Zaman, implemented by F. James with the name
C         RCARRY in 1991, and later improved by Martin Luescher
C         in 1993 to produce "Luxury Pseudorandom Numbers".
C     Fortran 77 coded by F. James, 1993
C          
C       references:
C  M. Luscher, Computer Physics Communications  79 (1994) 100
C  F. James, Computer Physics Communications 79 (1994) 111
C
C   LUXURY LEVELS.
C   ------ ------      The available luxury levels are:
C
C  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
C           and Zaman, very long period, but fails many tests.
C  level 1  (p=48): considerable improvement in quality over level 0,
C           now passes the gap test, but still fails spectral test.
C  level 2  (p=97): passes all known tests, but theoretically still
C           defective.
C  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
C           correlations have very small chance of being observed.
C  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
C
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for RANLUX:                                  ++
C!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
C!!!                   32-bit random floating point numbers between  ++
C!!!                   zero (not included) and one (also not incl.). ++
C!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
C!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
C!!!               which is integer between zero and MAXLEV, or if   ++
C!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++
C!!!               should be set to zero unless restarting at a break++
C!!!               point given by output of RLUXAT (see RLUXAT).     ++
C!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
C!!!               which can be used to restart the RANLUX generator ++
C!!!               at the current point by calling RLUXGO.  K1 and K2++
C!!!               specify how many numbers were generated since the ++
C!!!               initialization with LUX and INT.  The restarting  ++
C!!!               skips over  K1+K2*E9   numbers, so it can be long.++
C!!!   A more efficient but less convenient way of restarting is by: ++
C!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
C!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
C!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
C!!!                 32-bit integer seeds, to be used for restarting ++
C!!!      ISVEC must be dimensioned 25 in the calling program        ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit double precision (a-h,o-z)
      DIMENSION RVEC(LENV)
      DIMENSION SEEDS(24), ISEEDS(24), ISDEXT(25)
      PARAMETER (MAXLEV=4, LXDFLT=3)
      DIMENSION NDSKIP(0:MAXLEV)
      DIMENSION NEXT(24)
      PARAMETER (TWOP12=4096., IGIGA=1000000000,JSDFLT=314159265)
      PARAMETER (ITWO24=2**24, ICONS=2147483563)
      SAVE NOTYET, I24, J24, CARRY, SEEDS, TWOM24, TWOM12, LUXLEV
      SAVE NSKIP, NDSKIP, IN24, NEXT, KOUNT, MKOUNT, INSEED
      INTEGER LUXLEV
      LOGICAL NOTYET
      DATA NOTYET, LUXLEV, IN24, KOUNT, MKOUNT /.TRUE., LXDFLT, 0,0,0/
      DATA I24,J24,CARRY/24,10,0./
C                               default
C  Luxury Level   0     1     2   *3*    4
      DATA NDSKIP/0,   24,   73,  199,  365 /
Corresponds to p=24    48    97   223   389
C     time factor 1     2     3     6    10   on slow workstation
C                 1    1.5    2     3     5   on fast mainframe
C
C  NOTYET is .TRUE. if no initialization has been performed yet.
C              Default Initialization by Multiplicative Congruential
      IF (NOTYET) THEN
         NOTYET = .FALSE.
         JSEED = JSDFLT  
         INSEED = JSEED
C         WRITE(8,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',JSEED
         LUXLEV = LXDFLT
         NSKIP = NDSKIP(LUXLEV)
         LP = NSKIP + 24
         IN24 = 0
         KOUNT = 0
         MKOUNT = 0
C         WRITE(8,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',
C     +        LUXLEV,'      p =',LP
            TWOM24 = 1.
         DO 25 I= 1, 24
            TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
   25    CONTINUE
         TWOM12 = TWOM24 * 4096.
         DO 50 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         NEXT(I) = I-1
   50    CONTINUE
         NEXT(1) = 24
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
      ENDIF
C
C          The Generator proper: "Subtract-with-borrow",
C          as proposed by Marsaglia and Zaman,
C          Florida State University, March, 1989
C
      DO 100 IVEC= 1, LENV
      UNI = SEEDS(J24) - SEEDS(I24) - CARRY 
      IF (UNI .LT. 0.)  THEN
         UNI = UNI + 1.0
         CARRY = TWOM24
      ELSE
         CARRY = 0.
      ENDIF
      SEEDS(I24) = UNI
      I24 = NEXT(I24)
      J24 = NEXT(J24)
      RVEC(IVEC) = UNI
C  small numbers (with less than 12 "significant" bits) are "padded".
      IF (UNI .LT. TWOM12)  THEN
         RVEC(IVEC) = RVEC(IVEC) + TWOM24*SEEDS(J24)
C        and zero is forbidden in case someone takes a logarithm
         IF (RVEC(IVEC) .EQ. 0.)  RVEC(IVEC) = TWOM24*TWOM24
      ENDIF
C        Skipping to luxury.  As proposed by Martin Luscher.
      IN24 = IN24 + 1
      IF (IN24 .EQ. 24)  THEN
         IN24 = 0
         KOUNT = KOUNT + NSKIP
         DO 90 ISK= 1, NSKIP
         UNI = SEEDS(J24) - SEEDS(I24) - CARRY
         IF (UNI .LT. 0.)  THEN
            UNI = UNI + 1.0
            CARRY = TWOM24
         ELSE
            CARRY = 0.
         ENDIF
         SEEDS(I24) = UNI
         I24 = NEXT(I24)
         J24 = NEXT(J24)
   90    CONTINUE
      ENDIF
  100 CONTINUE
      KOUNT = KOUNT + LENV
      IF (KOUNT .GE. IGIGA)  THEN
         MKOUNT = MKOUNT + 1
         KOUNT = KOUNT - IGIGA
      ENDIF
      RETURN
C
C           Entry to input and float integer seeds from previous run
      ENTRY RLUXIN(ISDEXT)
         TWOM24 = 1.
         DO 195 I= 1, 24
         NEXT(I) = I-1
         TWOM24 = TWOM24 * 0.5
  195    CONTINUE
         NEXT(1) = 24
         TWOM12 = TWOM24 * 4096.
C      WRITE(8,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
C      WRITE(8,'(5X,5I12)') ISDEXT
      DO 200 I= 1, 24
      SEEDS(I) = REAL(ISDEXT(I))*TWOM24
  200 CONTINUE
      CARRY = 0.
      IF (ISDEXT(25) .LT. 0)  CARRY = TWOM24
      ISD = IABS(ISDEXT(25))
      I24 = MOD(ISD,100)
      ISD = ISD/100
      J24 = MOD(ISD,100)
      ISD = ISD/100
      IN24 = MOD(ISD,100)
      ISD = ISD/100
      LUXLEV = ISD
        IF (LUXLEV .LE. MAXLEV) THEN
          NSKIP = NDSKIP(LUXLEV)
C          WRITE (8,'(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ',
C     +                         LUXLEV
        ELSE  IF (LUXLEV .GE. 24) THEN
          NSKIP = LUXLEV - 24
C          WRITE (8,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:',LUXLEV
        ELSE
          NSKIP = NDSKIP(MAXLEV)
C          WRITE (8,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ',LUXLEV
          LUXLEV = MAXLEV
        ENDIF
      INSEED = -1
      RETURN
C
C                    Entry to ouput seeds as integers
      ENTRY RLUXUT(ISDEXT)
      DO 300 I= 1, 24
         ISDEXT(I) = INT(SEEDS(I)*TWOP12*TWOP12)
  300 CONTINUE
      ISDEXT(25) = I24 + 100*J24 + 10000*IN24 + 1000000*LUXLEV
      IF (CARRY .GT. 0.)  ISDEXT(25) = -ISDEXT(25)
      RETURN
C
C                    Entry to output the "convenient" restart point
      ENTRY RLUXAT(LOUT,INOUT,K1,K2)
      LOUT = LUXLEV
      INOUT = INSEED
      K1 = KOUNT
      K2 = MKOUNT
      RETURN
C
C                    Entry to initialize from one or three integers
      ENTRY RLUXGO(LUX,INS,K1,K2)
         IF (LUX .LT. 0) THEN
            LUXLEV = LXDFLT
         ELSE IF (LUX .LE. MAXLEV) THEN
            LUXLEV = LUX
         ELSE IF (LUX .LT. 24 .OR. LUX .GT. 2000) THEN
            LUXLEV = MAXLEV
C            WRITE (8,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ',LUX
         ELSE
            LUXLEV = LUX
            DO 310 ILX= 0, MAXLEV
              IF (LUX .EQ. NDSKIP(ILX)+24)  LUXLEV = ILX
  310       CONTINUE
         ENDIF
      IF (LUXLEV .LE. MAXLEV)  THEN
         NSKIP = NDSKIP(LUXLEV)
C         WRITE(8,'(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :',
C     +        LUXLEV,'     P=', NSKIP+24
      ELSE
          NSKIP = LUXLEV - 24
C          WRITE (8,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:',LUXLEV
      ENDIF
      IN24 = 0
      IF(INS.LT.0)INS=-INS
C      IF (INS .LT. 0)  WRITE (6,'(A)')
C     +   ' Illegal initialization by RLUXGO, negative input seed'
      IF (INS .GT. 0)  THEN
        JSEED = INS
C        WRITE(8,'(A,3I12)') ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',
C     +      JSEED, K1,K2
      ELSE
        JSEED = JSDFLT
C        WRITE(8,'(A)')' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'
      ENDIF
      INSEED = JSEED
      NOTYET = .FALSE.
      TWOM24 = 1.
         DO 325 I= 1, 24
           TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
  325    CONTINUE
      TWOM12 = TWOM24 * 4096.
         DO 350 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         NEXT(I) = I-1
  350    CONTINUE
      NEXT(1) = 24
      I24 = 24
      J24 = 10
      CARRY = 0.
      IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
C        If restarting at a break point, skip K1 + IGIGA*K2
C        Note that this is the number of numbers delivered to
C        the user PLUS the number skipped (if luxury .GT. 0).
      KOUNT = K1
      MKOUNT = K2
      IF (K1+K2 .NE. 0)  THEN
        DO 500 IOUTER= 1, K2+1
          INNER = IGIGA
          IF (IOUTER .EQ. K2+1)  INNER = K1
          DO 450 ISK= 1, INNER
            UNI = SEEDS(J24) - SEEDS(I24) - CARRY
            IF (UNI .LT. 0.)  THEN
               UNI = UNI + 1.0
               CARRY = TWOM24
            ELSE
               CARRY = 0.
            ENDIF
            SEEDS(I24) = UNI
            I24 = NEXT(I24)
            J24 = NEXT(J24)
  450     CONTINUE
  500   CONTINUE
C         Get the right value of IN24 by direct calculation
        IN24 = MOD(KOUNT, NSKIP+24)
        IF (MKOUNT .GT. 0)  THEN
           IZIP = MOD(IGIGA, NSKIP+24)
           IZIP2 = MKOUNT*IZIP + IN24
           IN24 = MOD(IZIP2, NSKIP+24)
        ENDIF
C       Now IN24 had better be between zero and 23 inclusive
        IF (IN24 .GT. 23) THEN
C           WRITE (8,'(A/A,3I11,A,I5)')  
C     +    '  Error in RESTARTING with RLUXGO:','  The values', INS,
C     +     K1, K2, ' cannot occur at luxury level', LUXLEV
C           IN24 = 0
        ENDIF
      ENDIF
      RETURN
      END
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     ***************************************************************
c
      subroutine ncoord(unit,dchar,asymp)
c
c     Reads in a new set of coordinates
c
      use mpi
      implicit double precision (a-h,m-z)
      character*30 unit,dchar,dummy
      parameter (ixlen=1000)
      dimension imass(ixlen),xmass(ixlen) 
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xeo,xk,xn,xmv,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(ixlen),k_const(ixlen),tcharge
      common/coordinates/fx(ixlen),fy(ixlen),fz(ixlen),
     ?ox(ixlen),oy(ixlen),oz(ixlen)
      common/ljparameters/eolj(ixlen),rolj(ixlen),eox4(ixlen),
     ?ro6lj(ixlen),ro12lj(ixlen),dro6(ixlen),dro12(ixlen)
      common/ljparameters2/alphai(ixlen),Ni(ixlen),Ai(ixlen),Gi(ixlen),
     1 RijStar(ixlen),eij(ixlen)
      common/mmffN2/alphaN2,NiN2,AiN2,GiN2,RN2Star,MMFF_B,MMFF_beta
      common/mmffHe/alphaHe,NiHe,AiHe,GiHe,RHeStar
      common/hsparameters/rhs(ixlen),rhs2(ixlen)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
      common/scaling/enersca,distsca
      common/mpicon/imyrank,inprocs,imp_per_node,inp_per_node,s_time
      common/runparams/itn,inp,imp,igas
c
      read(9,'(a30)',end=100) dummy
100   continue
c
      do 2000 iatom=1,inatom
       if(dchar.eq.'calc') then
        read(9,*) fx(iatom),fy(iatom),fz(iatom),
     1  ximass,pcharge(iatom)
       else
        read(9,*) fx(iatom),fy(iatom),fz(iatom),ximass
       endif
       imass(iatom)=nint(ximass)
       if(unit.eq.'au') then
        fx(iatom)=fx(iatom)*0.52917706d0
        fy(iatom)=fy(iatom)*0.52917706d0
        fz(iatom)=fz(iatom)*0.52917706d0
       endif
2000  continue
c
      mx=0.d0
      do 2021 iatom=1,inatom
       mx=mx+xmass(iatom)
2021  continue
c
      if(mx.ne.m2) then
       write(8,624)
       if(imyrank.eq.0)close(8)
       call MPI_FINALIZE(ierr)
       stop
      endif
c
      do 2030 iatom=1,inatom
       rhs2(iatom)=rhs(iatom)*rhs(iatom)
       eox4(iatom)=4.d0*eolj(iatom)
       ro2lj=rolj(iatom)*rolj(iatom)
       ro6lj(iatom)=ro2lj*ro2lj*ro2lj
       ro12lj(iatom)=ro6lj(iatom)*ro6lj(iatom)
       dro6(iatom)=6.d0*ro6lj(iatom)
       dro12(iatom)=12.d0*ro12lj(iatom)
2030  continue
c
      fxo=0.d0
      fyo=0.d0
      fzo=0.d0
      do 2009 iatom=1,inatom
       fxo=fxo+(fx(iatom)*xmass(iatom))
       fyo=fyo+(fy(iatom)*xmass(iatom))
       fzo=fzo+(fz(iatom)*xmass(iatom))
2009  continue
      fxo=fxo/m2
      fyo=fyo/m2
      fzo=fzo/m2
      do 2010 iatom=1,inatom
       fx(iatom)=(fx(iatom)-fxo)*1.d-10
       fy(iatom)=(fy(iatom)-fyo)*1.d-10
       fz(iatom)=(fz(iatom)-fzo)*1.d-10
2010  continue
c
      do 3000 iatom=1,inatom
       ox(iatom)=fx(iatom)
       oy(iatom)=fy(iatom)
       oz(iatom)=fz(iatom)
3000  continue
c
c     determine structural asymmetry parameter
c
      theta=0.d0
      asymp=0.d0
      do 5050 igamma=0,360,2
       do 5000 iphi=0,180,2
        gamma=dble(igamma)/cang
        phi=dble(iphi)/cang
        call rotate
        xyzsum=0.d0
        yzsum=0.d0
        do 5005 iatom=1,inatom
         xyz=dsqrt(fx(iatom)**2+fy(iatom)**2+fz(iatom)**2)
         yz=dsqrt(fy(iatom)**2+fz(iatom)**2)
         xyzsum=xyzsum+xyz
         yzsum=yzsum+yz
5005    continue
        hold=((pi/4.d0)*xyzsum)/yzsum
        if(hold.gt.asymp) asymp=hold
5000   continue
5050  continue
c
624   format(1x,'masses do not add up')
c
      return
      end
c
c     ***************************************************************
c
