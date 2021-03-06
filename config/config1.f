       program config
       implicit none
       integer type,A,cent2,naccept,d,rec_sw,recover_sw,recover_rec
       integer atom_num,ens_num,par_num,samebc,MAXGEN,RESETSEED,asw
       integer i,min,num,ic,z,MAXIT,MUPDATE,TAKESTEP,ndim_test2,FIX_GFAC
       integer ndim_test,residues_num(10000),j1,res_dim,load3,GBIN,DG_SW
       integer ch_st(10000),ch_end(10000),pdb_rec,dim3,tmp,gauss
       integer ibin,N0,FIXTEMP_T,FIXSTEP_T,FIXBOTH_T,GMAX_SW,DG_RECOVER
       integer STEPMAX1,MCMAX1,STEPMAX2,MCMAX2,STEPMAX3,MCMAX3,sw_eular1
       integer BH_SW1,BH_SW2,BH_SW3,ERANGE_SW,TEMP_EX,sw_eular2
       integer sw_ensnum1,sw_ensnum2,sw_ensnum3,sw_parnum1,sw_parnum2,
     &sw_parnum3,sw_same1,sw_same2,sw_same3
       real*8 aa,bb,cc,dd,ee,ff,gg,hh,ii,LL,ene,enef,enef2,dx,bc,simp,cg
       real*8 epsilon,zeta,p,q,rzero,epsilon0,zeta0,p0,q0,rzero0
       real*8 epsilon2,zeta2,p2,q2,rzero2,epsilon3,zeta3,p3,q3,rzero3
       real*8 DGUESS,ASTEP,STEP,GTOL,accrat,rec_ratio,GMAX
       real*8 recover_ratio,CUTOFF,THETA,CUTOFF_ANGLE,EPS,EPS2,EPS3
       real*8 EPS4,ACOIL,BCOIL,CCOIL,DCOIL,ECOIL,FCOIL,GCOIL,EXC_D
       real*8 THETA1,dih_step,a1,b1,c1,d1,e1,f1,g1,h1
       real*8 moment(20),frange,arange,eular_fac1,eular_fac2
       real*8 BETA0,TEMPFD,EMAX,EMIN,SCALEFAC,TEMPI,TEMPF,GFAC
       real*8 mva1,mvb1,mvc1,mvd1,mve1,mvf1,mvg1,mvh1,mvi1,mvj1,rms_fac1
       real*8 mva2,mvb2,mvc2,mvd2,mve2,mvf2,mvg2,mvh2,mvi2,mvj2,rms_fac2
       real*8 mva3,mvb3,mvc3,mvd3,mve3,mvf3,mvg3,mvh3,mvi3,mvj3,rms_fac3
       real*8 rms_start,rms_end,rms_rate,DEMIN,DEMAX
       real*8 MOVE_FAC1,MOVE_FAC2,STEP_MIN,STEP_MAX,ASTEP_MIN,ASTEP_MAX
       real*8 TEMP_MIN,TEMP_MAX
       real*8 eam_epsilon,eam_epsilon1,eam_epsilon2,eam_epsilon3
       real*8 eam_a,eam_a1,eam_a2,eam_a3,eam_C,eam_C1,eam_C2,eam_C3
       real*8 eam_n,eam_n1,eam_n2,eam_n3,eam_m,eam_m1,eam_m2,eam_m3
       character load*4,kind0*3,kind*2,kind2*2,kind3*2
       character load2*2,kind4*1,kind5*1,prot(10000)*3
       d=1
c------------------------ Application ----------------------------
       write(*,*)
       write(*,*)
       write(*,*)
       write(*,*) "                                             "
       write(*,*) "   MPI-PTMBHPGA For VARIOUS POTENTAIL        "
       write(*,*) "          by K.L Wu & P.J Hsu                "
       write(*,*) "                                             "
       write(*,*) "             VERSION  v2.4                   "
       write(*,*) "                                             "
       write(*,*)
       write(*,*)
c-------------------- Background Steps ----------------------------          
       open(1,file="atom_num.dat",status="old")
       read(1,*) atom_num,num
       close(1)
       write(*,*) "Please input the particle number:"
       read(*,*) atom_num
       enef=1.D0+(3.0D0*atom_num/17.77153175D0)**(1.0D0/3.0D0)
       enef=enef*2.0D0**(1.0D0/6.0D0)
       enef=enef*enef
       if(num.ne.0) num=atom_num-1
1      write(*,*) "Previous data is needed when you are using seeding"
       write(*,*) "(If you type 'y',the Seeding number is",num,")"
       write(*,*) "Load the previous configuration ? (y/n)"
       read(*,*) load
c----------------------- Start to load previous config. ----------
       open(10,file="data.dat")
       read(10,*) LL,ene,enef,enef2,simp,cg,gauss,d
       read(10,*) ens_num,min,aa,bb,cc,dd,ee,ff,gg,hh,ii,par_num,samebc
     &,MAXGEN
       read(10,*) dx,bc,asw,sw_eular1,sw_eular2,eular_fac1,eular_fac2
       read(10,*) MAXIT,MUPDATE,TAKESTEP,RESETSEED,DGUESS,ASTEP,STEP,
     &GTOL,GMAX,GMAX_SW
       read(10,*) rec_sw,rec_ratio,pdb_rec,dim3
       read(10,*) recover_sw,recover_rec,recover_ratio
       read(10,*) 
     &ibin,BETA0,N0,TEMPFD,EMAX,EMIN,naccept,accrat
       read(10,*) FIXTEMP_T,FIXSTEP_T,FIXBOTH_T,SCALEFAC,TEMPI,TEMPF
     &,TEMP_EX
       read(10,*) mva1,mvb1,mvc1,mvd1,mve1,mvf1,mvg1,mvh1,mvi1,mvj1
       read(10,*) mva2,mvb2,mvc2,mvd2,mve2,mvf2,mvg2,mvh2,mvi2,mvj2
       read(10,*) mva3,mvb3,mvc3,mvd3,mve3,mvf3,mvg3,mvh3,mvi3,mvj3
       read(10,*) GBIN,GFAC,FIX_GFAC,DG_SW,DG_RECOVER
       read(10,*) STEPMAX1,MCMAX1,STEPMAX2,MCMAX2,STEPMAX3,MCMAX3
       read(10,*) BH_SW1,BH_SW2,BH_SW3,DEMIN,DEMAX,ERANGE_SW
       read(10,*) rms_start,rms_end,rms_rate,rms_fac1,rms_fac2,rms_fac3
       read(10,*) sw_ensnum1,sw_parnum1,sw_same1
       read(10,*) sw_ensnum2,sw_parnum2,sw_same2
       read(10,*) sw_ensnum3,sw_parnum3,sw_same3
       read(10,*) MOVE_FAC1,MOVE_FAC2,TEMP_MIN,TEMP_MAX
       read(10,*) STEP_MIN,STEP_MAX,ASTEP_MIN,ASTEP_MAX
       close(10)
C   INITIALIZE THE RECOVER FUNCTION
       recover_sw=0
       open(25,file="type.dat",status="old")
       read(25,*) type,cent2
       close(25)
       open(26,file="alloy.dat",status="old")
       read(26,*) A
       close(26)
       open(30,file="parameter1.dat",status="old")
       read(30,*) kind,epsilon,zeta,p,q,rzero
       close(30)
       open(31,file="parameter2.dat",status="old")
       read(31,*) kind2,epsilon2,zeta2,p2,q2,rzero2
       close(31)
       open(32,file="parameter3.dat",status="old")
       read(32,*) kind3,epsilon3,zeta3,p3,q3,rzero3
       close(32)
       open(33,file="protein1.dat",status="old")
       read(33,*)CUTOFF,THETA,CUTOFF_ANGLE,EPS,EPS2,EPS3,EPS4,ACOIL,
     &BCOIL,CCOIL,DCOIL,ECOIL
       close(33)
       open(34,file="protein2.dat",status="old")
       ndim_test=0
11     read(34,*,end=22)
       ndim_test=ndim_test+1
       goto 11
22     rewind(34)
       do j1=1,ndim_test
         read(34,*) prot(j1),residues_num(j1)
       enddo
       close(34)
       if(type.eq.5)then
           open(50,file="eam_aa.dat")
           read(50,*) eam_a1,eam_epsilon1,eam_C1,eam_n1,eam_m1
           close(50)
       endif
       open(37,file="chima1.dat",status="old")
       read(37,*) FCOIL,GCOIL,THETA1,EXC_D,dih_step
       read(37,*) a1,b1,c1,d1,e1,f1,g1,h1
       open(38,file="chima2.dat",status="old")
       ndim_test2=0
21     read(38,*,end=23)
       ndim_test2=ndim_test2+1
       goto 21
23     rewind(38)
       do j1=1,ndim_test2
         read(38,*) ch_st(j1),ch_end(j1)
       enddo
       close(38)
       close(37)
       if (load.eq."y")then
         open(18,file="atom_num.dat",status="old")
         write(18,*) atom_num,num
         close(18)
         if (min.eq.5)then
           enef=1.D0+(3.0D0*atom_num/17.77153175D0)**(1.0D0/3.0D0)
           enef=enef*2.0D0**(1.0D0/6.0D0)
           enef=enef*enef
           if (d.eq.1) LL=enef
         endif
         open(12,file="data.dat",status="replace")
         write(12,*) LL,ene,enef,enef2,simp,cg,gauss,d
         write(12,*) ens_num,min,aa,bb,cc,dd,ee,ff,gg,hh,ii,par_num,
     &samebc,MAXGEN
         write(12,*) dx,bc,asw,sw_eular1,sw_eular2,eular_fac1,eular_fac2
         write(12,*) MAXIT,MUPDATE,TAKESTEP,RESETSEED,DGUESS,ASTEP,STEP,
     &GTOL,GMAX,GMAX_SW
         write(12,*) rec_sw,rec_ratio,pdb_rec,dim3
         write(12,*) recover_sw,recover_rec,recover_ratio
         write(12,*)
     &ibin,BETA0,N0,TEMPFD,EMAX,EMIN,naccept,accrat
         write(12,*) FIXTEMP_T,FIXSTEP_T,FIXBOTH_T,SCALEFAC,TEMPI,TEMPF
     &,TEMP_EX
         write(12,*) mva1,mvb1,mvc1,mvd1,mve1,mvf1,mvg1,mvh1,mvi1,mvj1
         write(12,*) mva2,mvb2,mvc2,mvd2,mve2,mvf2,mvg2,mvh2,mvi2,mvj2
         write(12,*) mva3,mvb3,mvc3,mvd3,mve3,mvf3,mvg3,mvh3,mvi3,mvj3
         write(12,*) GBIN,GFAC,FIX_GFAC,DG_SW,DG_RECOVER
         write(12,*) STEPMAX1,MCMAX1,STEPMAX2,MCMAX2,STEPMAX3,MCMAX3
         write(12,*) BH_SW1,BH_SW2,BH_SW3,DEMIN,DEMAX,ERANGE_SW
       write(12,*) rms_start,rms_end,rms_rate,rms_fac1,rms_fac2,rms_fac3
         write(12,*) sw_ensnum1,sw_parnum1,sw_same1
         write(12,*) sw_ensnum2,sw_parnum2,sw_same2
         write(12,*) sw_ensnum3,sw_parnum3,sw_same3
         write(12,*) MOVE_FAC1,MOVE_FAC2,TEMP_MIN,TEMP_MAX
         write(12,*) STEP_MIN,STEP_MAX,ASTEP_MIN,ASTEP_MAX
         close(12)

         write(*,*) "@ ",atom_num," particles"
         write(*,*) "@ seeding number= ",num,",(0 for no seeding)"
         if (gauss.eq.1) then
           write(*,*) "@ using Gaussian Random number"
         endif
         if(min.eq.1)then
           write(*,*) "@ using Simplex minimization"
         else if(min.eq.2)then
           write(*,*) "@ using Conjugate Gradient minimization"
         else if(min.eq.3)then
           write(*,*) "@ using Conjugate Gradient then Simplex"
         else if(min.eq.4)then
           write(*,*) "@ using Simplex then Conjugate Gradient"   
         else
           write(*,*) "@ using LBFGS minimization,MAXIT=",MAXIT
           if(TAKESTEP.eq.1)then
             write(*,*) "@ Take step!!"
           else
             write(*,*) "@ Don't take step!"
           endif
           if(RESETSEED.eq.1)then
             write(*,*) "@ Reset to seeding number!!"
           else
             write(*,*) "@ Don't reset to seeding number!!"
           endif
         endif
         if(cent2.eq.1)then
           write(*,*) "@ move coords to center of mass!"
         else
           write(*,*) "@ Do not move coords to center of mass!"
         endif
         write(*,*) "@ ",ens_num," ensembles for global GA"
         write(*,*) "@ ",sw_ensnum1," ensembles for parallel GA run 1"
         write(*,*) "@ ",sw_ensnum2," ensembles for parallel GA run 2"
         write(*,*) "@ ",sw_ensnum3," ensembles for parallel GA run 3"
         if (sw_eular1.eq.1)then 
           write(*,*) "@ Use eularangle tran for GGA,factor=",eular_fac1
         else
           write(*,*) "@ Do not use eularangle tran for GGA"
         endif
         if (sw_eular2.eq.1)then
           write(*,*) "@ Use eularangle tran for PGA,factor=",eular_fac2
         else
           write(*,*) "@ Do not use eularangle tran for PGA"
         endif
         write(*,*) "@ The ratio of GA:"
         write(*,*) "@ Inversion:    ",aa 
         write(*,*) "@ Arithmetic:   ",bb
         write(*,*) "@ Geometic:     ",cc
         write(*,*) "@ Crossing:     ",dd
         write(*,*) "@ 2-Point:      ",ee
         write(*,*) "@ 3N Mutation : ",ff
         write(*,*) "@ Moment op 1:  ",gg
         write(*,*) "@ Moment op 2:  ",hh
         write(*,*) "@ DIH move op 1:",ii
         write(*,*) "@ The ratio of Step move:"
         write(*,*) "@ 3N Particles Moving,  1st=",mva1,",2nd=",mva2,
     &",3rd=",mva3
         write(*,*) "@ 3N Particles Mutation,1st=",mvb1,",2nd=",mvb2,
     &",3rd=",mvb3
         write(*,*) "@ N Dimension Moving,   1st=",mvc1,",2nd=",mvc2,
     &",3rd=",mvc3
         write(*,*) "@ 3N Dihedral Moving,   1st=",mve1,",2nd=",mve2,
     &",3rd=",mve3
         write(*,*) "@ RMS Transition,       1st=",mvf1,",2nd=",mvf2,
     &",3rd=",mvf3
         write(*,*) "@ RMS SCALE,            1st=",rms_fac1,",2nd=",
     &rms_fac2,",3rd=",rms_fac3
         write(*,*) "@ Inversion,            1st=",mvd1,",2nd=",mvd2,
     &",3rd=",mvd3
         write(*,*) "@ Arithmetic,           1st=",mvg1,",2nd=",mvg2,
     &",3rd=",mvg3
         write(*,*) "@ Geometic,             1st=",mvh1,",2nd=",mvh2,
     &",3rd=",mvh3
         write(*,*) "@ Crossing,             1st=",mvi1,",2nd=",mvi2,
     &",3rd=",mvi3
         write(*,*) "@ 2-Point Crossover,    1st=",mvj1,",2nd=",mvj2,
     &",3rd=",mvj3
         write(*,*) "@",par_num," parents for global GA"
         write(*,*) "@",sw_parnum1," parents for parallel GA run 1"
         write(*,*) "@",sw_parnum2," parents for parallel GA run 2"
         write(*,*) "@",sw_parnum3," parents for parallel GA run 3"
         write(*,*) "@ same=",samebc," for global GA"
         write(*,*) "@ same=",sw_same1," for parallel GA run 1"
         write(*,*) "@ same=",sw_same2," for parallel GA run 2"
         write(*,*) "@ same=",sw_same3," for parallel GA run 3"
         write(*,*) "@ Max gens=",MAXGEN
         write(*,*) "@ Multi canonical max number for 1st Run=",MCMAX1
        write(*,*) "@ Max MC steps for each sweep for 1st Run=",STEPMAX1
         write(*,*) "@ Multi canonical max number for 2nd Run=",MCMAX2
        write(*,*) "@ Max MC steps for each sweep for 2nd Run=",STEPMAX2
         write(*,*) "@ Multi canonical max number for 3rd Run=",MCMAX3
        write(*,*) "@ Max MC steps for each sweep for 3rd Run=",STEPMAX3
         if (d.eq.1) then
           write(*,*) "@ Using close pack condition"
         endif
         write(*,*) "@ Width of initial box is ",LL
         write(*,*) "@ BC in func=",enef
         write(*,*) "@ Accuracy parameter in func=",enef2
         write(*,*) "@ The energy sort accuracy is ",ene
         write(*,*) "@ The dX in dfunc is ",dx
         write(*,*) "@ BC in dfunc is ",bc
         write(*,*) "@ Tolerance of Simplex is ",simp
         write(*,*) "@ Tolerance of Conjugate Gradient is ",cg  
         write(*,*) "@ Tolerance of LBFGS(GTOL)=",GTOL
         if(GMAX_SW.eq.1)then
      write(*,*) "@ Using RMS force for sloppy convergence criterion"
           write(*,*) "@ Threshold criterion=",GMAX
         else if(GMAX_SW.eq.2)then
      write(*,*) "@ The sloppy convergence criterion for RMS force=",
     &GMAX
         else if(GMAX_SW.eq.3)then
      write(*,*) "@ Use multi-convergence criterion from ",rms_start,
     &" to ",rms_end
           write(*,*) "@ The variation rate is ",rms_rate
         endif
         write(*,*) "@ MCSCRH UPDATE =",MUPDATE
         write(*,*) "@ Takestep in bmin:Astep=",ASTEP,",Step=",STEP
         write(*,*) "@ Max Step=",STEP_MAX,",Min Step=",STEP_MIN
         write(*,*) "@ Max Astep=",ASTEP_MAX,",Min Astep=",ASTEP_MIN
         write(*,*) "@ Max TEMP=",TEMP_MAX,",Min TEMP=",TEMP_MIN
         write(*,*) "@ The adjustment of the Acceptance checking is +-",
     &MOVE_FAC1
         write(*,*) "@ Input DGUESS=",DGUESS
         if(DG_SW.eq.1)then
           write(*,*) "@ Do not use Guided function"
         else
           write(*,*) "@ Use Guided function"
           write(*,*) "@ # of Guided Function BINS=",GBIN
         endif
         if (FIX_GFAC.eq.1)then
           write(*,*) "@ Using temperature for Guided function=",GFAC
         else
           write(*,*) "@ The acceptance probability for Guided fun is"
     &,GFAC
         endif
         write(*,*) "@ # of Multicanonical BH BINS=",ibin
       if(BH_SW1.EQ.1.or.BH_SW1.eq.4.or.BH_SW2.eq.1.or.BH_SW3.eq.1)then
           write(*,*) "@ INITIAL BETA0=1/(",TEMPF,"-",TEMPFD,")=",BETA0
           if(ERANGE_SW.eq.0)then
             write(*,*) "@ ENERGY RANGE FROM ",EMAX," TO ",EMIN
           else if(ERANGE_SW.eq.1)then
             write(*,*) "@ Using Auto Energy Range(Half)"
             write(*,*) "@ TDE from ",DEMIN," to ",DEMAX
           else if(ERANGE_SW.eq.2)then
             write(*,*) "@ Using Multi-Canonical Annealing!"
             write(*,*) "@ TDE from ",DEMIN," to ",DEMAX
           endif
         else
           write(*,*) "@ Using MonteCarlo Boltzmann distribution"
         endif
      write(*,*)"@ Every ",naccept," generations check acceptance ratio"
      write(*,*)"@ Acceptance ratio=",accrat
      write(*,*) "@ The adjustment of the Acceptance Ratio is +-",
     &MOVE_FAC2
      write(*,*) "@ Temperature start from ",TEMPI," to ",TEMPF,
     & ",multiply by ",SCALEFAC
         if(TEMP_EX.eq.1)then
           write(*,*) "@ Use temperature range for all ensembles!"
         else if(TEMP_EX.eq.2)then
           write(*,*) "@ Use temperature range for the children!"
         else
           write(*,*) "ERROR!!"
           return
         endif
         if(FIXTEMP_T.eq.1.or.FIXBOTH_T.eq.1)then
           write(*,*) "@ FIX TEMPERATURE!"
         else
           write(*,*) "@ DO NOT FIX TEMPERATURE!"
         endif
         if(FIXSTEP_T.eq.1.or.FIXBOTH_T.eq.1)then
           write(*,*) "@ FIX STEPS!"
         else
           write(*,*) "@ DO NOT FIX STEPS!"
         endif
         if(rec_sw.eq.0)then
           write(*,*) "@ Normal recording"
         endif
         if(rec_sw.eq.1)then
           if (pdb_rec.eq.1)then
             write(*,*) "@ Step recording in PDB file!"
           else
             write(*,*) "@ Step recording in seperate coords file!"
           endif
           write(*,*) "@ Step recording,recording ratio=",rec_ratio
         endif
         if(rec_sw.eq.2)then
           write(*,*) "@ Record nothing!"
         endif
         if(recover_rec.eq.1)then
           write(*,*) "@ Start recover,recoverratio=",recover_ratio
         endif
         if (type.eq.1)then
           write(*,*) "@    e0       c0        p       q      r0"
           write(*,*) "@ ",epsilon,",",zeta,",",p,",",q,",",rzero
           write(*,*) "@ For ",kind
         endif
         if (type.eq.2)then
           write(*,*) "@    e0       c0        p       q      r0"
           write(*,*) "@ ",epsilon,",",zeta,",",p,",",q,",",rzero
           write(*,*) "@ For ",A," ",kind
           write(*,*) "@    e0       c0        p       q      r0"
      write(*,*) "@ ",epsilon2,",",zeta2,",",p2,",",q2,",",rzero2
           write(*,*) "@ For ",atom_num-A," ",kind2     
           write(*,*) "@    e0       c0        p       q      r0"
      write(*,*) "@ ",epsilon3,",",zeta3,",",p3,",",q3,",",rzero3
           write(*,*) "@ For ",kind3       
         endif
         if (type.eq.3)then
           write(*,*) "@ For Beta-Sheets"
           write(*,*) "@ ",(prot(j1),residues_num(j1),j1=1,ndim_test)
         endif
         if (type.eq.4)then
           write(*,*) "@ For Alpha-Helices"
           write(*,*) "@ ",(prot(j1),residues_num(j1),j1=1,ndim_test)
         endif
c=======For EAM Potential=====================================================================
         if (type.eq.5)then
c           write(*,*) "@  eam_a  eam_e  eam_C  eam_n  eam_m"
c           write(*,*) "@ ",eam_a1,",",eam_epsilon,",",eam_C1,","
c     &,eam_n1,",",eam_m1
c           write(*,*) "@ For ",kind
           write(*,*) "EAM Potential Parameter"
           write(*,*) eam_a1,eam_epsilon1,eam_C1,eam_n1,eam_m1
           write(*,*) "======================"
         endif
c================================================================================
         if (type.eq.7)then
           write(*,*) "@ For Chen and Imamura Folding N"
           write(*,*) "@ The Neutral bond sets:(Start-end)"
        write(*,*) "@ ",("(",ch_st(j1),ch_end(j1),")",j1=1,ndim_test2)
           write(*,*) "@ Lo=",FCOIL,",Epsilon=",GCOIL,",ANGLE=",THETA1
     &,",Exclusive diameter=",EXC_D
         endif
         if (type.eq.8)then
           write(*,*) "@ For Chen and ImamuraFolding 3N"
           write(*,*) "@ The Neutral bond sets:(Start-end)"
        write(*,*) "@ ",("(",ch_st(j1),ch_end(j1),")",j1=1,ndim_test2)
           write(*,*) "@ Lo=",FCOIL,",Epsilon=",GCOIL,",ANGLE=",THETA1
     &,",Exclusive diameter=",EXC_D
         endif
         if (asw.eq.1)then
           write(*,*) "@ Using analytic gradient"
         else
           write(*,*) "@ Using numerical gradient"
         endif
         write(*,*)
c------------------------ end to load previous config. -----------
       else if(load.eq."n")then
2        write(*,*) "Do you want to: (choose a number)"
         write(*,*)  
         write(*,*) "0.Change recover method.."
         write(*,*) "1.Change Seeding number.. "
         write(*,*) "2.Change minimization method.."
         write(*,*) "3.Change record/backup method.."
         write(*,*) "4.Change Program accuracy.."
         write(*,*) "5.Change ensemble parameters.. "
         write(*,*) "6.Change the ratio of moving.."
         writE(*,*) "7.Change energy parameters (for Many body)"
         write(*,*) "8.Quit without save and see previous data"
         write(*,*) "9.Quit and save !"
         write(*,*) "10.On-line help !"
         read(*,*) ic
c---------------------- Online help -------------------------------
****IC=0
         if(ic.eq.0)then
3          write(*,*) "Do you want to recover your job?"
           write(*,*) "Please make sure you set the right material and"
           write(*,*) "particle number(VERY IMPORTANT)"
       write(*,*) "Note:This function can only recover the latest job"
       write(*,*) "You may check 8 to see the previous parameters first"
           write(*,*) "Attention!!"
           write(*,*) "Before you run the main program(GA),please"
          write(*,*) "backup the recover files,once you run the program"
          write(*,*) "The backup file will lost!"
           write(*,*) "Recover the latest job?(y/n)"
           read(*,*) load
           if(load.eq."y")then
             write(*,*) "You have turned on the recover mode"
             write(*,*) "Good Luck!!!" 
             recover_sw=1
           else if(load.eq."n")then
             recover_sw=0
           else
             call ERROR
             goto 3
           endif
4          write(*,*) "Start from the guided data?"
           write(*,*) "It require the previous guided data"
           write(*,*) "Start from the Guided data(y/n)(default=n)"
           read(*,*) load
           if(load.eq."y")then
             DG_RECOVER=1
           else if(load.eq."n")then
             DG_RECOVER=0
           else
             call ERROR
             goto 4
           endif
           goto 2
****IC=10         
         else if(ic.eq.10)then
           write(*,*) "              Welcome to online help !          "
           write(*,*)
           write(*,*) "If you are not sure how to change the accuracy"
           write(*,*) "Please input the default value"
           write(*,*)
     &"Program will produce several files: (# is particle number)"
           write(*,*) "#.dat    =>final result !"
           write(*,*) "#.pos    =>particle position !"
           write(*,*) "#.xyz    =>for xmakemol !"
           write(*,*) "#.rec    =>the latest result !"
           write(*,*)
     &"#_both.sed =>time data, conj. grad., simplex with seeding "
           write(*,*) "#_simp.sed =>time data,simplex with seeding !"
           write(*,*)
     &"#_conj.sed =>time data,conjugate gradient with seeding !"
           write(*,*)
     &"#_both.sed.ini =>time data, initial time recording !"
           write(*,*)
     &"#_both.sed.genper =>time data,mean time per generation !"
           write(*,*)
     &"#_both.nos =>like '#_both.sed' but without seeding !"
           write(*,*)
           write(*,*) 
     &"Tip1: You can change parents number and same number"
           write(*,*) "Tip2: Please refer 'enrgpmt.dat'"
           pause
           write(*,*) "              END OF HELP                     "
           goto 2   
***IC=1                                                      
         else if(ic.eq.1)then
           write(*,*) "Input Seeding number (0 for no seeding)"
           num=atom_num-1
           write(*,*) "Input the seeding number: (default=",num,")"
           read(*,*) num
           goto 2
***IC=2
         else if(ic.eq.2)then
24         write(*,*) "Minimization method to use ? (default=5)"
           write(*,*) "1. Only Simplex"
           write(*,*) "2. Only Conjugate Gradient"
           write(*,*) "3. Conjugate Gradient then Simplex"
           write(*,*) "4. Simplex then Conjugate Gradient"    
           write(*,*) "5. LBFGS"
           write(*,*) "6. Do not use any local minimum method!"
           read(*,*) min
        if(min.eq.1.or.min.eq.2.or.min.eq.3.or.min.eq.4.or.min.eq.5)then
           else
             call ERROR
             goto 24
           endif
           if(min.eq.2.or.min.eq.3.or.min.eq.4)then
             enef=500.
           endif
25         write(*,*) "1.move coords to center of mass!"
           write(*,*) "2.Do not move coords to center of mass!"
           write(*,*) "(default=2)"
           read(*,*) cent2
           if(cent2.eq.1.or.cent2.eq.2)then
           else
             call ERROR
             goto 25
           endif
           goto 2
***IC=3  
         else if(ic.eq.3)then
26         write(*,*) "Select a recording mode:(default=0)"
           write(*,*) "0.Normal recording(For lowest rec purpose)"
           write(*,*) "1.Step recording(For step evolution)"
           write(*,*) "2.Record nothing(faster)!"
           read(*,*) rec_sw
           if (rec_sw.eq.0.or.rec_sw.eq.1.or.rec_sw.eq.2)then
           else
             call ERROR
             goto 26
           endif
           if(rec_sw.ne.2.and.rec_sw.ne.1.and.rec_sw.ne.0) then
             write(*,*) "You have choose a invalid number,reset to 0"
             rec_sw=0
           endif
           if(rec_sw.eq.1)then
29           write(*,*) "1.Step record in PDB Animation file"
             write(*,*) "2.Step record in seperate files"
             read(*,*) pdb_rec
             if (pdb_rec.eq.1.or.pdb_rec.eq.2)then
             else
               call ERROR
               goto 29
             endif
             write(*,*) "Input the recording ratio:(default=1.D-6)"
             read(*,*) rec_ratio
           endif
32      write(*,*) "Do you want to backup the parents(y/n)?(default=y)"
         write(*,*) "Note:This may increase the computing time!"
           read(*,*) load
           if(load.eq."y")then
             recover_rec=1
             write(*,*) "Please input the recover ratio(default=1.D-4)"
             read(*,*) recover_ratio
           else if(load.eq."n")then
             recover_rec=0
             recover_ratio=1.D-4
           else
             call ERROR
             goto 32
           endif
           goto 2
***IC=4
         else if(ic.eq.4)then
58         write(*,*) "If you are not sure for some parameter.."
           write(*,*) "Please type in the default value !"
           write(*,*) "1.Change local minima serching parameters"
           write(*,*) "2.Change the Monte Carlo depent parameters"
           write(*,*) "3.Change the Guided Function parameters"
           write(*,*) "4.Back to main menu"
           read(*,*) tmp
           if (tmp.eq.1) then
             write(*,*) "Input energy sort accuracy: (default=1.D-12)"
             read(*,*) ene
             if(min.eq.2.or.min.eq.3.or.min.eq.4)then
               write(*,*) "Input func accuracy: (default=500.D0)"
               read(*,*) enef
               write(*,*) "Input func fit parameter: (default=1.D0)"
               read(*,*) enef2
             endif
             if(min.ne.1)then
41             write(*,*) "Do you want to:"
               write(*,*) "1.Use analytic gradient"
               write(*,*) "2.Use numerical gradient"
               write(*,*) "(If the analytic gradient is not available"
               write(*,*) "please choose 2)"
               read(*,*) asw
               if (asw.eq.1)then
               else if(asw.eq.2)then
                 write(*,*) "Input dX in grad: (default=7.D-5)"
                 read(*,*) dx
            write(*,*) "Input the reduced factor in grad:(default=1.D0)"
                 read(*,*) bc  
               else
                 call ERROR
                 goto 41
               endif
             endif
             if (min.ne.2.and.min.ne.5)then
               write(*,*) "Input tolerance of simplex (default=1.D-08)"
               write(*,*) "for simplex method" 
               read(*,*) simp
               write(*,*) "Input func accuracy: (default=500)"
               read(*,*) enef
             endif
             if (min.ne.1.and.min.ne.5)then
           write(*,*) "Input the tolerance of Conjugate(default=1.D-08)"
               write(*,*) "for conjugate gradient"
               read(*,*) cg
             endif
             if(min.eq.5)then
           write(*,*) "Input the tolerance of LBFGS(GTOL,default=0.9D0)"
               write(*,*) "This number should between 0.1 and 0.9"
               read(*,*) GTOL
               write(*,*) "Input the MCSRCH update(default=4)"
               read(*,*) MUPDATE
               write(*,*) "Input DGUESS in LBFGS(default=0.1D0)"
               read(*,*) DGUESS
42      write(*,*) "1.To use RMS force for sloppy convergence criterion"
        write(*,*) "2.Input sloppy convergence criterion"
        write(*,*) "3.Use variate convergence criterion"
        write(*,*) "(default=1)"
               read(*,*) GMAX_SW
               if(GMAX_SW.eq.1)then
      write(*,*)"Input threshold convergence criterion(default=0.001D0)"
                 read(*,*) GMAX
               else if(GMAX_SW.eq.2)then
        write(*,*) "Input sloppy convergence criterion(default=0.001D0)"
                 read(*,*) GMAX
               else if(GMAX_SW.eq.3)then
                 write(*,*) "Please input the range of the Convergence"
                 write(*,*) "From:(default=0.1D0)"
                 read(*,*) rms_start
                 write(*,*) "To:(default=0.00001D0)"
                 read(*,*) rms_end
                 write(*,*) "Input the variation rate(defaults=1.D0)"
                 read(*,*) rms_rate
               else 
                 call ERROR
                 goto 42
               endif
               write(*,*) "1.To use Gaussian Random Generator"
               write(*,*) "2.To use normal Random Generator(default)"
               read(*,*) gauss
43             write(*,*) "Would you like to take step?"
               write(*,*) "1.Take step!!"
               write(*,*) "2.Don't take step!!(default)"
               read(*,*) TAKESTEP
               if(TAKESTEP.eq.1)then
                 write(*,*) "Input the takestep parameter:"
                 write(*,*) "Input astep(default=0.4D0)"
                 read(*,*) ASTEP
                 write(*,*) "Input step(default=0.36D0)"
                 read(*,*) STEP
               else if(TAKESTEP.eq.2)then
               else
                 call ERROR
                 goto 43
               endif
               if(num.ne.0)then
30               write(*,*) "Would you like to reset to seeding number?"
                 write(*,*) "1.Reset to seeding number!!"
                 write(*,*) "2.Don't reset to seeding number!!(default)"
                 read(*,*) RESETSEED
                 if(RESETSEED.eq.1.or.RESETSEED.eq.2)then
                 else
                   call ERROR
                   goto 30
                 endif
               else
                 RESETSEED=2
               endif
             endif  
           goto 58
           else if (tmp.eq.2)then
59          write(*,*) "1.Change Basin-Hopping or PGA step numbers"
             write(*,*) "2.Change Multicanonical BH parameters"
             write(*,*) "3.Change Random move parameters"
             write(*,*) "4.Change Acceptance Ratio parameters"
             write(*,*) "5.Back to previous selection"
             read(*,*) tmp
             if (tmp.eq.1)then
51             write(*,*) "For Multi-canonical Monte Carlo,MC must>1"
               write(*,*) "1st MC RUN:"
               write(*,*) "1.Use Multi-canonical Basin-Hopping"
               write(*,*) "2.Use Canonical Basin-Hopping(default)"
               write(*,*) "3.USE Parallel GA"
               write(*,*) "4.USE Multi-canonical Annealing"
               read(*,*) BH_SW1
               if (BH_SW1.eq.1)then
         write(*,*) "Monte Carlo step number for 1st Run:(default=100)"
                 read(*,*) STEPMAX1
         write(*,*)"Multi-Canonical MC number for 1st Run:(default=30)"
                 read(*,*) MCMAX1
               else if(BH_SW1.eq.2.or.BH_SW1.eq.4)then
         write(*,*) "Monte Carlo step number for 1st Run:(default=1000)"
                 read(*,*) STEPMAX1
                 MCMAX1=1
               else if(BH_SW1.eq.3) then
      write(*,*)"Parallel GA generation number for 1st Run(default=500)"
                 read(*,*) STEPMAX1
                 MCMAX1=1
               else
                 call ERROR
                 goto 51
               endif
               ERANGE_SW=0
52             write(*,*) "2nd MC RUN:"
               if (BH_SW1.eq.4)then
                 ERANGE_SW=2
                 BH_SW2=5       !    For Aneealing purpose,reweighting!!
                 EMIN=0.D0
                 EMAX=0.D0
          write(*,*) "Monte Carlo step number for 2nd Run:(default=50)"
                 read(*,*) STEPMAX2
          write(*,*)"Multi-Canonical MC number for 2nd Run:(default=30)"
                 read(*,*) MCMAX2
               else
                 write(*,*) "1.Use Multi-canonical Basin-Hopping"
                 write(*,*) "2.Use Canonical Basin-Hopping(default)"
                 write(*,*) "3.Use Parallel GA"
                 read(*,*) BH_SW2
               endif
               if (BH_SW2.eq.1)then
           write(*,*) "Monte Carlo step number for 2nd Run:(default=50)"
                 read(*,*) STEPMAX2
          write(*,*)"Multi-Canonical MC number for 2nd Run:(default=30)"
                 read(*,*) MCMAX2
               else if(BH_SW2.eq.2)then
           write(*,*) "Monte Carlo step number for 2nd Run:(default=1)"
                 read(*,*) STEPMAX2
                 MCMAX2=1
               else if(BH_SW2.eq.3) then
      write(*,*)"Parallel GA generation number for 2nd Run(default=500)"
                 read(*,*) STEPMAX2
                 MCMAX2=2
               else if(BH_SW2.eq.4.or.BH_SW2.eq.5)then
               else
                 call ERROR
                 goto 52
               endif
53             write(*,*) "3rd MC RUN:"
               write(*,*) "1.Use Multi-canonical Basin-Hopping"
               write(*,*) "2.Use Canonical Basin-Hopping(default)"
               write(*,*) "3.Use Parallel GA" 
               read(*,*) BH_SW3
               if (BH_SW3.eq.1)then
          write(*,*) "Monte Carlo step number for 3rd Run:(default=100)"
               else if(BH_SW3.eq.2)then
         write(*,*) "Monte Carlo step number for 3rd Run:(default=1000)"
               else if(BH_SW3.eq.3)then
      write(*,*)"Parallel GA generation number for 3rd Run(default=500)"
               else
                 call ERROR
                 goto 53
               endif
               read(*,*) STEPMAX3
               if (BH_SW3.eq.1)then
          write(*,*)"Multi-Canonical MC number for 3rd Run:(default=30)"
                 read(*,*) MCMAX3
               else
                 MCMAX3=1
               endif
               write(*,*) "THE RANGE OF INTEGRAL OF THE ENERGY:"
               if(BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1)then
                 write(*,*) "0.INPUT ENERGY HISTOGRAM RANGE(DEFAULT)"
                 write(*,*) "1.USING AUTO MMC ENERGY RANGE(Half)"
                 read(*,*) ERANGE_SW
                 if(ERANGE_SW.eq.0)then
                   write(*,*) "UPPER ENERGY(EMAX):(default=?)"
                   read(*,*) EMAX
                   write(*,*) "LOWWER ENERGY(EMIN):(default=?)"
                   read(*,*) EMIN
                 endif
                 if(ERANGE_SW.eq.1.or.ERANGE_SW.eq.2)then
                   write(*,*) "Input the width of the energy histogram"
                   write(*,*) "From (default=2.D0)"
                   read(*,*) DEMIN
                   write(*,*) "To (default=5.D0)"
                   read(*,*) DEMAX
                 endif
               endif
               goto 59             
             else if(tmp.eq.2)then
      if(BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.BH_SW1.eq.4)then
               write(*,*) "You are using Multi-Canonical Basin-Hopping!"
                 if(BH_SW1.eq.4)then
                   write(*,*) "Using MMC Annealing!"
                 endif
                 write(*,*) "Please Input parameters:"
                 write(*,*) "# of Multicanonical BH BINS:(default=100)"
                 read(*,*) ibin
               else
              write(*,*) "You are using Single-Canonical Basin-Hopping!"
               endif
56        write(*,*) " Input initial temperature from:(default=0.08D0)"
               read(*,*) tempi
               write(*,*) " To :(default=0.8D0)"
               read(*,*) tempf
          write(*,*) "Please select what this range use for?(default=1)"
               write(*,*)"1.Use the temperature range for all ensembles"
               write(*,*)"2.Use the temperature range for the children"
               read(*,*) TEMP_EX
               write(*,*) " Input the quenching rate:(default=1.D0)"
               read(*,*) SCALEFAC
      if(BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.BH_SW1.eq.4)then
        write(*,*) "INITIAL TEMPFD:(default=BETA0=1/(",tempf,"-
     &TEMPFD))"
                 write(*,*) "(default TEMPFD=",(tempf-1000.D0),"D0)"
                 read(*,*) TEMPFD
                 BETA0=1.D0/(tempf-TEMPFD)
                 write(*,*) "BETA0=",BETA0
                 if (tempi-TEMPFD.le.0.D0) then
                   write(*,*) " ATTATION!!"
                   write(*,*) " Your current max temperature is ",tempf
                   write(*,*) " TEMPFD=",TEMPFD
                   write(*,*) " It is out of range!"
                   write(*,*) " Please select again!"
                   goto 56
                 endif
               endif
               goto 59
             else if(tmp.eq.3)then
               write(*,*) "VT>'ASTEP'*Vmin,ASTEP= (default=0.4D0)"
               read(*,*) ASTEP
               write(*,*) "The lower limit of ASTEP (default=1.D-4)"
               read(*,*) ASTEP_MIN
               write(*,*) "The upper limit of ASTEP (default<=1.D0)"
               read(*,*) ASTEP_MAX
           write(*,*) "The maximum change of any moving(default=0.36D0)"
               read(*,*) STEP
               write(*,*) "The lower limit of STEP (default=1.D-4)"
               read(*,*) STEP_MIN
               write(*,*) "The upper limit of STEP (default=",
     &dsqrt(enef)/2.D0,")"
               read(*,*) STEP_MAX
               write(*,*) "The temperature are now=",TEMPI," to ",TEMPF
               write(*,*) "The lower limit of TEMP (default=1.D-4)"
               read(*,*) TEMP_MIN
               write(*,*) "The upper limit of TEMP (default=100.D0)"
               read(*,*) TEMP_MAX
           write(*,*) "Input the adjustment of the Step,Astep and TEMP"
               write(*,*) "(default=1.005D0)"
               read(*,*) MOVE_FAC1
               goto 59
             else if(tmp.eq.4)then
      write(*,*) "The interval at which the acceptance ratio is checked"
               write(*,*) "(default=1)"
               read(*,*) naccept
               write(*,*) "Acceptance ratio for MC exploration:"
               write(*,*) "Please Note!!"
               write(*,*) "This value should dynamically change"
            write(*,*) "Each case should have its own Acceptance ratio"
               write(*,*)
     &"You have to check the hit rate(rate) while running the program!"
               write(*,*) "(default=0.5D0)"
               read(*,*) accrat
               write(*,*) "Input the adjustment of the Acceptance ratio"
               write(*,*) "(default=0.005D0)"
               read(*,*) MOVE_FAC2
               write(*,*) "1.Fix temperature(default)"
               write(*,*) "2.Do not fix temperature"
               read(*,*) FIXTEMP_T
               write(*,*) "1.Fix steps"
               write(*,*) "2.Do not fix steps(default)"
               read(*,*) FIXSTEP_T
               if (FIXTEMP_T.eq.1.and.FIXSTEP_T.eq.1) then
                 FIXBOTH_T=1
               else
                 FIXBOTH_T=2
               endif
               goto 59
             else
               goto 58
             endif
           else if(tmp.eq.3)then
             write(*,*) "1.Do not use Guided Function(default)"
             write(*,*) "2.Use Guided function"
             read(*,*) DG_SW
             if(DG_SW.ne.1)then
               write(*,*) "Please Input the Guided function resolution"
               write(*,*) "(default=50)"
               read(*,*) GBIN
      write(*,*) "1.Use temperature for Guided function,Now the T="
     &,TEMPF
      write(*,*) "2.Do not use temperature for Guided function(default)"
               read(*,*) FIX_GFAC
               if(FIX_GFAC.eq.1) then
                 GFAC=TEMPF
               else 
        write(*,*) "Input the acceptance probability(default=0.001D0)"
                 read(*,*) GFAC
               endif
             endif
             goto 58
           else 
             goto 2
           endif
***IC=5
         else if(ic.eq.5)then
           write(*,*) "1.Input the parameters of the Global GA"
           write(*,*) "2.Input the parameters of the Parallel GA"
           read(*,*) tmp
           if(tmp.eq.1)then
             write(*,*) "Input ensemble number: (default=20)"
             read(*,*) ens_num
178          write(*,*) "Input parents number in GA: (default=15)"
             read(*,*) par_num
             if(par_num.gt.ens_num)then
               write(*,*) "Unreasonable parent number !"
               goto 178
             endif
             write(*,*) "There will be ",ens_num," -",par_num,
     &"=",(ens_num-par_num)," children!!"
             write(*,*) 
176          write(*,*) "Input same number in GA(default=",ens_num+1,")"
             write(*,*) "Same number must between 2"," and ",ens_num+1
             read(*,*) samebc
             if(samebc.lt.2)then
               write(*,*) "Unreasonable same number !"
               goto 176
             endif
             write(*,*) "Input the max generation in GA: (default=500)"
             read(*,*) MAXGEN
             if(samebc.gt.ens_num)then
        write(*,*)"It will never stop until reach",MAXGEN," generations"
              write(*,*)
             endif
           else if(tmp.eq.2)then
             write(*,*) "Input the parameters of the Parallel GA"
           write(*,*) "1.Use the same parameters for all Runs(default)"
             write(*,*) "2.Use individual parameters for each Runs"
             read(*,*) tmp
80       write(*,*)"Input the ensemble number for 1st Run(default=20)"
             write(*,*) "(Max limitation=30)"
             read(*,*) sw_ensnum1
         write(*,*)"Parallel GA parents number for 1st Run(default=15)"
             read(*,*) sw_parnum1
         write(*,*)"Parallel GA same number for 1st Run(default=",
     &sw_ensnum1-sw_parnum1,")"
             read(*,*) sw_same1
             if(sw_parnum1.gt.sw_ensnum1.or.sw_same1.gt.sw_parnum1)then
               write(*,*) "WARNING!! INPUT ERROR,PLEASE INPUT AGAIN!!"
               goto 80
             endif
             if(tmp.eq.1)then
               write(*,*) "Use the same parameters!"
               sw_ensnum2=sw_ensnum1
               sw_parnum2=sw_parnum1
               sw_same2=sw_same1
               sw_ensnum3=sw_ensnum1
               sw_parnum3=sw_parnum1
               sw_same3=sw_same1
             else if(tmp.eq.2)then
81       write(*,*)"Input the ensemble number for 2nd Run(default=20)"
               write(*,*) "(Max limitation=30)"
               read(*,*) sw_ensnum2
         write(*,*)"Parallel GA parents number for 2nd Run(default=15)"
               read(*,*) sw_parnum2
         write(*,*)"Parallel GA same number for 2nd Run(default=",
     &sw_ensnum2-sw_parnum2,")"
               read(*,*) sw_same1
              if(sw_parnum2.gt.sw_ensnum2.or.sw_same2.gt.sw_parnum2)then
                 write(*,*) "WARNING!! INPUT ERROR,PLEASE INPUT AGAIN!!"
                 goto 81
               endif
82       write(*,*)"Input the ensemble number for 3rd Run(default=20)"
               write(*,*) "(Max limitation=30)"
               read(*,*) sw_ensnum3
         write(*,*)"Parallel GA parents number for 3rd Run(default=15)"
               read(*,*) sw_parnum3
         write(*,*)"Parallel GA same number for 3rd Run(default=",
     &sw_ensnum3-sw_parnum3,")"
               read(*,*) sw_same3
              if(sw_parnum3.gt.sw_ensnum3.or.sw_same3.gt.sw_parnum3)then
                 write(*,*) "WARNING!! INPUT ERROR,PLEASE INPUT AGAIN!!"
                 goto 82
               endif
             endif
           endif
           goto 2
***IC=6
         else if(ic.eq.6)then
57         write(*,*) "1.Input the ratio of the genetic operator"
           write(*,*) "2.Input the ratio of the step moving"
           write(*,*) "3.Back to main menu"
           read(*,*) tmp
           if (tmp.eq.1) then
             write(*,*) "Input the ratio of the genetic operator:"
             write(*,*) "Input the ratio of Inversion: (default=5)"
             read(*,*) aa
             write(*,*) "Input the ratio of Arithmetic: (default=1)"
             read(*,*) bb 
             write(*,*) "Input the ratio of Geometic: (default=1)"
             read(*,*) cc
             write(*,*) "Input the ratio of Crossing: (default=5)"
             read(*,*) dd
             write(*,*) "Input the ratio of 2-Point: (default=5)"
             read(*,*) ee
             write(*,*) "Input the ratio of 3N Mutation: (default=0)"
             read(*,*) ff
             write(*,*) "Input the ratio of Moment op 1:(default=0)"
             read(*,*) gg
             write(*,*) "Input the ratio of Moment op 2:(default=0)"
             read(*,*) hh
             write(*,*) "Input the ratio of DIH move op 1:(default=0)"
             read(*,*) ii
             write(*,*) "1.Use eularangle transformation"
            write(*,*) "2.Do not use eularangle transformation(default)"
             read(*,*) sw_eular1
             if(sw_eular1.eq.1)then
               write(*,*) "Input the factor of the eularangle trans"
               write(*,*) "(Eular angle=Pei*factor*rnd number)"
               write(*,*) "(default=0.5D0)"
               read(*,*) eular_fac1
             endif
             goto 57
           else if (tmp.eq.2) then
60           write(*,*) "Input the ratio of the step moving:"
             write(*,*) "1.Use the same ratios for all MC runs(default)"
             write(*,*) "2.Use seperate ratio for each MC runs"
             read(*,*) tmp
             write(*,*) "Input the ratio for 1st MC run"
        write(*,*) "Input the ratio of 3N particles moving :(default=5)"
             read(*,*) mva1
      write(*,*) "Input the ratio of 3N particles mutation :(default=0)"
             read(*,*) mvb1
        write(*,*) "Input the ratio of N dimension moving :(default=0)"
             read(*,*) mvc1
        write(*,*) "Input the ratio of 3N dihedral moving :(default=0)"
             read(*,*) mve1
        write(*,*) "Input the ratio of RMS transition :(default=0)"
             read(*,*) mvf1
             if (mvf1.ne.0.D0)then
               write(*,*) "Input the RMS scale :(default=0.01D0)"
               read(*,*) rms_fac1
             endif
             write(*,*) "Input the ratio of Inversion:(default=5)"
             read(*,*) mvd1
             write(*,*) "Input the ratio of Arithmetic:(default=1)"
             read(*,*) mvg1
             write(*,*) "Input the ratio of Geometic:(default=1)"
             read(*,*) mvh1
             write(*,*) "Input the ratio of Crossing:(default=5)"
             read(*,*) mvi1
           write(*,*) "Input the ratio of 2-Point Crossover:(default=5)"
             read(*,*) mvj1
             if(tmp.eq.1)then
               mva2=mva1
               mva3=mva1
               mvb2=mvb1
               mvb3=mvb1
               mvc2=mvc1
               mvc3=mvc1
               mvd2=mvd1
               mvd3=mvd1
               mve2=mve1
               mve3=mve1
               mvf2=mvf1
               mvf3=mvf1
               mvg2=mvg1
               mvg3=mvg1
               mvh2=mvh1
               mvh3=mvh1
               mvi2=mvi1
               mvi3=mvi1
               mvj2=mvj1
               mvj3=mvj1
               rms_fac2=rms_fac1
               rms_fac3=rms_fac1
             else if(tmp.eq.2)then
               write(*,*) "Input the ratio for 2nd MC run"
      write(*,*) "Input the ratio of 3N particles moving :(default=5)"
               read(*,*) mva2
      write(*,*) "Input the ratio of 3N particles mutation :(default=0)"
               read(*,*) mvb2
        write(*,*) "Input the ratio of N dimension moving :(default=0)"
               read(*,*) mvc2
        write(*,*) "Input the ratio of 3N dihedral moving :(default=0)"
               read(*,*) mve2
        write(*,*) "Input the ratio of RMS transition :(default=0)"
               read(*,*) mvf2
               if (mvf2.ne.0.D0)then
                 write(*,*) "Input the RMS scale :(default=0.01D0)"
                 read(*,*) rms_fac2
               endif
               write(*,*) "Input the ratio of Inversion:(default=5)"
               read(*,*) mvd2
               write(*,*) "Input the ratio of Arithmetic:(default=1)"
               read(*,*) mvg2
               write(*,*) "Input the ratio of Geometic:(default=1)"
               read(*,*) mvh2
               write(*,*) "Input the ratio of Crossing:(default=5)"
               read(*,*) mvi2
           write(*,*) "Input the ratio of 2-Point Crossover:(default=5)"
               read(*,*) mvj2
               write(*,*) "Input the ratio for 3rd MC run.(Final run)"
      write(*,*) "Input the ratio of 3N particles moving :(default=5)"
               read(*,*) mva3
      write(*,*) "Input the ratio of 3N particles mutation :(default=0)"
               read(*,*) mvb3
        write(*,*) "Input the ratio of N dimension moving :(default=0)"
               read(*,*) mvc3
        write(*,*) "Input the ratio of 3N dihedral moving :(default=0)"
               read(*,*) mve3
        write(*,*) "Input the ratio of RMS transition :(default=0)"
               read(*,*) mvf3
               if (mvf3.ne.0.D0)then
                 write(*,*) "Input the RMS scale :(default=0.01D0)"
                 read(*,*) rms_fac3
               endif
               write(*,*) "Input the ratio of Inversion:(default=5)"
               read(*,*) mvd3
               write(*,*) "Input the ratio of Arithmetic:(default=1)"
               read(*,*) mvg3
               write(*,*) "Input the ratio of Geometic:(default=1)"
               read(*,*) mvh3
               write(*,*) "Input the ratio of Crossing:(default=5)"
               read(*,*) mvi3
           write(*,*) "Input the ratio of 2-Point Crossover:(default=5)"
               read(*,*) mvj3
             else
               call ERROR
               goto 60 
             endif
             write(*,*) "1.Use eularangle transformation"
            write(*,*) "2.Do not use eularangle transformation(default)"
             read(*,*) sw_eular2
             if(sw_eular2.eq.1)then
               write(*,*) "Input the factor of the eularangle trans"
               write(*,*) "(Eular angle=Pei*factor*rnd number)"
               write(*,*) "(default=0.5D0)"
               read(*,*) eular_fac2
             endif
             goto 57
           else
             goto 2
           endif
***IC=7
         else if(ic.eq.7)then
           write(*,*) "Which one do you want to change ? :"
           write(*,*) "1.Change the Initial condition :"
           write(*,*) "2.Change the parameter of energy :"
           read(*,*) i
c---------------------
           if(i.eq.1)then
71          write(*,*) "Do you want to use the close_pack radius? (y/n)"
             write(*,*) "(default='y',)"
             read(*,*) load
c---------------------
             if(load.eq."y")then
               d=1
             else if(load.eq."n")then
              write(*,*) "Input the initial width of box:(default=5.5D)"
               write(*,*) "This value will change with particle number!"
               read(*,*) LL
               LL=LL*(real(atom_num)**(1.D0/3.D0))
               d=0
             else
               call ERROR
               goto 71
             endif
c---------------------
             goto 2
           else if(i.eq.2)then
72           write(*,*)" 1.Many-body potential"
             write(*,*)" 2.Alloys with many-body potential"
             write(*,*)" 3.Protein-folding:Beta folding"
             writE(*,*)" 4.Protein-folding:Alpha folding"
             write(*,*)" 5.EAM Potential  "
             write(*,*)" 6.Universal protein folding"
             write(*,*)" 7.Chen and Imamura's protein folding(N)"
             write(*,*)" 8.Chen and Imamura's protein folding(3N)"
             read(*,*) type
             if(type.eq.1.or.type.eq.2.or.type.eq.3.or.type.eq.4
     &.or.type.eq.5.or.type.eq.6.or.type.eq.7.or.type.eq.8)then
             else
               call ERROR
               goto 72
             endif
C            <<<START MANY-BODY>>>
             if(type.eq.1.or.type.eq.2)then
               dim3=0
               open(77,file="enrgpmt.dat",status="old")
               z=0
               kind0="mwz"
               write(*,*) "@      e0       c0      p    q    r0"
               do while(kind0.ne."end")
                 z=z+1
                 read(77,*) kind0,epsilon0,zeta0,p0,q0,rzero0
                 if(kind0.ne."end")then
                   kind=kind0
                   epsilon=epsilon0
                   zeta=zeta0
                   p=p0
                   q=q0
                   rzero=rzero0
                   write(*,*) kind,epsilon,zeta,p,q,rzero
                 endif
                 if(z.eq.10.or.z.eq.20.or.z.eq.30.or.z.eq.40.or.z
     &.eq.50.or.z.eq.60.or.z.eq.70.or.z.eq.80.or.z.eq.90
     &.or.z.eq.100.or.z.eq.110)then
                   pause
                 endif
               enddo
               close(77)
               if(type.eq.1)then
                 write(*,*) "Please input a materials: [default='Na']"
                 write(*,*) "Or you can type 'none' to input directly"
                 read(*,*) load  
               endif
               if(type.eq.2)then 
                 write(*,*) "Please input the name of A atom:"
                 read(*,*) load
                 write(*,*) "Please input the name of B atom:"
                 read(*,*) load2
                 write(*,*) "Please input the # of A atoms:"
                 read(*,*) A
               endif
               if(load.eq."none")then
                 write(*,*) "Input the parameter of the energy: (MC)"
                 write(*,*) "Input epsilon: [for NA:0.015948 (eV)]"
                 read(*,*) epsilon
                 write(*,*) "Input zeta: [for NA:0.29113 (eV)]"
                 read(*,*) zeta
                 write(*,*) "Input p: [for NA: p=10.13]"
                 read(*,*) p
                 write(*,*) "Input q: [for NA: q=1.30]"
                 read(*,*) q
                 write(*,*) "Input rzero: [for NA: r0=6.99 a.u]"
                 read(*,*) rzero
                 goto 2
               endif
c----------- Energy parameters DATABASE ! -------------
               kind0="mzw"
               open(68,file="enrgpmt.dat",status="old")
               kind4=load
               kind5=load2
               do while(kind0.ne."end")
                 read(68,*) kind0,epsilon0,zeta0,p0,q0,rzero0
                 if(kind0.eq.load)then
                   kind=kind0
                   epsilon=epsilon0
                   zeta=zeta0
                   p=p0
                   q=q0
                   rzero=rzero0
                   write(*,*) "@      e0       c0      p    q   r0"
                   write(*,*) kind,epsilon,zeta,p,q,rzero
                   if(type.eq.1) goto 2
                 endif
                 if(type.eq.2)then
                   if(kind0.eq.load2)then
                     kind2=kind0
                     epsilon2=epsilon0
                     zeta2=zeta0
                     p2=p0
                     q2=q0
                     rzero2=rzero0
                     write(*,*) "@      e0       c0      p    q   r0"
                     write(*,*) kind2,epsilon2,zeta2,p2,q2,rzero2
                   endif
                  if(kind0.eq.kind4//kind5.or.kind0.eq.kind5//kind4)then
                     kind3=kind0
                     epsilon3=epsilon0
                     zeta3=zeta0
                     p3=p0
                     q3=q0
                     rzero3=rzero0
                     write(*,*) "@      e0       c0      p    q   r0"
                     write(*,*) kind3,epsilon3,zeta3,p3,q3,rzero3
                     goto 2    
                   endif 
                 endif
                 if(kind0.eq."end")then
                   write(*,*) "Data doesn't exist !"
                   write(*,*) "Please type 'go' and try again !"
                   pause
                   close(68)
                   goto 2
                 endif
               enddo
               close(68)
             endif
C            <<<END OF MANY-BODY>>>
             if(type.eq.3.or.type.eq.4.or.type.eq.6)then
               dim3=0
31             write(*,*) "1.Edit residues"
               write(*,*) "2.Edit potential parameters"
               write(*,*) "3.See reference"
               read(*,*) load3
               if (load3.eq.3)then
                 write(*,*) "Input the residue sets:"
                 write(*,*) "ex:"
                 write(*,*) "total residues= 0,residue type:"
                 write(*,*) "(input residues type,ex:B)"
                 write(*,*) "#="
                 write(*,*) "(input residues length,ex:3)"
                 write(*,*) "total residues= 3,residue type:"
                 write(*,*) "(input next residues until stop)"
                 write(*,*) "etc......"
                 goto 31
               else if (load3.eq.1)then
                 write(*,*) "Starting input residues:"
                 res_dim=0
                 ndim_test=0
                 do while(res_dim.lt.atom_num)
                   ndim_test=ndim_test+1
                   write(*,*) "total residues=",res_dim,",residue type:"
                   read(*,*) prot(ndim_test)
                   write(*,*) "#="
                   read(*,*) residues_num(ndim_test)
                   if (prot(ndim_test).eq."LB")then
                     res_dim=res_dim+(residues_num(ndim_test)*2)
                   else if (prot(ndim_test).eq."BBL".or.prot(ndim_test)
     &.eq."BLB".or.prot(ndim_test).eq."LBB")then
                     res_dim=res_dim+(residues_num(ndim_test)*3)
                   else
                     res_dim=res_dim+residues_num(ndim_test)
                   endif
                 enddo
                 write(*,*) "total residues=",res_dim,",STOP INPUT!"
               else if (load3.eq.2)then
                 write(*,*) "Input cutoff length in LJ:(default=5.5D0)"
                 read(*,*) CUTOFF
                 write(*,*) "Input E in LJ:(default=1.D0)"
                 read(*,*) ECOIL
             write(*,*) "Input theta0 in Bond energy:(default=1.8326D0)"
                 read(*,*) THETA
                 write(*,*) 
     &"Input cutoff angle in Bond energy:(default=0.087266D0)"
                 read(*,*) CUTOFF_ANGLE
                 if(type.eq.6)then
                  write(*,*) "Input eps in Dihedral energy:(default=?)"
                  read(*,*) EPS
                  write(*,*) "Input eps2 in Dihedral energy:(default=?)"
                  read(*,*) EPS2
                  write(*,*) "Input eps3 in Dihedral energy:(default=?)"
                  read(*,*) EPS3
                  write(*,*) "Input eps4 in Dihedral energy:(default=?)"
                  read(*,*) EPS4
                 endif
                 if(type.eq.3.or.type.eq.4)then
                 write(*,*) "Input A in Dihedral energy:(default=1.2D0)"
                   read(*,*) ACOIL
                 write(*,*) "Input B in Dihedral energy:(default=1.2D0)"
                   read(*,*) BCOIL
                   write(*,*) "Input C in Dihedral energy:"
                   write(*,*) "(default:Alpha=1.2D0/Beta=0.D0)"
                   read(*,*) CCOIL
                 write(*,*) "Input D in Dihedral energy:(default=0.2D0)"
                   read(*,*) DCOIL
                 endif
               else
                 call ERROR
                 goto 31
               endif
               goto 2
             endif  
             if(type.eq.5)then
c                dim3=0
c               open(58,file="eamenrg.dat",status="old")
c               z=0
c              kind0="mwz"
c               write(*,*) "@  eam_a  eam_e  eam_C  eam_n  eam_m"
c               do while(kind0.ne."end")
c                 z=z+1
c                 read(58,*) kind0,eam_a1,eam_epsilon1,eam_C1,
c     &eam_n1,eam_m1
c                 if(kind0.ne."end")then
c                  kind=kind0
c                   eam_a=eam_a1
c                   eam_epsilon=eam_epsilon1
c                   eam_C=eam_C1
c                 eam_n=eam_n1c
c                 eam_m=eam_m1
c                  write(*,*) kind,eam_a,eam_epsilon,eam_C,eam_n,eam_m
c                endif
c               if(z.eq.10.or.z.eq.20.or.z.eq.30.or.z.eq.40.or.zc
c    &.eq.50.or.z.eq.60.or.z.eq.70.or.z.eq.80.or.z.eq.90
c    &.or.z.eq.100.or.z.eq.110)then
c                  pause
c                endif
c              enddo               close(58)
                 write(*,*) "Please input a materials: [default='Na']"
                 write(*,*) "Or you can type 'none' to input directly"
                 read(*,*) load 
               open(57,file="eamenergy.dat",status="old") 
               do while(kind0.ne."end")
                 read(57,*) kind0,eam_a,eam_epsilon,eam_C,
     &eam_n,eam_m
                 if(kind0.eq.load)then
                   kind=kind0
                   eam_a1=eam_a
                   eam_epsilon1=eam_epsilon
                   eam_C1=eam_C
                   eam_n1=eam_n
                   eam_m1=eam_m
                   write(*,*) "@  eam_a  eam_e  eam_C  eam_n  eam_m"
                   write(*,*) kind,eam_a1,eam_epsilon1,eam_C1,
     &eam_n1,eam_m1
                 endif
                 goto 2    
                 if(kind0.eq."end")then
                   write(*,*) "Data doesn't exist !"
                   write(*,*) "Please type 'go' and try again !"
                   pause
                   close(57)
                   goto 2
                 endif
               enddo
            endif
C      Chen amd Imamura's protein folding!           
47           if (type.eq.7.or.type.eq.8)then
               if(type.eq.7)then
                 dim3=1
               endif
               if(type.eq.8)then
                 dim3=0
               endif
74             write(*,*) "1.Edit residues" 
               write(*,*) "2.Edit potential parameters"
               read(*,*) load3
               if (load3.eq.1)then
             write(*,*) "Please input the number of neutral bonds set:"
               write(*,*) 
     &"(Could include several residues;default=0;No neutral bonds)"
                 read(*,*) ndim_test2
                 if (ndim_test2.eq.0)then
                   ndim_test2=1
                   ch_st(1)=0
                   ch_end(1)=0
                   goto 2
                 endif
                 write(*,*) "Please make sure no overlaping!"
                 write(*,*) 
                 do j1=1,ndim_test2
      write(*,*)"The ",j1,"'th neutral bonds set start from residue:"
                   read(*,*) ch_st(j1)
                   write(*,*) "End from:"
                   read(*,*) ch_end(j1)
                   if (ch_st(j1).gt.ch_end(j1))then
                     write(*,*) "ERROR!!"
                   write(*,*) "Start number is greater than end number"
                     write(*,*) "Please type again!!"
                     write(*,*) 
                     goto 47
                   endif
                 enddo
                 do j1=1,ndim_test2
                   write(*,*) 
     &j1,"'th neutral set start from ",ch_st(j1),",end from ",ch_end(j1)
                 enddo
               else if(load3.eq.2)then
       write(*,*)"Please input Lo,the bond distance:(default=1.5D0)"
                 read(*,*) FCOIL
       write(*,*)"Please input exclusive diameter:(default=1.D0)"
                 read(*,*) EXC_D  
       write(*,*) "Please input epsilon,the energy scale:(default=1.D0)"
                 read(*,*) GCOIL
       write(*,*) "Please input the bond angle:(default=105.D0)"
                 read(*,*) THETA1
                 THETA1=(3.14159D0*THETA1)/180.D0
       write(*,*) "Your bond-angle is ",THETA1," rad"
               else
                 call ERROR
                 goto 74
               endif
               goto 2
             endif
c---------------
           else
             write(*,*) "Please make a selection again!"
             goto 2
           endif
***IC=9
         else if(ic.eq.9)then
           if (min.eq.5)then
             enef=1.D0+(3.0D0*atom_num/17.77153175D0)**(1.0D0/3.0D0)
             enef=enef*2.0D0**(1.0D0/6.0D0)
             enef=enef*enef
             if (d.eq.1) LL=enef
           endif                     
           write(*,*) "@ ",atom_num," particles"
           write(*,*) "@ seeding number= ",num,",(0 for no seeding)"
           if(gauss.eq.1)then
             write(*,*) "@ using Gaussian Random number"
           endif
           if(min.eq.1)then
             write(*,*) "@ using Simplex minimization"
           else if(min.eq.2)then
             write(*,*) "@ using Conjugate Gradient minimization"
           else if(min.eq.3)then
             write(*,*) "@ using Conjugate Gradient then Simplex"
           else if(min.eq.4)then
             write(*,*) "@ using Simplex then Conjugate Gradient"
           else
             write(*,*) "@ using LBFGS minimization,MAXIT=",MAXIT
             if(TAKESTEP.eq.1)then
               write(*,*) "@ Take step!!"
             else
               write(*,*) "@ Don't take step!!"
             endif
             if(RESETSEED.eq.1)then
               write(*,*) "@ Reset to seeding number!!"
             else
               write(*,*) "@ Don't reset to seeding number!!"
             endif
             if(cent2.eq.1)then
               write(*,*) "@ move coords to center of mass!"
             else
               write(*,*) "@ Do not move coords to center of mass!"
             endif
           endif
           write(*,*) "@ ",ens_num," ensembles for global GA"
           write(*,*) "@ ",sw_ensnum1," ensembles for parallel GA run 1"
           write(*,*) "@ ",sw_ensnum2," ensembles for parallel GA run 2"
           write(*,*) "@ ",sw_ensnum3," ensembles for parallel GA run 3"
           if (sw_eular1.eq.1)then
           write(*,*) "@ Use eularangle tran for GGA,factor=",eular_fac1
           else
             write(*,*) "@ Do not use eularangle tran for GGA"
           endif
           if (sw_eular2.eq.1)then
           write(*,*) "@ Use eularangle tran for PGA,factor=",eular_fac2
           else
             write(*,*) "@ Do not use eularangle tran for PGA"
           endif
           write(*,*) "@ The ratio of GA:"
           write(*,*) "@ Inversion:    ",aa
           write(*,*) "@ Arithmetic:   ",bb
           write(*,*) "@ Geometic:     ",cc
           write(*,*) "@ Crossing:     ",dd
           write(*,*) "@ 2-Point:      ",ee
           write(*,*) "@ 3N Mutation:  ",ff
           write(*,*) "@ Moment op 1:  ",gg
           write(*,*) "@ Moment op 2:  ",hh
           write(*,*) "@ DIH move op 1:",ii
           write(*,*) "@ The ratio of Step move:"
           write(*,*) "@ 3N Particles Moving,1st=",mva1,",2nd=",mva2,
     &",3rd=",mva3
         write(*,*) "@ 3N Particles Mutation,1st=",mvb1,",2nd=",mvb2,
     &",3rd=",mvb3
         write(*,*) "@ N Dimension Moving,   1st=",mvc1,",2nd=",mvc2,
     &",3rd=",mvc3
         write(*,*) "@ 3N Dihedral Moving,   1st=",mve1,",2nd=",mve2,
     &",3rd=",mve3
         write(*,*) "@ RMS Transition,       1st=",mvf1,",2nd=",mvf2,
     &",3rd=",mvf3
         write(*,*) "@ RMS SCALE,            1st=",rms_fac1,",2nd=",
     &rms_fac2,",3rd=",rms_fac3
           write(*,*) "@ Inversion,          1st=",mvd1,",2nd=",mvd2,
     &",3rd=",mvd3
           write(*,*) "@ Arithmetic,         1st=",mvg1,",2nd=",mvg2,
     &",3rd=",mvg3
           write(*,*) "@ Geometic,           1st=",mvh1,",2nd=",mvh2,
     &",3rd=",mvh3
           write(*,*) "@ Crossing,           1st=",mvi1,",2nd=",mvi2,
     &",3rd=",mvi3
         write(*,*) "@ 2-Point Crossover,    1st=",mvj1,",2nd=",mvj2,
     &",3rd=",mvj3
           write(*,*) "@",par_num," parents for global GA"
           write(*,*) "@",sw_parnum1," parents for parallel GA run 1"
           write(*,*) "@",sw_parnum2," parents for parallel GA run 2"
           write(*,*) "@",sw_parnum3," parents for parallel GA run 3"
           write(*,*) "@ Same number=",samebc," for global GA"
           write(*,*) "@ Same number=",sw_same1," for parallel GA run 1"
           write(*,*) "@ Same number=",sw_same2," for parallel GA run 2"
           write(*,*) "@ Same number=",sw_same3," for parallel GA run 3"
           write(*,*) "@ Max gens=",MAXGEN
           write(*,*) "@ Multi canonical max number for 1st Run=",MCMAX1
        write(*,*) "@ Max MC steps for each sweep for 1st Run=",STEPMAX1
           write(*,*) "@ Multi canonical max number for 2nd Run=",MCMAX2
        write(*,*) "@ Max MC steps for each sweep for 2nd Run=",STEPMAX2
           write(*,*) "@ Multi canonical max number for 3rd Run=",MCMAX3
        write(*,*) "@ Max MC steps for each sweep for 3rd Run=",STEPMAX3
           if (d.eq.1)then
             write(*,*) "@ Using close pack condition"
           endif
           write(*,*) "@ Width of initial box is ",LL
           write(*,*) "@ BC in func=",enef
           write(*,*) "@ Accuracy parameter in func=",enef2
           write(*,*) "@ The energy sort accuracy is ",ene
           write(*,*) "@ The dX in dfunc is ",dx
           write(*,*) "@ BC in dfunc is ",bc                     
           write(*,*) "@ Tolerance of Simplex is ",simp       
           write(*,*) "@ Tolerance of Conjugate Gradient is ",cg
           write(*,*) "@ Tolerance of LBFGS(GTOL)=",GTOL
           if(GMAX_SW.eq.1)then
      write(*,*) "@ Using RMS force for sloppy convergence criterion"
             write(*,*) "@ Threshold criterion=",GMAX
           else if(GMAX_SW.eq.2)then
      write(*,*) "@ The sloppy convergence criterion for RMS force=",
     &GMAX
           else if(GMAX_SW.eq.3)then
      write(*,*) "@ The multi-convergence criterion from ",rms_start,
     &" to ",rms_end
             write(*,*) "@ The variation rate is ",rms_rate
           endif
           write(*,*) "@ MCSCRH UPDATE =",MUPDATE
           write(*,*) "@ Takestep in bmin:Astep=",ASTEP,",Step=",STEP
           write(*,*) "@ Max Step=",STEP_MAX,",Min Step=",STEP_MIN
           write(*,*) "@ Max Astep=",ASTEP_MAX,",Min Astep=",ASTEP_MIN
           write(*,*) "@ Max TEMP=",TEMP_MAX,",Min TEMP=",TEMP_MIN
        write(*,*) "@ The adjustment of the Acceptance checking is +-",
     &MOVE_FAC1

           write(*,*) "@ Input DGUESS=",DGUESS
           if(DG_SW.eq.1)then
             write(*,*) "@ Do not use Guided function"
           else
             write(*,*) "@ Use Guided function"
             write(*,*) "@ # of Guided Function BINS=",GBIN
           endif
           if (FIX_GFAC.eq.1)then
             write(*,*) "@ Using temperature for Guided function=",GFAC
           else
             write(*,*) "@ The acceptance probability for Guided fun is"
     &,GFAC
           endif
           write(*,*) "@ # of Multicanonical BH BINS=",ibin
        if(BH_SW1.eq.1.or.BH_SW1.eq.4.or.BH_SW2.eq.1.or.BH_SW3.eq.1)then
            write(*,*) "@ INITIAL BETA0=1/(",TEMPF,"-",TEMPFD,")=",BETA0
             if(ERANGE_SW.eq.0)then
               write(*,*) "@ ENERGY RANGE FROM ",EMAX," TO ",EMIN
             else if(ERANGE_SW.eq.1)then
               write(*,*) "@ Using Auto Energy Range(Half)"
               write(*,*) "@ TDE from ",DEMIN," to ",DEMAX
             else if(ERANGE_SW.eq.2)then
               write(*,*) "@ Using Multi-Canonical Annealing BH"
               write(*,*) "@ TDE from ",DEMIN," to ",DEMAX
             endif
           else
             write(*,*) "@ Using MonteCarlo Boltzmann distribution"
           endif
      write(*,*)"@ Every ",naccept," generations check acceptance ratio"
           write(*,*)"@ Acceptance ratio=",accrat
           write(*,*) "@ The adjustment of the Acceptance Ratio is +-",
     &MOVE_FAC2    
      write(*,*) "@ Temperature start from ",TEMPI," to ",TEMPF,
     &",multiply by ",SCALEFAC
           if(TEMP_EX.eq.1)then
             write(*,*) "@ Use temperature range for all ensembles!"
           else if(TEMP_EX.eq.2)then
             write(*,*) "@ Use temperature range for the children!"
           else
             write(*,*) "ERROR!!"
             return
           endif
           if(FIXTEMP_T.eq.1.or.FIXBOTH_T.eq.1)then
             write(*,*) "@ FIX TEMPERATURE!"
           else
             write(*,*) "@ DO NOT FIX TEMPERATURE!"
           endif
           if(FIXSTEP_T.eq.1.or.FIXBOTH_T.eq.1)then
             write(*,*) "@ FIX STEPS!"
           else  
             write(*,*) "@ DO NOT FIX STEPS!"
           endif
           if(rec_sw.eq.0)then
             write(*,*) "@ Normal recording"
           endif
           if(rec_sw.eq.1)then
             if (pdb_rec.eq.1)then
               write(*,*) "@ Step recording in PDB file!"
             else
               write(*,*) "@ Step recording in seperate coords file!"
             endif
             write(*,*) "@ Step recording,recording ratio=",rec_ratio
           endif
           if(rec_sw.eq.2)then
             write(*,*) "@ Record nothing!"
           endif
           if(recover_rec.eq.1)then
             write(*,*) "@ Start recover,recover ratio=",recover_ratio
           endif
           if (type.eq.1)then
             write(*,*) "@     e0        c0       p        q    r0"
             write(*,*) "@ ",epsilon,",",zeta,",",p,",",q,",",rzero
             write(*,*) "@ For ",kind
           endif
           if (type.eq.2)then
             write(*,*) "@     e0        c0       p        q    r0"
             write(*,*) "@ ",epsilon,",",zeta,",",p,",",q,",",rzero
             write(*,*) "@ For ",A," ",kind
             write(*,*) "@     e0        c0       p        q    r0"
         write(*,*) "@ ",epsilon2,",",zeta2,",",p2,",",q2,",",rzero2
             write(*,*) "@ For ",atom_num-A," ",kind2         
             write(*,*) "@     e0        c0       p        q    r0"
         write(*,*) "@ ",epsilon3,",",zeta3,",",p3,",",q3,",",rzero3
             write(*,*) "@ For ",kind3
           endif
           if (type.eq.3)then
             write(*,*) "@ For Beta-Sheets"
             write(*,*) "@ ",(prot(j1),residues_num(j1),j1=1,ndim_test)
           endif
           if (type.eq.4)then
             write(*,*) "@ For Alpha-Helices"
             write(*,*) "@ ",(prot(j1),residues_num(j1),j1=1,ndim_test)
           endif
           if (type.eq.5)then
              write(*,*) "@  eam_a  eam_e  eam_C  eam_n  eam_m"
              write(*,*) "@ ",eam_a1,",",eam_epsilon1,",",eam_C1,","
     &,eam_n1,",",eam_m1
              write(*,*) "@ For ",kind
           endif
           if (type.eq.7)then
             write(*,*) "@ For Chen and Imamura Folding N"
             write(*,*) "@ The Neutral bond sets:(Start-end)"
          write(*,*) "@ ",("(",ch_st(j1),ch_end(j1),")",j1=1,ndim_test2)
             write(*,*) "@ Lo=",FCOIL,",Epsilon=",GCOIL,",ANGLE=",THETA1
     &,",Exclusive diameter=",EXC_D
           endif
           if (type.eq.8)then
             write(*,*) "@ For Chen and Imamura Folding 3N"
             write(*,*) "@ The Neutral bond sets:(Start-end)"
          write(*,*) "@ ",("(",ch_st(j1),ch_end(j1),")",j1=1,ndim_test2)
             write(*,*) "@ Lo=",FCOIL,",Epsilon=",GCOIL,",ANGLE=",THETA1
     &,",Exclusive diameter=",EXC_D
           endif


           if (asw.eq.1)then
             write(*,*) "@ Using analytic gradient"
           else
             write(*,*) "@ Using numerical gradient"
           endif
           write(*,*)
           write(*,*) "                                         "
           write(*,*) "   Quit and Save the config. file ? (y/n)"
           write(*,*) "                                         "
           read(*,*) load
           write(*,*)
           if(load.eq."y")then
             open(20,file="data.dat")
             write(20,*) LL,ene,enef,enef2,simp,cg,gauss,d
       write(20,*) ens_num,min,aa,bb,cc,dd,ee,ff,gg,hh,ii,par_num,samebc
     &,MAXGEN
         write(20,*) dx,bc,asw,sw_eular1,sw_eular2,eular_fac1,eular_fac2
      write(20,*) MAXIT,MUPDATE,TAKESTEP,RESETSEED,DGUESS,ASTEP,STEP,
     &GTOL,GMAX,GMAX_SW
             write(20,*) rec_sw,rec_ratio,pdb_rec,dim3
             write(20,*) recover_sw,recover_rec,recover_ratio
             write(20,*)
     &ibin,BETA0,N0,TEMPFD,EMAX,EMIN,naccept,accrat
         write(20,*) FIXTEMP_T,FIXSTEP_T,FIXBOTH_T,SCALEFAC,TEMPI,TEMPF
     &,TEMP_EX
           write(20,*) mva1,mvb1,mvc1,mvd1,mve1,mvf1,mvg1,mvh1,mvi1,mvj1
           write(20,*) mva2,mvb2,mvc2,mvd2,mve2,mvf2,mvg2,mvh2,mvi2,mvj2
           write(20,*) mva3,mvb3,mvc3,mvd3,mve3,mvf3,mvg3,mvh3,mvi3,mvj3
             write(20,*) GBIN,GFAC,FIX_GFAC,DG_SW,DG_RECOVER 
             write(20,*) STEPMAX1,MCMAX1,STEPMAX2,MCMAX2,STEPMAX3,MCMAX3
             write(20,*) BH_SW1,BH_SW2,BH_SW3,DEMIN,DEMAX,ERANGE_SW
       write(20,*) rms_start,rms_end,rms_rate,rms_fac1,rms_fac2,rms_fac3
             write(20,*) sw_ensnum1,sw_parnum1,sw_same1
             write(20,*) sw_ensnum2,sw_parnum2,sw_same2
             write(20,*) sw_ensnum3,sw_parnum3,sw_same3
             write(20,*) MOVE_FAC1,MOVE_FAC2,TEMP_MIN,TEMP_MAX
             write(20,*) STEP_MIN,STEP_MAX,ASTEP_MIN,ASTEP_MAX
             close(20)
             open (13,file="atom_num.dat",status="old")
             write(13,*) atom_num,num
             close(13)
             open(10,file="type.dat",status="old")
             write(10,*) type,cent2
             close(10)
             open(11,file="alloy.dat",status="old")
             write(11,*) A
             close(11)
             open(34,file="parameter1.dat",status="old")
             write(34,*) kind,epsilon,zeta,p,q,rzero
             close(34)
             if(type.eq.2)then
               open(35,file="parameter2.dat",status="old")
               write(35,*) kind2,epsilon2,zeta2,p2,q2,rzero2
               open(36,file="parameter3.dat",status="old")
               write(36,*) kind3,epsilon3,zeta3,p3,q3,rzero3
               close(35)
               close(36)
             endif
             if(type.eq.3.or.type.eq.4)then
               open(37,file="protein1.dat",status="old")
               write(37,*)
     &CUTOFF,THETA,CUTOFF_ANGLE,EPS,EPS2,EPS3,EPS4,ACOIL,BCOIL
     &,CCOIL,DCOIL,ECOIL
               open(38,file="protein2.dat",status="old")
               do i=1,ndim_test
                 write(38,*) prot(i),residues_num(i)
               enddo
               close(37)
               close(38)
             endif
             if(type.eq.5)then
               open(50,file="eam_aa.dat")
               write(50,*) eam_a1,eam_epsilon1,eam_C1,eam_n1,eam_m1
               close(50)
             endif
             if (type.eq.7.or.type.eq.8)then
               open(41,file="chima1.dat",status="old")
               write(41,*) FCOIL,GCOIL,THETA1,EXC_D,dih_step
               write(41,*) a1,b1,c1,d1,e1,f1,g1,h1
               open(42,file="chima2.dat",status="old")
               do j1=1,ndim_test2
                 write(42,*) ch_st(j1),ch_end(j1)
               enddo
               close(42)
               close(41)
             endif
             write(*,*) "config. dat has been saved !"
           else if(load.eq."n")then
             write(*,*)
             write(*,*) "Load the previous configuration and quit !"
           else
             call ERROR
             goto 2
           endif
***IC=8
         else if(ic.eq.8)then
           open(30,file="data.dat")
           read(30,*) LL,ene,enef,enef2,simp,cg,gauss,d           
        read(30,*) ens_num,min,aa,bb,cc,dd,ee,ff,gg,hh,ii,par_num,samebc
     &,MAXGEN
         read(30,*) dx,bc,asw,sw_eular1,sw_eular2,eular_fac1,eular_fac2
      read(30,*) MAXIT,MUPDATE,TAKESTEP,RESETSEED,DGUESS,ASTEP,STEP,
     &GTOL,GMAX,GMAX_SW
           read(30,*) rec_sw,rec_ratio,pdb_rec,dim3
           read(30,*) recover_sw,recover_rec,recover_ratio
           read(30,*)
     &ibin,BETA0,N0,TEMPFD,EMAX,EMIN,naccept,accrat
           read(30,*) FIXTEMP_T,FIXSTEP_T,FIXBOTH_T,SCALEFAC,TEMPI,TEMPF
     &,TEMP_EX
           read(30,*) mva1,mvb1,mvc1,mvd1,mve1,mvf1,mvg1,mvh1,mvi1,mvj1
           read(30,*) mva2,mvb2,mvc2,mvd2,mve2,mvf2,mvg2,mvh2,mvi2,mvj2
           read(30,*) mva3,mvb3,mvc3,mvd3,mve3,mvf3,mvg3,mvh3,mvi3,mvj3
           read(30,*) GBIN,GFAC,FIX_GFAC,DG_SW,DG_RECOVER
           read(30,*) STEPMAX1,MCMAX1,STEPMAX2,MCMAX2,STEPMAX3,MCMAX3
           read(30,*) BH_SW1,BH_SW2,BH_SW3,DEMIN,DEMAX,ERANGE_SW
        read(30,*) rms_start,rms_end,rms_rate,rms_fac1,rms_fac2,rms_fac3
           read(30,*) sw_ensnum1,sw_parnum1,sw_same1
           read(30,*) sw_ensnum2,sw_parnum2,sw_same2
           read(30,*) sw_ensnum3,sw_parnum3,sw_same3
           read(30,*) MOVE_FAC1,MOVE_FAC2,TEMP_MIN,TEMP_MAX
           read(30,*) STEP_MIN,STEP_MAX,ASTEP_MIN,ASTEP_MAX
           close(30)

           open(18,file="atom_num.dat",status="old")
           read(18,*) atom_num,num
           close(18)
           open(25,file="type.dat",status="old")
           read(25,*) type,cent2
           close(25)
           open(26,file="alloy.dat",status="old")
           read(26,*) A
           close(26)
           open(34,file="parameter1.dat",status="old")
           read(34,*) kind,epsilon,zeta,p,q,rzero
           close(34)
           if(type.eq.2)then
             open(35,file="parameter2.dat",status="old")
             read(35,*) kind2,epsilon2,zeta2,p2,q2,rzero2
             open(36,file="parameter3.dat",status="old")
             read(36,*) kind3,epsilon3,zeta3,p3,q3,rzero3
             close(35)
             close(36)
           endif
           if(type.eq.3.or.type.eq.4)then
             open(37,file="protein1.dat",status="old")
             read(37,*)
     &CUTOFF,THETA,CUTOFF_ANGLE,EPS,EPS2,EPS3,EPS4,ACOIL,BCOIL
     &,CCOIL,DCOIL,ECOIL
             open(38,file="protein2.dat",status="old")
             ndim_test=0
55           read(38,*,end=66)
             ndim_test=ndim_test+1
             goto 55
66           rewind(38)
             do j1=1,ndim_test
               read(38,*) prot(j1),residues_num(j1)
             enddo
             close(37)
             close(38)
           endif 
           if(type.eq.5)then
               open(50,file="eam_aa.dat")
               read(50,*) eam_a1,eam_epsilon1,eam_C1,eam_n1,eam_m1
               close(50)
           endif
           if (type.eq.7.or.type.eq.8)then
             open(41,file="chima1.dat",status="old")
             read(41,*) FCOIL,GCOIL,THETA1,EXC_D,dih_step
             read(41,*) a1,b1,c1,d1,e1,f1,g1,h1
             open(42,file="chima2.dat",status="old")
             ndim_test2=0
28           read(42,*,end=27)
             ndim_test2=ndim_test2+1
             goto 28
27           rewind(42)
             do j1=1,ndim_test2
               read(42,*) ch_st(j1),ch_end(j1)
             enddo
             close(42)
             close(41)
           endif
           write(*,*) "Load the previous configuration !"
           write(*,*) "@ ",atom_num," particles"
           write(*,*) "@ seeding number= ",num,",(0 for no seeding)"
           if (gauss.eq.1)then
             write(*,*) "@ using Gaussian Random number"
           endif
           if(min.eq.1)then
             write(*,*) "@ using Simplex minimization"
           else if(min.eq.2)then
             write(*,*) "@ using Conjugate Gradient minimization"
           else if(min.eq.3)then
             write(*,*) "@ using Conjugate Gradient then Simplex"
           else if(min.eq.4)then
             write(*,*) "@ using Simplex then Conjugate Gradient"
           else
             write(*,*) "@ using LBFGS minimization,MAXIT=",MAXIT
             if(TAKESTEP.eq.1)then 
               write(*,*) "@ Take step!!" 
             else
               write(*,*) "@ Don't take step!"
             endif
             if(RESETSEED.eq.1)then
               write(*,*) "@ Reset to seeding number!!"
             else
               write(*,*) "@ Don't reset to seeding number!!"
             endif
             if(cent2.eq.1)then
               write(*,*) "@ move coords to center of mass!"
             else
               write(*,*) "@ Do not move coords to center of mass!"
             endif
           endif
           write(*,*) "@ ",ens_num," ensembles for global GA"
           write(*,*) "@ ",sw_ensnum1," ensembles for parallel GA run 1"
           write(*,*) "@ ",sw_ensnum2," ensembles for parallel GA run 2"
           write(*,*) "@ ",sw_ensnum3," ensembles for parallel GA run 3"
           if (sw_eular1.eq.1)then
           write(*,*) "@ Use eularangle tran for GGA,factor=",eular_fac1
           else
             write(*,*) "@ Do not use eularangle tran for GGA"
           endif
           if (sw_eular2.eq.1)then
           write(*,*) "@ Use eularangle tran for PGA,factor=",eular_fac2
           else
             write(*,*) "@ Do not use eularangle tran for PGA"
           endif
           write(*,*) "@ The ratio of GA:"
           write(*,*) "@ Inversion:    ",aa
           write(*,*) "@ Arithmetic:   ",bb
           write(*,*) "@ Geometic:     ",cc
           write(*,*) "@ Crossing:     ",dd
           write(*,*) "@ 2-Point:      ",ee
           write(*,*) "@ 3N Mutation:  ",ff
           write(*,*) "@ Moment op 1:  ",gg
           write(*,*) "@ Moment op 2:  ",hh
           write(*,*) "@ DIH move op 1:",ii
           write(*,*) "@ The ratio of Step move:"
           write(*,*) "@ 3N Particles Moving,1st=",mva1,",2nd=",mva2,
     &",3rd=",mva3
         write(*,*) "@ 3N Particles Mutation,1st=",mvb1,",2nd=",mvb2,
     &",3rd=",mvb3
         write(*,*) "@ N Dimension Moving,   1st=",mvc1,",2nd=",mvc2,
     &",3rd=",mvc3
         write(*,*) "@ 3N Dihedral Moving,   1st=",mve1,",2nd=",mve2,
     &",3rd=",mve3
         write(*,*) "@ RMS Transition,       1st=",mvf1,",2nd=",mvf2,
     &",3rd=",mvf3
         write(*,*) "@ RMS SCALE,            1st=",rms_fac1,",2nd=",
     &rms_fac2,",3rd=",rms_fac3
           write(*,*) "@ Inversion,          1st=",mvd1,",2nd=",mvd2,
     &",3rd=",mvd3
           write(*,*) "@ Arithmetic,         1st=",mvg1,",2nd=",mvg2,
     &",3rd=",mvg3
           write(*,*) "@ Geometic,           1st=",mvh1,",2nd=",mvh2,
     &",3rd=",mvh3
           write(*,*) "@ Crossing,           1st=",mvi1,",2nd=",mvi2,
     &",3rd=",mvi3
         write(*,*) "@ 2-Point Crossover,    1st=",mvj1,",2nd=",mvj2,
     &",3rd=",mvj3
           write(*,*) "@",par_num," parents for global GA"
           write(*,*) "@",sw_parnum1," parents for parallel GA run 1"
           write(*,*) "@",sw_parnum2," parents for parallel GA run 2" 
           write(*,*) "@",sw_parnum3," parents for parallel GA run 3" 
           write(*,*) "@ Same number=",samebc," for global GA"
           write(*,*) "@ Same number=",sw_same1," for parallel GA run 1"
           write(*,*) "@ Same number=",sw_same2," for parallel GA run 2"
           write(*,*) "@ Same number=",sw_same3," for parallel GA run 3"
           write(*,*) "@ Max gens=",MAXGEN
           write(*,*) "@ Multi canonical max number for 1st Run=",MCMAX1
        write(*,*) "@ Max MC steps for each sweep for 1st Run=",STEPMAX1
           write(*,*) "@ Multi canonical max number for 2nd Run=",MCMAX2
        write(*,*) "@ Max MC steps for each sweep for 2nd Run=",STEPMAX2
           write(*,*) "@ Multi canonical max number for 3rd Run=",MCMAX3
        write(*,*) "@ Max MC steps for each sweep for 3rd Run=",STEPMAX3
           if (d.eq.1)then
             write(*,*) "@ Using close pack condition"
           endif
           write(*,*) "@ Width of initial box is",LL
           write(*,*) "@ BC in func=",enef
           write(*,*) "@ Accuracy parameter in func=",enef2
           write(*,*) "@ The energy sort accuracy is ",ene 
           write(*,*) "@ The dX in dfunc is ",dx               
           write(*,*) "@ BC in dfunc is ",bc   
           write(*,*) "@ Tolerance of Simplex is ",simp
           write(*,*) "@ Tolerance of Conjugate Gradient is ",cg
           write(*,*) "@ Tolerance of LBFGS(GTOL)=",GTOL
           if(GMAX_SW.eq.1)then
      write(*,*) "@ Using RMS force for sloppy convergence criterion"
             write(*,*) "@ Threshold criterion=",GMAX
           else if(GMAX_SW.eq.2)then
      write(*,*) "@ The sloppy convergence criterion for RMS force=",
     &GMAX
           else if(GMAX_SW.eq.3)then
      write(*,*) "@ Use multi-convergence criterion from ",rms_start,
     &" to ",rms_end
             write(*,*) "@ The variation rate is ",rms_rate
           endif
           write(*,*) "@ MCSCRH UPDATE =",MUPDATE  
           write(*,*) "@ Takestep in bmin:Astep=",ASTEP,",Step=",STEP
           write(*,*) "@ Max Step=",STEP_MAX,",Min Step=",STEP_MIN
           write(*,*) "@ Max Astep=",ASTEP_MAX,",Min Astep=",ASTEP_MIN
           write(*,*) "@ Max TEMP=",TEMP_MAX,",Min TEMP=",TEMP_MIN
         write(*,*) "@ The adjustment of the Acceptance checking is +-",
     &MOVE_FAC1
           write(*,*) "@ Input DGUESS=",DGUESS
           if(DG_SW.eq.1.)then
             write(*,*) "@ Do not use Guided function"
           else
             write(*,*) "@ Use Guided function"
             write(*,*) "@ # of Guided Function BINS=",GBIN 
           endif
           if (FIX_GFAC.eq.1)then
             write(*,*) "@ Using temperature for Guided function=",GFAC
           else
             write(*,*) "@ The acceptance probability for Guided fun is"
     &,GFAC 
           endif
           write(*,*) "@ # of Multicanonical BH BINS=",ibin
        if(BH_SW1.eq.1.or.BH_SW1.eq.4.or.BH_SW2.eq.1.or.BH_SW3.eq.1)then
            write(*,*) "@ INITIAL BETA0=1/(",TEMPF,"-",TEMPFD,")=",BETA0
             if(ERANGE_SW.eq.0)then
               write(*,*) "@ ENERGY RANGE FROM ",EMAX," TO ",EMIN
             else if(ERANGE_SW.eq.1)then
               write(*,*) "@ Using Auto Energy Range(Half)"
               write(*,*) "@ TDE from ",DEMIN," to ",DEMAX
             else if(ERANGE_SW.eq.2)then
               write(*,*) "@ Using Multi-Canonical Annealing!"
               write(*,*) "@ TDE from ",DEMIN," to ",DEMAX
             endif
           else
             write(*,*) "@ Using MonteCarlo Boltzmann distribution"
           endif
      write(*,*)"@ Every ",naccept," generations check acceptance ratio"
           write(*,*)"@ Acceptance ratio=",accrat
           write(*,*) "@ The adjustment of the Acceptance Ratio is +-",
     &MOVE_FAC2
      write(*,*) "@ Temperature start from ",TEMPI," to ",TEMPF,
     &",multiply by ",SCALEFAC
           if(TEMP_EX.eq.1)then
             write(*,*) "@ Use temperature range for all ensembles!"
           else if(TEMP_EX.eq.2)then
             write(*,*) "@ Use temperature range for the children!"
           else
             write(*,*) "ERROR!!"
             return
           endif
           if(FIXTEMP_T.eq.1.or.FIXBOTH_T.eq.1)then
             write(*,*) "@ FIX TEMPERATURE!"
           else
             write(*,*) "@ DO NOT FIX TEMPERATURE!"
           endif
           if(FIXSTEP_T.eq.1.or.FIXBOTH_T.eq.1)then
             write(*,*) "@ FIX STEPS!"
           else
             write(*,*) "@ DO NOT FIX STEPS!"
           endif                                       
           if(rec_sw.eq.0)then
             write(*,*) "@ Normal recording"
           endif
           if(rec_sw.eq.1)then
             if (pdb_rec.eq.1)then
               write(*,*) "@ Step recording in PDB file!"
             else
               write(*,*) "@ Step recording in seperate coords file!"
             endif
             write(*,*) "@ Step recording,recording ratio=",rec_ratio
           endif
           if(rec_sw.eq.2)then
             write(*,*) "@ Record nothing!"
           endif
           if(recover_rec.eq.1)then
             write(*,*) "@ Start recover,recover ratio=",recover_ratio
           endif
           if (type.eq.1)then
             write(*,*) "@    e0        c0        p    q     r0"
             write(*,*) "@ ",epsilon,",",zeta,",",p,",",q,",",rzero  
             write(*,*) "@ For ",kind
           endif
           if (type.eq.2)then
             write(*,*) "@    e0        c0        p    q     r0"
             write(*,*) "@ ",epsilon,",",zeta,",",p,",",q,",",rzero
             write(*,*) "@ For ",A," ",kind   
             write(*,*) "@    e0        c0        p    q     r0"
        write(*,*) "@ ",epsilon2,",",zeta2,",",p2,",",q2,",",rzero2
             write(*,*) "@ For ",atom_num-A," ",kind2           
             write(*,*) "@    e0        c0        p    q     r0"
        write(*,*) "@ ",epsilon3,",",zeta3,",",p3,",",q3,",",rzero3
             write(*,*) "@ For ",kind3    
           endif
           if (type.eq.3)then
             write(*,*) "@ For Beta-Sheets"
             write(*,*) "@ ",(prot(j1),residues_num(j1),j1=1,ndim_test)
           endif
           if (type.eq.4)then
             write(*,*) "@ For Alpha-Helices"
             write(*,*) "@ ",(prot(j1),residues_num(j1),j1=1,ndim_test)
           endif
           if (type.eq.5)then
              write(*,*) "@  eam_a  eam_e  eam_C  eam_n  eam_m"
              write(*,*) "@ ",eam_a1,",",eam_epsilon1,",",eam_C1,","
     &,eam_n1,",",eam_m1
              write(*,*) "@ For ",kind
           endif
           if (type.eq.7)then
             write(*,*) "@ For Chen and Imamura Folding N"
             write(*,*) "@ The Neutral bond sets:(Start-end)"
          write(*,*) "@ ",("(",ch_st(j1),ch_end(j1),")",j1=1,ndim_test2)
      write(*,*) "@ Lo=",FCOIL,",Epsilon=",GCOIL,",ANGLE=",THETA1,
     &",Exclusive diameter=",EXC_D
           endif
           if (type.eq.8)then
             write(*,*) "@ For Chen and Imamura Folding 3N"
             write(*,*) "@ The Neutral bond sets:(Start-end)"
          write(*,*) "@ ",("(",ch_st(j1),ch_end(j1),")",j1=1,ndim_test2)
             write(*,*) "@ Lo=",FCOIL,",Epsilon=",GCOIL,",ANGLE=",THETA1
     &,",Exclusive diameter=",EXC_D
           endif
           if (asw.eq.1)then
             write(*,*) "@ Using analytic gradient"
           else
             write(*,*) "@ Using numerical gradient"
           endif
         else
           call ERROR
           goto 2
         endif
       else
         call ERROR
         goto 1
       endif          
       write(*,*) "--------------------------------------------------"
       write(*,*) " To run mpi: mpirun -np ## a.out &                " 
       write(*,*) " which                  ## is node number         "
       write(*,*) "             MPI-PTMBHPGA V2.4                    "
       write(*,*) "             update:2003/6/21                     "
       write(*,*) "               By Po-Jen Hsu                      "
       write(*,*) "--------------------------------------------------" 
       stop
       end
 
       subroutine ERROR
       write(*,*)
       write(*,*) "WARING!!!!!!!!!!!!!!!!!!!!!!!!!"
       write(*,*) "Can not recognize the input parameter!"
       write(*,*) "Please select again!"
       write(*,*)
       return
       end 
