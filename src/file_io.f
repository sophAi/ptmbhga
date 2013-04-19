      subroutine file_io(IBIN,ens_num)
      implicit none
      include "mpif.h"
      integer nproc,myid
      integer ens_num,par_num,MAXGEN,samebc,num,asw,dim3,d,G_SW
      integer BH_SW1,BH_SW2,BH_SW3,ERANGE_SW,sw_eular1,sw_eular2,PRE_BH
      integer sw_ensnum1,sw_parnum1,sw_same1
      integer sw_ensnum2,sw_parnum2,sw_same2
      integer sw_ensnum3,sw_parnum3,sw_same3
      integer MUPDATE,MAXIT,FALSEIT,type,rec_sw,recover_rec
      integer TAKESTEP,RESETSEED,recover_sw,min,pdb_rec,DG_SW,temp_ex
      integer IBIN,N0,MCMAX,STEPMAX,J1,naccept,FIX_GFAC,GBIN,
     &rms_tran_sw
      integer FIXBOTH_T,FIXTEMP_T,FIXSTEP_T,GAUSS_T,GMAX_SW,DG_RECOVER
      integer STEPMAX1,MCMAX1,STEPMAX2,MCMAX2,STEPMAX3,MCMAX3
      real*8 BETA0,TEMPFD,DBETA,EMAX,EMIN,EMINGL,STEP_TEMP,ASTEP_TEMP
      real*8 ESYS,SCALEFAC,DG,GFAC,GMAX,eular_fac1,eular_fac2
      real*8 simp,cg,bc,recover_ratio,TEMP,TEMPI,TEMPF,TEMP_TEMP
      real*8 aa,bb,cc,dd,ee,ff,gg,hh,ii,LL,ene,enef,enef2,dx
      real*8 aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii
      real*8 DGUESS,ASTEP,STEP,GTOL,rec_ratio,accrat,DEMIN,DEMAX
      real*8 mva1,mvb1,mvc1,mvd1,mve1,mvf1,mvg1,mvh1,mvi1,mvj1,
     &rms_fac1
      real*8 mva2,mvb2,mvc2,mvd2,mve2,mvf2,mvg2,mvh2,mvi2,mvj2,
     &rms_fac2
      real*8 mva3,mvb3,mvc3,mvd3,mve3,mvf3,mvg3,mvh3,mvi3,mvj3,
     &rms_fac3
      real*8 mv1_1,mv2_1,mv3_1,mv4_1,mv5_1,mv6_1,mv7_1,mv8_1,mv9_1,
     &mv10_1,rzero_max
      real*8 mv1_2,mv2_2,mv3_2,mv4_2,mv5_2,mv6_2,mv7_2,mv8_2,mv9_2,
     &mv10_2
      real*8 mv1_3,mv2_3,mv3_3,mv4_3,mv5_3,mv6_3,mv7_3,mv8_3,mv9_3,
     &mv10_3
      real*8 rms_dx1,rms_dx2,rms_start,rms_end,rms_rate
      real*8 MOVE_FAC1,MOVE_FAC2,TEMP_MIN,TEMP_MAX
      real*8 STEP_MIN,STEP_MAX,ASTEP_MIN,ASTEP_MAX
      logical FIXTEMP,FIXSTEP,GAUSS
      character name*6
      common nproc,myid
      common /name1/name
      common /seeding/ num
      common /same1/ samebc,MAXGEN
      common /accuracy1/ LL
      common /accuracy2/ ene
      common /accuracy3/ enef,enef2
      common /accuracy4/ dx,bc
      common /accuracy5/ simp
      common /accuracy6/ cg
      common /anatsw/ asw
      common /bmin1/ DGUESS,ASTEP,STEP,GTOL
      common /bmin2/ MAXIT,MUPDATE,TAKESTEP,RESETSEED,FALSEIT
      common /bmin3/ GMAX,GMAX_SW
      common /engtype/ type
      common /rec_sw1/ rec_sw
      common /rec_sw2/ rec_ratio
      common /recover1/ recover_sw,recover_rec
      common /recover2/ recover_ratio
      common /mpigen1/ aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii
      common /mpigen2/ par_num
      common /pdbrec/ pdb_rec,dim3
*  <<Multicanonical MC parameters>>
      common /MULTIC1/ EMAX,EMIN
      common /MULTIC2/ ERANGE_SW,PRE_BH
      common /UPPER/ BETA0 ,N0 ,STEPMAX,MCMAX
      common /GLOBAL/ EMINGL
      common /ANNEAL/ TEMPFD
      common /STEP1/ FIXTEMP,FIXSTEP
      common /TEMP1/ SCALEFAC,TEMP
      common /TEMP2/ TEMPI,TEMPF
      common /TEMP3/ temp_ex
      common /accept1/ naccept
      common /accept2/ accrat
      common /accept3/ MOVE_FAC1,MOVE_FAC2
      common /accept4/ STEP_MIN,STEP_MAX,ASTEP_MIN,ASTEP_MAX
      common /accept5/ TEMP_MIN,TEMP_MAX
      common /accept6/ STEP_TEMP,ASTEP_TEMP,TEMP_TEMP
      common /move1/ mv1_1,mv2_1,mv3_1,mv4_1,mv5_1,mv6_1,mv7_1,mv8_1,
     &mv9_1,mv10_1,rms_fac1
      common /move2/ mv1_2,mv2_2,mv3_2,mv4_2,mv5_2,mv6_2,mv7_2,mv8_2,
     &mv9_2,mv10_2,rms_fac2
      common /move3/ mv1_3,mv2_3,mv3_3,mv4_3,mv5_3,mv6_3,mv7_3,mv8_3,
     &mv9_3,mv10_3,rms_fac3
      common /gauss1/ GAUSS
      common /alloy2/ rzero_max
      common /guide1/ DG,DG_SW
      common /guide2/ GFAC,FIX_GFAC
      common /guide3/ DG_RECOVER
      common /min1/ min
      common /gabh1/ STEPMAX1,MCMAX1,STEPMAX2,MCMAX2,STEPMAX3,MCMAX3
      common /BH2/ BH_SW1,BH_SW2,BH_SW3
      common /BH3/ DEMIN,DEMAX
      common /RMS2/ rms_start,rms_end,rms_dx1,rms_dx2,rms_rate
      common /RMS4/ rms_tran_sw
      common /SW_GA1/ sw_ensnum1,sw_parnum1,sw_same1,sw_ensnum2,
     &sw_parnum2,sw_same2,sw_ensnum3,sw_parnum3,sw_same3
      common /eularangle1/ sw_eular1
      common /eularangle2/ sw_eular2
      common /eularangle3/ eular_fac1
      common /eularangle4/ eular_fac2
      common /WEIGHT1/ G_SW
      rms_tran_sw=0
      FALSEIT=0
      open(10,file="data.dat",status="old")
      read(10,*) LL,ene,enef,enef2,simp,cg,GAUSS_T,d      !accuracy
      read(10,*) ens_num,min,aa,bb,cc,dd,ee,ff,gg,hh,ii,par_num,samebc, 
     &MAXGEN                                  !config
      read(10,*) dx,bc,asw,sw_eular1,sw_eular2,eular_fac1,eular_fac2
      read(10,*) MAXIT,MUPDATE,TAKESTEP,RESETSEED,DGUESS,ASTEP,
     &STEP,GTOL,GMAX,GMAX_SW                  !bmin
      read(10,*) rec_sw,rec_ratio,pdb_rec,dim3   !rec
      read(10,*) recover_sw,recover_rec,recover_ratio    !recover
      read(10,*)
     & IBIN,BETA0,N0,TEMPFD,EMAX,EMIN,naccept,accrat       !anneal
      read(10,*) FIXTEMP_T,FIXSTEP_T,FIXBOTH_T,SCALEFAC,TEMPI,TEMPF,
     &temp_ex                          !temperature
      read(10,*) mva1,mvb1,mvc1,mvd1,mve1,mvf1,mvg1,mvh1,mvi1,mvj1
      read(10,*) mva2,mvb2,mvc2,mvd2,mve2,mvf2,mvg2,mvh2,mvi2,mvj2
      read(10,*) mva3,mvb3,mvc3,mvd3,mve3,mvf3,mvg3,mvh3,mvi3,mvj3
      read(10,*) GBIN,GFAC,FIX_GFAC,DG_SW,DG_RECOVER      
      read(10,*) STEPMAX1,MCMAX1,STEPMAX2,MCMAX2,STEPMAX3,MCMAX3
      read(10,*) BH_SW1,BH_SW2,BH_SW3,DEMIN,DEMAX,ERANGE_SW,PRE_BH,G_SW
      read(10,*) rms_start,rms_end,rms_rate,rms_fac1,rms_fac2,rms_fac3
      read(10,*) sw_ensnum1,sw_parnum1,sw_same1
      read(10,*) sw_ensnum2,sw_parnum2,sw_same2
      read(10,*) sw_ensnum3,sw_parnum3,sw_same3
      read(10,*) MOVE_FAC1,MOVE_FAC2,TEMP_MIN,TEMP_MAX
      read(10,*) STEP_MIN,STEP_MAX,ASTEP_MIN,ASTEP_MAX
      close(10)
      rms_dx1=(rms_start-rms_end)/dble(ens_num)
      rms_dx2=(rms_start-rms_end)/dble(ens_num-par_num)
      STEPMAX=STEPMAX1
      MCMAX=MCMAX1
      TEMP=TEMPI
      TEMP_TEMP=TEMP
      STEP_TEMP=STEP
      ASTEP_TEMP=ASTEP 
****************************************************
C   THIS IS VERY IMPORTANT FOR CLOSE PACK CONDITION"
C   Input your boundary condition here!"
      if (type.eq.1.or.type.eq.2)then
        enef=enef*(rzero_max**2)
        LL=enef
      endif
****************************************************
      aaa=aa/(aa+bb+cc+dd+ee+ff+gg+hh+ii)
      bbb=bb/(aa+bb+cc+dd+ee+ff+gg+hh+ii)
      ccc=cc/(aa+bb+cc+dd+ee+ff+gg+hh+ii)
      ddd=dd/(aa+bb+cc+dd+ee+ff+gg+hh+ii)
      eee=ee/(aa+bb+cc+dd+ee+ff+gg+hh+ii)
      fff=ff/(aa+bb+cc+dd+ee+ff+gg+hh+ii)
      ggg=gg/(aa+bb+cc+dd+ee+ff+gg+hh+ii)
      hhh=hh/(aa+bb+cc+dd+ee+ff+gg+hh+ii)
      iii=ii/(aa+bb+cc+dd+ee+ff+gg+hh+ii)
      mv1_1=mva1/(mva1+mvb1+mvc1+mvd1+mve1+mvf1+mvg1+mvh1+mvi1+mvj1)
      mv2_1=mvb1/(mva1+mvb1+mvc1+mvd1+mve1+mvf1+mvg1+mvh1+mvi1+mvj1)
      mv3_1=mvc1/(mva1+mvb1+mvc1+mvd1+mve1+mvf1+mvg1+mvh1+mvi1+mvj1)
      mv4_1=mvd1/(mva1+mvb1+mvc1+mvd1+mve1+mvf1+mvg1+mvh1+mvi1+mvj1)
      mv5_1=mve1/(mva1+mvb1+mvc1+mvd1+mve1+mvf1+mvg1+mvh1+mvi1+mvj1)
      mv6_1=mvf1/(mva1+mvb1+mvc1+mvd1+mve1+mvf1+mvg1+mvh1+mvi1+mvj1)   
      mv7_1=mvg1/(mva1+mvb1+mvc1+mvd1+mve1+mvf1+mvg1+mvh1+mvi1+mvj1)
      mv8_1=mvh1/(mva1+mvb1+mvc1+mvd1+mve1+mvf1+mvg1+mvh1+mvi1+mvj1)
      mv9_1=mvi1/(mva1+mvb1+mvc1+mvd1+mve1+mvf1+mvg1+mvh1+mvi1+mvj1)
      mv10_1=mvj1/(mva1+mvb1+mvc1+mvd1+mve1+mvf1+mvg1+mvh1+mvi1+mvj1)
      mv1_2=mva2/(mva2+mvb2+mvc2+mvd2+mve2+mvf2+mvg2+mvh2+mvi2+mvj2)
      mv2_2=mvb2/(mva2+mvb2+mvc2+mvd2+mve2+mvf2+mvg2+mvh2+mvi2+mvj2)
      mv3_2=mvc2/(mva2+mvb2+mvc2+mvd2+mve2+mvf2+mvg2+mvh2+mvi2+mvj2)
      mv4_2=mvd2/(mva2+mvb2+mvc2+mvd2+mve2+mvf2+mvg2+mvh2+mvi2+mvj2)
      mv5_2=mve2/(mva2+mvb2+mvc2+mvd2+mve2+mvf2+mvg2+mvh2+mvi2+mvj2)
      mv6_2=mvf2/(mva2+mvb2+mvc2+mvd2+mve2+mvf2+mvg2+mvh2+mvi2+mvj2)
      mv7_2=mvg2/(mva2+mvb2+mvc2+mvd2+mve2+mvf2+mvg2+mvh2+mvi2+mvj2)
      mv8_2=mvh2/(mva2+mvb2+mvc2+mvd2+mve2+mvf2+mvg2+mvh2+mvi2+mvj2)
      mv9_2=mvi2/(mva2+mvb2+mvc2+mvd2+mve2+mvf2+mvg2+mvh2+mvi2+mvj2)
      mv10_2=mvj2/(mva2+mvb2+mvc2+mvd2+mve2+mvf2+mvg2+mvh2+mvi2+mvj2)
      mv1_3=mva3/(mva3+mvb3+mvc3+mvd3+mve3+mvf3+mvg3+mvh3+mvi3+mvj3)
      mv2_3=mvb3/(mva3+mvb3+mvc3+mvd3+mve3+mvf3+mvg3+mvh3+mvi3+mvj3)
      mv3_3=mvc3/(mva3+mvb3+mvc3+mvd3+mve3+mvf3+mvg3+mvh3+mvi3+mvj3)
      mv4_3=mvd3/(mva3+mvb3+mvc3+mvd3+mve3+mvf3+mvg3+mvh3+mvi3+mvj3)
      mv5_3=mve3/(mva3+mvb3+mvc3+mvd3+mve3+mvf3+mvg3+mvh3+mvi3+mvj3)
      mv6_3=mvf3/(mva3+mvb3+mvc3+mvd3+mve3+mvf3+mvg3+mvh3+mvi3+mvj3)
      mv7_3=mvg3/(mva3+mvb3+mvc3+mvd3+mve3+mvf3+mvg3+mvh3+mvi3+mvj3)
      mv8_3=mvh3/(mva3+mvb3+mvc3+mvd3+mve3+mvf3+mvg3+mvh3+mvi3+mvj3)
      mv9_3=mvi3/(mva3+mvb3+mvc3+mvd3+mve3+mvf3+mvg3+mvh3+mvi3+mvj3)
      mv10_3=mvj3/(mva3+mvb3+mvc3+mvd3+mve3+mvf3+mvg3+mvh3+mvi3+mvj3)
      IF(FIXTEMP_T.EQ.1)THEN
        FIXTEMP=.TRUE.
      ELSE
        FIXTEMP=.FALSE.
      ENDIF
      IF(FIXSTEP_T.EQ.1)THEN
        FIXSTEP=.TRUE.
      ELSE
        FIXSTEP=.FALSE.
      ENDIF
C      IF(FIXBOTH_T.EQ.1)THEN
C        FIXBOTH=.TRUE.
C        FIXTEMP=.TRUE.
C        FIXSTEP=.TRUE.
C      ELSE
C        FIXBOTH=.FALSE.
C      ENDIF
      IF(GAUSS_T.EQ.1)THEN
        GAUSS=.TRUE.
      ELSE
        GAUSS=.FALSE.
      ENDIF
      EMINGL=999999.D0
C      DE=(EMAX-EMIN)/DBLE(IBIN)
      if (FIX_GFAC.eq.1) GFAC=TEMP
      if (myid.eq.0)then
        write(*,*) "@ You use ",nproc," node(s)"
        write(*,*) "@ seeding number= ",num,",(0 for no seeding)"
        write(*,*) "@ ",ens_num," ensembles for global GA"
        write(*,*) "@ ",par_num," parents ,",ens_num-par_num," childs"
        write(*,*) "@ Same number=",samebc," for global GA"
        if (BH_SW1.eq.3.or.BH_SW2.eq.3.or.BH_SW3.eq.3)then
        write(*,*)"@ ",sw_ensnum1," ensembles for parallel GA 1st Run"
        write(*,*)"@ ",sw_ensnum2," ensembles for parallel GA 2nd Run"
        write(*,*)"@ ",sw_ensnum3," ensembles for parallel GA 3rd Run"
        write(*,*)"@ ",sw_parnum1," parents ,",sw_ensnum1-sw_parnum1,
     &" childs for 1st Run "
        write(*,*)"@ ",sw_parnum2," parents ,",sw_ensnum2-sw_parnum2,
     &" childs for 2nd Run "
        write(*,*)"@ ",sw_parnum3," parents ,",sw_ensnum3-sw_parnum3,
     &" childs for 3rd Run "
      write(*,*)"@ Same number=",sw_same1," for parallel GA for 1st Run"
      write(*,*)"@ Same number=",sw_same2," for parallel GA for 2nd Run"
      write(*,*)"@ Same number=",sw_same3," for parallel GA for 3rd Run"
        endif
        write(*,*) "@",MAXGEN," Max generations BH"
        write(*,*) "@ 1st Run,",STEPMAX1," Steps,",MCMAX1," MMC"
        write(*,*) "@ 2nd Run,",STEPMAX2," Steps,",MCMAX2," MMC"
        write(*,*) "@ 3rd Run,",STEPMAX3," Steps,",MCMAX3," MMC"
        write(*,*) "@ MC Select:1=MMC,2=SMC,3=PGA,4 or 5=MSA(Annealing)"
        write(*,*) "@ 1st=",BH_SW1,",2nd=",BH_SW2,",3rd=",BH_SW3
       if(BH_SW1.eq.1.or.BH_SW1.eq.4.or.BH_SW2.eq.1.or.BH_SW3.eq.1) then
          write(*,*) "@ INITIAL BETA0=1/(",TEMPF,"-",TEMPFD,")=",BETA0
          if(ERANGE_SW.eq.2)then
            write(*,*) "@ Using Multi-Canonical Annealing BH"
            write(*,*) "@ TDE form ",DEMIN," to ",DEMAX
          else if(ERANGE_SW.eq.1)then
            write(*,*) "@ USING AUTO ENERGY RANGE(HALF)"
            write(*,*) "@ TDE from ",DEMIN," to ",DEMAX
          else if(ERANGE_SW.eq.0)then
            if(G_SW.eq.1)then
              write(*,*) "@ Using Berg weighting Multi-Canonical BH"
            else
              write(*,*) "@ Using usual Multi-Canonical BH"
            endif
            write(*,*) "@ ENERGY RANGE FROM ",EMAX," TO ",EMIN
          endif
        else
          write(*,*) "@ Using Single-Canonical BH or Parallel GA"
        endif
        if (min.eq.1)then
          write(*,*) "@ using Simplex minimization"
        else if(min.eq.2)then
          write(*,*) "@ using Conjugate Gradient minimization"
        else if(min.eq.3)then  
          write(*,*) "@ using Conjugate Gradient then Simplex"
        else if(min.eq.4)then
          write(*,*) "@ using Simplex then Conjugate Gradient"
        else
          write(*,*) "@ using LBFGS minimization,MAXIT=",MAXGEN
          if(GMAX_SW.eq.1)then
            write(*,*) "@ Use RMS force for Sloppy convergence"
          else if(GMAX_SW.eq.2)then
            write(*,*) "@ sloppy convergence criterion=",GMAX
          else
            write(*,*) "@ Use variate convergence criterion from,",
     &rms_start," to ",rms_end
          endif
          if (TAKESTEP.eq.1)then
            write(*,*) "@ Take step!!"
          else
            write(*,*) "@ Don't take step!!"
          endif
          if (RESETSEED.eq.1)then
            write(*,*) "@ Reset to seeding number!!"
          else
            write(*,*) "@ Don't reset to seeding number!!"
          endif
        endif
        if(FIXSTEP)then
          write(*,*) "@ Fix step and astep"
        endif
        if(FIXTEMP)then
          write(*,*) "@ Fix temperature"
        endif
        if (d.eq.1)then
          write(*,*) "@ Using close pack initial condition:",LL
        else
          write(*,*) "@ Do not use close pack initial condition:",LL
        endif
        if (temp_ex.eq.1)then
          write(*,*) "@ Temp=,",TEMPI," to ",TEMPF," for all ensembles"
        else if(temp_ex.eq.2)then
          write(*,*) "@ Temp=,",TEMPI," to ",TEMPF," for all children"
        endif
        write(*,*) "@ Radius=",enef," ,quenching rate=",SCALEFAC
      write(*,*) "@ The ratio of Genetic Operator"
        write(*,*) "@ Inversion:  ",aa,",Arithmetic: ",bb
        write(*,*) "@ Geometic:   ",cc,".Crossing:   ",dd
        write(*,*) "@ 2-Point:    ",ee,",3N Mutation:",ff
        write(*,*) "@ Moment op 1:",gg,",Moment op 2:",hh
        write(*,*) "@ DIH MOVE-1: ",ii
        write(*,*) "@ The ratio of the Random Move:(1st/2nd/3rd)"
        write(*,*) "@ 3N Moving:",mva1,mva2,mva3," 3N Mutation:    "
     &,mvb1,mvb2,mvb3
        write(*,*) "@ N Moving:  ",mvc1,mvc2,mvc3," 3N Dih  Moving:"
     &,mve1,mve2,mve3
        write(*,*) "@ RMS Moving:",mvf1,mvf2,mvf3," RMS Scale:     "
     &,rms_fac1,rms_fac2,rms_fac3
        write(*,*) "@ Inversion: ",mvd1,mvd2,mvd3," Arithmetic:    "
     &,mvg1,mvg2,mvg3
        write(*,*) "@ Geometic:  ",mvh1,mvh2,mvh3," Crossing:      "
     &,mvi1,mvi2,mvi3
        write(*,*) "@ 2-Point Crossover:",mvj1,mvj2,mvj3
        if (rec_sw.eq.0)then           
          write(*,*) "@ Normal recording!"
        endif
        if (rec_sw.eq.1)then
          if (pdb_rec.eq.1)then
      write(*,*) "@ Step record in PDB Animation,rec ratio=",rec_ratio
          else
      write(*,*) "@ Step record in seperate files,rec ratio=",rec_ratio
          endif
        endif
        if (rec_sw.eq.2)then
          write(*,*) "@ Record nothing!"
        endif
        if (sw_eular1.eq.1)then
        write(*,*)"@ Turn on eularangle tran for GGA,factor=",eular_fac1
        endif
        if (sw_eular2.eq.1)then
        write(*,*)"@ Turn on rularangle tran for PGA,factor=",eular_fac2
        endif
C       WRITE INFORMATION FILE 
        open(92,file=name//".inf",status="unknown")
        write(92,*) "<<<Config information v2.8>>>"
        write(92,*) "@ Last update:2004/9/10"
        write(92,*) "@ For",name
        write(92,*) "@ seeding number= ",num,",(0 for no seeding)"
        if (min.eq.1)then
          write(92,*) "@ using Simplex minimization"
        else if(min.eq.2)then  
          write(92,*) "@ using Conjugate Gradient minimization"
        else if(min.eq.3)then  
          write(92,*) "@ using Conjugate Gradient then Simplex"
        else if(min.eq.4)then
          write(92,*) "@ using Simplex then Conjugate Gradient"
        else
          write(92,*) "@ using LBFGS minimization,MAXIT",MAXIT
          if(TAKESTEP.eq.1)then
            write(92,*) "@ Take step!!"
          else
            write(92,*) "@ Don't take step!!"
          endif
          if(RESETSEED.eq.1)then
            write(92,*) "@ Reset to seeding number!!"
          else
            write(92,*) "@ Don't reset to seeding number!!"
          endif
        endif
        if(FIXSTEP)then
          write(92,*) "@ Fix step and astep"
        endif
        if(FIXTEMP)then
          write(92,*) "@ Fix temperature"
        endif
        write(92,*) "@ ",ens_num," ensembles for global GA"
        write(92,*) "@ ",sw_ensnum1," ensembles for parallel GA 1st Run"
        write(92,*) "@ ",sw_ensnum2," ensembles for parallel GA 2nd Run"
        write(92,*) "@ ",sw_ensnum3," ensembles for parallel GA 3rd Run"
        if (sw_eular1.eq.1)then
       write(92,*)"@ Turn on eularangle tran for GGA,factor=",eular_fac1
        else
          write(92,*) "@ Turn off eularangle transformation for GGA"
        endif
        if (sw_eular2.eq.1)then
       write(92,*)"@ Turn on eularangle tran for PGA,factor=",eular_fac2
        else
          write(92,*) "@ Turn off eularangle transformation for PGA"
        endif
        write(92,*) "@ The ratio of GA:"
        write(92,*) "@ Inversion:  ",aa
        write(92,*) "@ Arithmetic: ",bb
        write(92,*) "@ Geometic:   ",cc
        write(92,*) "@ Crossing:   ",dd
        write(92,*) "@ 2-Point:    ",ee
        write(92,*) "@ 3N Mutation:",ff
        write(92,*) "@ Moment op 1:",gg
        write(92,*) "@ Moment op 2:",hh
        write(92,*) "@ DIH MOVE-1: ",ii
        write(92,*) "@ The ratio of moving:    1st 2nd 3rd"
        write(92,*) "@ 3N Particles Moving:  ",mva1,mva2,mva3
        write(92,*) "@ 3N Particles Mutation:",mvb1,mvb2,mvb3
        write(92,*) "@ N Dimension Moving:   ",mvc1,mvc2,mvc3
        write(92,*) "@ 3N Dihedral Moving:   ",mve1,mve2,mve3
        write(92,*) "@ RMS transition:       ",mvf1,mvf2,mvf3
        write(92,*)"@ RMS transition Scale: ",rms_fac1,rms_fac2,rms_fac3
        write(92,*) "@ Inversion:            ",mvd1,mvd2,mvd3
        write(92,*) "@ Arithmetic:           ",mvg1,mvg2,mvg3
        write(92,*) "@ Geometic:             ",mvh1,mvh2,mvh3
        write(92,*) "@ Crossing:             ",mvi1,mvi2,mvi3
        write(92,*) "@ 2-Point Crossover:    ",mvj1,mvj2,mvj3
        write(92,*) "@",par_num," parents for global GA"
        write(92,*) "@",sw_parnum1," parents for parallel GA 1st Run"
        write(92,*) "@",sw_parnum2," parents for parallel GA 2nd Run"
        write(92,*) "@",sw_parnum3," parents for parallel GA 3rd Run"
        write(92,*) "@ Same number=",samebc," for global GA"
        write(92,*) "@ Same number=",sw_same1," for parallel GA 1st Run"
        write(92,*) "@ Same number=",sw_same2," for parallel GA 2nd Run"
        write(92,*) "@ Same number=",sw_same3," for parallel GA 3rd Run"
        write(92,*) "@ Max gens=",MAXGEN
        write(92,*)"@ MC Select:1=MMC,2=SMC,3=PGA,4 or 5=MSA(Annealing)"
        write(92,*) "@ 1st=",BH_SW1,",2nd=",BH_SW2,",3rd=",BH_SW3
        write(92,*) "@ Multi canonical max number for 1st Run=",MCMAX1
       write(92,*) "@ Max MC steps for each sweep for 1st Run=",STEPMAX1
        write(92,*) "@ Multi canonical max number for 2nd Run=",MCMAX2
       write(92,*) "@ Max MC steps for each sweep for 2nd Run=",STEPMAX2
        write(92,*) "@ Multi canonical max number for 3rd Run=",MCMAX3
       write(92,*) "@ Max MC steps for each sweep for 3rd Run=",STEPMAX3
        if (d.eq.1)then
          write(92,*) "@ Using close pack conditon"
        endif
        write(92,*) "@ Width of initial box is",LL
        write(92,*) "@ BC in func=",enef
        write(92,*) "@ Accuracy parameter in func=",enef2
        write(92,*) "@ The energy sort accuracy is ",ene
        write(92,*) "@ The dX in dfunc is ",dx
        write(92,*) "@ BC in dfunc is ",bc
        write(92,*) "@ Tolerance of Simplex is ",simp
        write(92,*) "@ Tolerance of Conjugate Gradient is ",cg
        write(92,*) "@ Tolerance of LBFGS(GTOL)=",GTOL
        write(92,*) "@ MCSCRH UPDATE =",MUPDATE
        write(92,*) "@ Takestep in bmin:Astep=",ASTEP,",Step=",STEP
        write(92,*) "@ Max Step=",STEP_MAX,",Min Step=",STEP_MIN
        write(92,*) "@ Max Astep=",ASTEP_MAX,",Min Astep=",ASTEP_MIN
        write(92,*) "@ Max TEMP=",TEMP_MAX,",Min TEMP=",TEMP_MIN
        write(92,*) "@ The adjustment of the Acceptance checking is +-",
     &MOVE_FAC1
        write(92,*) "@ Input DGUESS=",DGUESS
        if(GMAX_SW.eq.1)then
      write(92,*) "@ Using RMS force for sloppy convergence criterion"
           write(92,*) "@ Threshold criterion=",GMAX
         else if(GMAX_SW.eq.2)then
      write(92,*) "@ The sloppy convergence criterion for RMS force=",
     &GMAX
         else if(GMAX_SW.eq.3)then
      write(92,*) "@ Use multi-convergence criterion from ",rms_start,
     &" to ",rms_end
           write(92,*) "@ The variation rate is ",rms_rate
        endif
        if(DG_SW.eq.1)then
          write(92,*) "@ Do not use Guided function"
        else
          write(92,*) "@ Use Guided function"
        endif
        if (FIX_GFAC.eq.1)then
          write(92,*) "@ Using temperature for Guided function=",GFAC
        else
          write(92,*) "@ The acceptance probability for Guided fun is"
     &,GFAC
        endif
        write(92,*) "@ # of Multicanonical BH BINS=",ibin
        if(BH_SW1.eq.1.or.BH_SW1.eq.4.or.BH_SW2.eq.1.or.
     &BH_SW3.eq.1) then
          write(92,*) "@ INITIAL BETA0=1/(",TEMPF,"-",
     &TEMPFD,")=",BETA0
          if(ERANGE_SW.eq.2)then
            write(92,*) "@ Using Multi-Canonical Annealing BH"
            write(92,*) "@ TDE from ",DEMIN," to ",DEMAX
            write(92,*) "@ DE from ",DEMIN/dble(ibin)," to ",
     &DEMAX/dble(ibin)
          else if(ERANGE_SW.eq.1)then
            if(G_SW.eq.1)then
              write(92,*) "@ Using Berg Multi-Canonical BH"
            else
               write(92,*) "@ Using usual Multi-Canonical BH"
            endif
            write(92,*) "@ Using auto energy range(half)"
            write(92,*) "@ TDE from ",DEMIN, " to ",DEMAX
            write(92,*) "@ DE from ",DEMIN/dble(ibin)," to ",
     &DEMAX/dble(ibin)
          else if(ERANGE_SW.eq.0)then
            write(92,*) "@ Using usual Multi-Canonical BH" 
            write(92,*) "@ ENERGY RANGE FROM ",EMAX," TO ",EMIN
          endif
        else
          write(92,*) "@ Using MonteCarlo Boltzmann distribution"
        endif
      write(92,*)"@ Every ",naccept,
     &" generations check acceptance ratio"
        write(92,*)"@ Acceptance ratio=",accrat
        write(92,*) "@ The adjustment of the Acceptance Ratio is +-",
     &MOVE_FAC2
        write(92,*) "@ Temperature start FROM ",TEMPI," to ",TEMPF,
     &",multiply by ",SCALEFAC
        if (temp_ex.eq.1)then
          write(92,*) "@ Use temperature range for all ensembles!"
        else if(temp_ex.eq.2)then
          write(92,*) "@ Use temperature range for the children!"
        else
          write(*,*) "ERROR!!"
          return
        endif
        if(FIXTEMP_T.eq.1.or.FIXBOTH_T.eq.1)then
          write(92,*) "@ FIX TEMPERATURE!"
        else
          write(92,*) "@ DO NOT FIX TEMPERATURE!"
        endif
        if(FIXSTEP_T.eq.1.or.FIXBOTH_T.eq.1)then
          write(92,*) "@ FIX STEPS!"
        else
          write(92,*) "@ DO NOT FIX STEPS!"
        endif
        if (asw.eq.1)then
           write(92,*) "@ Using analytic gradient"
         else 
           write(92,*) "@ Using numerical gradient"
         endif

        if(rec_sw.eq.0)then
          write(92,*) "@ Normal recording"
        endif
        if(rec_sw.eq.1)then
          if (pdb_rec.eq.1)then
            write(92,*) "@ Step recording in PDB file!"
          else
            write(92,*) "@ Step recording in seperate coords file!" 
          endif  
          write(92,*) "@ Step recording,recording ratio=",rec_ratio
        endif
        if(rec_sw.eq.2)then
          write(92,*) "@ Record nothing!"
        endif
        if(recover_rec.eq.1) then
          write(92,*) "@ Start backup,backup ratio=",recover_ratio
        endif
******INPUT THE POTENTIAL DETIALS********************
        close(92)
      endif
      return
      end
