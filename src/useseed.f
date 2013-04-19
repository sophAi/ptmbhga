       subroutine useseed(node,ndim,ibin,ens_num,seed,
     &outx,upot)
       implicit none
       include "mpif.h"
C===============================
       integer myid,nproc,node,ndim,ibin,ens_num,par_num,seed,min,MAXGEN
       integer same,samebc,almin,num,ndim2,dim3,type,MULMAX,STEPMAX,N0
       real*8 same_temp(1),almin_temp(1),np_temp(ens_num)
       real*8 mv1,mv2,mv3,mv4,mv5,mv6,mv7,mv8,mv9,mv10,rms_fac
       real*8 mv1_1,mv2_1,mv3_1,mv4_1,mv5_1,mv6_1,mv7_1,mv8_1,mv9_1,
     &mv10_1,rms_fac1
       real*8 mv1_2,mv2_2,mv3_2,mv4_2,mv5_2,mv6_2,mv7_2,mv8_2,mv9_2,
     &mv10_2,rms_fac2
       real*8 mv1_3,mv2_3,mv3_3,mv4_3,mv5_3,mv6_3,mv7_3,mv8_3,mv9_3,
     &mv10_3,rms_fac3
       integer A,GMAX_SW,G_SW,PRE_BH
       common nproc,myid
       common/seeding/ num
       common/same1/ samebc,MAXGEN
       common/mpigen2/par_num
       common/engtype/ type
       common/alloy1/ A
       common/min1/ min
       common/move1/ mv1_1,mv2_1,mv3_1,mv4_1,mv5_1,mv6_1,mv7_1,mv8_1,
     &mv9_1,mv10_1,rms_fac1
       common/move2/ mv1_2,mv2_2,mv3_2,mv4_2,mv5_2,mv6_2,mv7_2,mv8_2,
     &mv9_2,mv10_2,rms_fac2
       common/move3/ mv1_3,mv2_3,mv3_3,mv4_3,mv5_3,mv6_3,mv7_3,mv8_3,
     &mv9_3,mv10_3,rms_fac3
       common/move4/ mv1,mv2,mv3,mv4,mv5,mv6,mv7,mv8,mv9,mv10
C===============================
       integer J0,J1,J2,J3,J4,J5,J6,J7,mpierror
       integer np(ens_num),gen
       integer recover_rec,naccept,rec_sw,recover_sw,pdb_rec,
     &DG_RECOVER
       common/pdbrec/ pdb_rec,dim3
       common/accept1/ naccept
C===============================
       integer ist1,iend1,ist7,iend7
       integer gdisp1(0:node),gdisp2(0:node),gdisp3(0:node),
     &gdisp4(0:node),gdisp5(0:node),gdisp6(0:node),gdisp7(0:node),
     &gcount1(0:node),gcount2(0:node),gcount3(0:node),gcount4(0:node),
     &gcount5(0:node),gcount6(0:node),gcount7(0:node),mykount1,mykount2,
     &mykount3,mykount4,mykount5,mykount6,mykount7
       integer p1,p2,p3,p4,r1,r2,r3,r4,r5
       common/seeding2/ p1,p2,p3,p4
C==============================
       integer ex(ens_num*2),min_s,min_t,BH_SW1,BH_SW2,BH_SW3,BH_SW
       integer ken(ibin),ERANGE_SW,temp_ex,sw_ensnum,sw_parnum,sw_same
       integer STEPMAX1,MULMAX1,STEPMAX2,MULMAX2,STEPMAX3,MULMAX3
       integer sw_ensnum1,sw_ensnum2,sw_ensnum3,sw_parnum1,sw_parnum2,
     &sw_parnum3,sw_same1,sw_same2,sw_same3,view_recursion
       real*8 DEMIN,DEMAX,TDE,
C     &tde_node(0:ens_num/node+2),
      
     &tde_node(0:100),tde_total(ens_num)
       real*8 ex_temp(ens_num*2),ex_node(0:ens_num*2),
     &sw_node(0:ens_num),sw_temp(ens_num)
       real*8 outx(6000),x_dih(6000),grad(ndim)
       real*8 xinib(ndim),
C     &x_node(0:ens_num*ndim/node+2),
     &x_node(0:3000),
     &x_total(ens_num*ndim),xini(ens_num,ndim),xw(ens_num*ndim),
     &x_temp(ens_num*ndim),xinib_temp(ndim),
     &xinib_temp2(ndim),xinib_temp3(ndim)
       real*8 upot,pot,
C     &pot_node(0:ens_num/node+2)
     &pot_node(0:3000)
     &,pot_total(ens_num)
     &,pot_temp(ens_num),lmin,esys2
       real*8 s(ibin),
C     &s_node(0:ens_num*ibin/node+2),
     &s_node(0:5000),
     &s_total(ens_num*ibin),s_temp(ens_num*ibin)
       real*8 beta(ibin),
C     &beta0_node(0:ens_num/node+2),
     &beta0_node(0:5000),
     &beta0_total(ens_num),beta0_temp(ens_num)
       real*8 g_all(ibin),g_all_out(ibin),
C     &g_all_node(0:ens_num*ibin/node+2),
     &g_all_node(0:5000),
     &g_all_total(ens_num*ibin),g_all_temp(ens_num*ibin)
       real*8 
C     &temp_node(0:ens_num/node+2),
     &temp_node(0:5000),temp_total(ens_num),
     &temp_temp(ens_num),temp2
       real*8 al(ibin)
       real*8 s_out(ibin),beta_out(ibin),al_out(ibin)
C==============================
       real*8 iini,fini,igen,fgen,bgen(MAXGEN),per,fit(ens_num)
       real*8 initime,endtime,rnd,esys
       real*8 accrat,EP,fitg(ens_num),RMS,rms_start,rms_end,rms_dx1,
     &rms_dx2,rms_rate,GMAX,GMAX_TEMP,EMC
       real*8 recover_ratio
       real*8 pot_rec,rec_ratio
       real*8 scalefac,temp,tempi,tempf,dtemp,tempfd
       character name*15,name2*4,name3*15,name6*6,pdb_name*6,s_name*8
       character name8*19,kinda*2,kindb*2,kindab*2,sed*3
       real*8 RATE,BETA0,EMAX,EMIN,DE,EMIN_total(ens_num),
C     &EMIN_node(0:ens_num/node+2)
     &EMIN_node(0:100),EMIN_temp(ens_num)
       real*8 ASTEP,STEP,DGUESS,GTOL
       logical cflag,recursion_sw,EVAP
       common/bmin1/ DGUESS,ASTEP,STEP,GTOL
       common/name1/ name6
       common/pdbname1/ pdb_name
C       common/accuracy1/ r
       common/kinda/ kinda
       common/kindb/ kindb
       common/kindab/ kindab
       common/accept2/ accrat
       common/TEMP1/ scalefac,temp
       common/TEMP2/ tempi,tempf
       common/TEMP3/ temp_ex
       common/anneal/ tempfd
       common/UPPER/BETA0,N0,STEPMAX,MULMAX
       common/recover1/ recover_sw,recover_rec
       common/rec_sw1/ rec_sw
       common/rec_sw2/ rec_ratio
       common/recover2/ recover_ratio
       common/gabh1/ STEPMAX1,MULMAX1,STEPMAX2,MULMAX2,STEPMAX3,MULMAX3
       common/guide3/ DG_RECOVER
       common/RMS1/ RMS
       common/RMS2/ rms_start,rms_end,rms_dx1,rms_dx2,rms_rate
       common/RMS3/ rms_fac
       common/bmin3/ GMAX,GMAX_SW
       common/MULTIC1/ EMAX,EMIN
       common/MULTIC2/ ERANGE_SW,PRE_BH
       common/MULTIC3/ DE
       common/BH1/ BH_SW
       common/BH2/ BH_SW1,BH_SW2,BH_SW3
       common/BH3/ DEMIN,DEMAX
       common/BH4/ TDE
       common/SW_GA0/ sw_ensnum,sw_parnum,sw_same
       common/SW_GA1/ sw_ensnum1,sw_parnum1,sw_same1,sw_ensnum2,
     &sw_parnum2,sw_same2,sw_ensnum3,sw_parnum3,sw_same3
       COMMON/VIEW1/ view_recursion
       COMMON/WEIGHT1/G_SW
130    format(I2,A1,I3,A1,I4,F11.7,A1,F11.7,A1,
     &F8.4,A1,F8.6,A1,F8.6,A1,F8.6)
131    format(A4,F15.10,F15.10,F15.10)
132    format(A12,F25.13,A5,F14.8,A9)
133    format(I5,A13,F25.13,A6,I2,A5)
134    format(A4,F15.10)
135    format(I5,A4,F14.8,A5,F25.13,A1,I2,A5)
136    format(I4,F25.13,I5)
137    format(F25.13,F15.10)
138    format(I4,F9.4,1X, F9.4, 1X, F9.4,2X,F9.4,1X,I6)
139    format(I6,F15.10,F15.10,F15.10)
140    format(A6,I3,F14.8,A5,F25.13,A1,I2,A5) 
141    format(F25.13,F25.13)
142    format(a4,I7, 2X,a2,2X,  a3, I6,4X,3(F8.3))
143    format(A10,F25.13,I6) 
144    format(A3)  
145    format(F25.13)
       if(ens_num.eq.1.or.MAXGEN.le.1)then
         view_recursion=1
       else
         view_recursion=0
       endif
       EVAP=.false.      
       initime=MPI_WTIME()*10.D0
       dtemp=(tempf-tempi)/dble(ens_num)
       upot=9999.D0
C==================Parallel=============================       
       do J0=0,node-1
         call startend(J0,node,1,ens_num,ist1,iend1)
         gdisp1(J0)=ist1-1
         gdisp2(J0)=(ist1-1)*ndim
         gdisp3(J0)=(ist1-1)*ibin
         gcount1(J0)=iend1-ist1+1
         gcount2(J0)=(iend1-ist1+1)*ndim
         gcount3(J0)=(iend1-ist1+1)*ibin
       enddo
       mykount1=gcount1(myid)
       mykount2=gcount2(myid)
       mykount3=gcount3(myid)
       pot_rec=0.D0
       ndim2=ndim
****NDIM2 for Carticen coords output,ndim for caculation recording*******
       if(myid.eq.0) then
         iini=MPI_WTIME()*10.D0
       endif
       if (num.ne.0) then
         sed="sed"
         p1=num/1000
         p2=num/100
         p2=mod(p2,10)
         p3=num/10
         p3=mod(p3,10)
         p4=mod(num,10)
       else
         sed="nos"
       endif
       call mpi_barrier(mpi_comm_world,mpierror)
       if (recover_sw.eq.1)then
         if (BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
           if (ERANGE_SW.eq.0)then
             TDE=EMAX-EMIN
             DE=TDE/dble(ibin)
           else if(ens_num.eq.1)then
             TDE=DEMIN
           else if(ERANGE_SW.eq.1)then
             TDE=DEMIN+(dabs(DEMIN-DEMAX)/dble(ens_num-1))*
     &dble(myid*gcount1(myid)+(J4-1))
           endif
         endif
         if (GMAX_SW.eq.3)then
           GMAX=(rms_end+rms_dx1*dble(myid*gcount1(myid)+(J4-1)))
     &*rms_rate
         endif
         temp=(tempi+dtemp*dble(myid*gcount1(myid)+(J4-1)))
         if(BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
           BETA0=1.D0/(temp-tempfd)
         write(*,*) "temp=",temp,",tempfd=",tempfd,",beta0=",beta0
           if(BETA0.le.0.D0)then
             BETA0=0.00001D0
           endif
         endif
         do J5=1,ibin
           ken(J5)=0.D0
           beta(J5)=BETA0
           al(J5)=0.D0
           s(J5)=0.D0
           g_all(J5)=0.D0
         enddo
         if (myid.eq.0)then
           write(*,*) "@ Recover mode switch on!!"
           write(*,*) "@ Reading recover file ",name6//"_recover.rec"
           open(49,file=name6//"_recover.rec",status="old")
           read (49,*) J7
           open(50,file=name6//"_recover.bak",status="replace")
           write(*,*) "@ Start from ",J7,"' generation!"
           write(50,*) J7
           do J4=1,ens_num
             do J5=1,ndim
               read(49,*) pot_total(J4),xini(J4,J5)
               write(50,137) pot_total(J4),xini(J4,J5)
             enddo
             if(ERANGE_SW.eq.1)then
               EMIN_total(J4)=pot_total(J4)-(TDE/2.D0)
             endif            
           enddo
           do J4=1,ens_num
             read(49,*) temp_total(J4)
             write(50,145) temp_total(J4)
           enddo
           if(BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
             do J4=1,ens_num
               do J5=1,ibin
                 J2=J2+1
                 read(49,*) beta0_total(J4),s_total((J4-1)*ibin+J5)
                 write(50,141) beta0_total(J4),s_total((J4-1)*ibin+J5)
               enddo
             enddo
             if(ERANGE_SW.ne.0)then
               do J4=1,ens_num
                 read(49,*) EMIN_total(J4),tde_total(J4),EMAX,EMIN
                 write(50,141) EMIN_total(J4),tde_total(J4),EMAX,EMIN
               enddo
             endif
             write(*,*) "@ The MBH parameters reading complete!!"
           endif
           close(49)
           close(50)
           write(*,*) "@ Reading All data complete!!"
         endif
       else
         J7=0
       endif
       call mpi_barrier(mpi_comm_world,mpierror)
       J1=0
       J2=0
       J3=0
       if (recover_sw.ne.1)then
         if (myid.eq.0)then
           write(*,*)
           write(*,*) "Start Caculating Ensemble"
           write(*,*)
         endif
         GMAX_TEMP=GMAX
         do J4=1,gcount1(myid)
***********The 1st Canonical Basin-Hopping**********************
           call dim_gen(xinib,ndim,seed)
C           call POTENTIAL(ndim,xinib,grad,esys,RMS,.false.,EVAP)
           STEPMAX=STEPMAX1
           MULMAX=MULMAX1
           BH_SW=BH_SW1
           sw_ensnum=sw_ensnum1
           sw_parnum=sw_parnum1
           sw_same=sw_same1
           mv1=mv1_1
           mv2=mv2_1
           mv3=mv3_1
           mv4=mv4_1
           mv5=mv5_1
           mv6=mv6_1
           mv7=mv7_1
           mv8=mv8_1
           mv9=mv9_1
           mv10=mv10_1
           rms_fac=rms_fac1
           if (BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
             if (ERANGE_SW.eq.0)then
               TDE=EMAX-EMIN
               DE=TDE/dble(ibin)
             else if(ens_num.eq.1)then
               TDE=DEMIN
             else if(ERANGE_SW.eq.1)then
               TDE=DEMIN+(dabs(DEMIN-DEMAX)/dble(ens_num-1))*
     &dble(myid*gcount1(myid)+(J4-1))
             endif
           endif
           if (GMAX_SW.eq.3)then
             GMAX=(rms_end+rms_dx1*dble(myid*gcount1(myid)+(J4-1)))
     &*rms_rate
           endif
           temp=(tempi+dtemp*dble(myid*gcount1(myid)+(J4-1)))
           if(BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
             BETA0=1.D0/(temp-tempfd)
C           write(*,*) "temp=",temp,",tempfd=",tempfd,",beta0=",beta0
             if(BETA0.le.0.D0)then
               BETA0=0.00001D0
             endif
           endif
           do J5=1,ibin
             ken(J5)=0.D0
             beta(J5)=BETA0
             al(J5)=0.D0
             s(J5)=0.D0
             g_all(J5)=0.D0
           enddo  
C           CALL step_move(xinib,xinib,NDIM,SEED,GBIN,
C     &GUIDE,GBIN_TOTAL)    
C          call POTENTIAL(ndim,xinib,grad,esys2,RMS,.false.,EVAP)   
c           write(*,*) "esys2=",esys2
           call LOCAL_MIN(xinib,pot,ndim,seed,cflag,min_s,min_t,
     &EVAP)
           do J5=1,ndim
             xinib_temp3(J5)=xinib(J5)
             xinib_temp2(J5)=xinib(J5)
             xinib_temp(J5)=xinib(J5)
           enddo
           esys=pot
           lmin=pot
           if (BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
             if(ERANGE_SW.eq.1)then
               EMIN=pot-(TDE/2.D0)
               EMAX=EMIN+TDE
             endif  
             if(ERANGE_SW.eq.3)then
               STEPMAX=PRE_BH
               MULMAX=1
               BH_SW=2
           write(*,*) "START PRE-BH PROCEDURE,NUMBER OF STEP=",PRE_BH
            call RECURSION(xinib_temp3,xinib_temp2,ndim,esys,lmin,
     &seed,RATE,ibin,beta,al,ken)
               BH_SW=BH_SW1
               STEPMAX=STEPMAX1
               MULMAX=MULMAX1
               
             endif
             write(*,*)"==> PC=",myid,",MAXSTEP=",STEPMAX,",TEMP=",TEMP,
     &",BETA0=",BETA0
             write(*,*) "EMIN=",EMIN,",EMAX=",EMAX,",TDE",TDE
           else
             write(*,*)"==> PC=",myid,",MAXSTEP=",STEPMAX,",TEMP=",TEMP
           endif
           RATE=0.D0
           do J5=1,MULMAX
             call RECURSION(xinib_temp,xinib_temp2,ndim,esys,lmin,
     &seed,RATE,ibin,beta,al,ken)
             if (lmin.lt.pot) then
               pot=lmin
               do J6=1,ndim
                 xinib(J6)=xinib_temp2(J6)
               enddo
             endif
             if (GMAX_SW.eq.1.and.RMS.le.GMAX) then
       write(*,*)"ENS=",J4,",RUN=",J5,",lmin=",esys,",gmin=",pot,
     *",RMS=",RMS
             else
       write(*,*)"ENS=",J4,",RUN=",J5,",lmin=",esys,",gmin=",pot,
     &",GMAX=",GMAX
             endif
             if (BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
               call NEWAB(ibin,pot,ken,beta,al,s,g_all)
             endif
           enddo
           call move(RATE)
33         do J5=1,ndim
             J1=J1+1
             x_node(J1)=xinib(J5)
           enddo
           J3=J3+1
           temp_node(J3)=temp
           if(BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
C             do J6=1,ibin
C               ken_test(J6)=ken(J6)
C             enddo
             do J5=1,ibin
               J2=J2+1
C               beta_node(J2)=beta(J5)
C               al_node(J2)=al(J5)
               s_node(J2)=s(J5)
               g_all_node(J2)=g_all(J5)
C               write(*,*) s(J5),al(J5),beta(J5)
             enddo
             beta0_node(J3)=BETA0
             if(ERANGE_SW.ne.0)then
               tde_node(J3)=TDE
               EMIN_node(J3)=EMIN
             endif
           endif
           pot_node(J3)=pot
           if(RATE.gt.accrat)then
      write(*,*) "ENS=",J4,",",RATE,"(rate) >",accrat,"(acc),STEP="
     &,STEP,",ASTEP=",ASTEP,",IN"
           else if(RATE.eq.accrat)then
      write(*,*) "ENS=",J4,",",RATE,"(rate) =",accrat,"(acc),STEP="
     &,STEP,",ASTEP=",ASTEP,",*1"
           else           
      write(*,*) "ENS=",J4,",",RATE,"(rate) <",accrat,"(acc),STEP="
     &,STEP,",ASTEP=",ASTEP,",DE"
           endif
         enddo
************Distribute all information for all nodes**************
         call mpi_barrier(mpi_comm_world,mpierror)
         call mpi_gatherv(pot_node(1),mykount1,mpi_real8,pot_total,
     &gcount1,gdisp1,mpi_real8,0,mpi_comm_world,mpierror)
         if (BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
           if(ERANGE_SW.ne.0)then
             call mpi_gatherv(tde_node(1),mykount1,mpi_real8,tde_total,
     &gcount1,gdisp1,mpi_real8,0,mpi_comm_world,mpierror)
            call mpi_gatherv(EMIN_node(1),mykount1,mpi_real8,EMIN_total,
     &gcount1,gdisp1,mpi_real8,0,mpi_comm_world,mpierror)
           endif
          call mpi_gatherv(beta0_node(1),mykount1,mpi_real8,beta0_total,
     &gcount1,gdisp1,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_gatherv(s_node(1),mykount3,mpi_real8,s_total,gcount3,
     &gdisp3,mpi_real8,0,mpi_comm_world,mpierror)
          call mpi_gatherv(g_all_node(1),mykount3,mpi_real8,
     &g_all_total,gcount3,gdisp3,mpi_real8,0,mpi_comm_world,mpierror)
         endif
         call mpi_gatherv(temp_node(1),mykount1,mpi_real8,temp_total,
     &gcount1,gdisp1,mpi_real8,0,mpi_comm_world,mpierror)
         call mpi_gatherv(x_node(1),mykount2,mpi_real8,x_total,gcount2,
     &gdisp2,mpi_real8,0,mpi_comm_world,mpierror)
         call mpi_barrier(mpi_comm_world,mpierror)
         if(myid.eq.0)then
           J1=0
           do J4=1,ens_num
             do J5=1,ndim
               J1=J1+1
               xini(J4,J5)=x_total(J1)
C              xinib(J5)=xini(J4,J5)
             enddo
           enddo
         endif
       endif
*************End of the 1st  Basin-Hopping*********************
       if(myid.eq.0)then
         call ene_sort(pot_total,same,ens_num,np)
         do J1=1,ndim
           outx(J1)=xini(np(1),J1)
         enddo
         do J1=1,ibin
           s_out(J1)=s_total((np(1)-1)*ibin+J1)
           g_all_out(J1)=g_all_total((np(1)-1)*ibin+J1)
         enddo
         upot=pot_total(np(1))
         call NEWAB(ibin,upot,ken,beta_out,al_out,s_out,g_all_out)
         fini=MPI_WTIME()*10.D0
         if (min.eq.1)then
           name=name6//"_simp."//sed
         else if(min.eq.2)then
           name=name6//"_conj."//sed
         else if(min.eq.3)then
           name=name6//"_cosp."//sed
         else if(min.eq.4)then
           name=name6//"_spco."//sed
         else
           name=name6//"_lbfg."//sed
         endif
         open(51,file=name,status="replace")
         open(52,file=name//".ini",status="replace")
         open(53,file=name//".genper",status="replace")
*****************<PDB ANIMATION RECORDING>*********************
         if (pdb_rec.eq.1) then
           open(61,file=name6//"_rec.pdb",status="replace")
           write(61,143) "HEADER    ",pot_total(np(1)),0
           j1=0
**************for dih angle,call dih_gen******************
*****dim3=1:call dih/dim3=0:origin-->>N**************
*****dim3=0 or 2 -->>3N***************************
           if (dim3.eq.1) then
             do J5=1,ndim
               xinib(J5)=xini(np(1),J5)
             enddo
             call dih_gen(ndim,xinib,x_dih)
             ndim2=ndim*3+9
           endif
           do J5=1,ndim2,3
             j1=j1+1
             if(dim3.eq.0)then
               x_dih(J5)=xini(np(1),J5)
               x_dih(J5+1)=xini(np(1),J5+1)
               x_dih(J5+2)=xini(np(1),J5+2)
             endif
             if((type.eq.2.or.type.eq.9).and.(J5.gt.(A*3)))then     
        write(61,142)'ATOM', J1,kindb, 'UNK',J1,(x_dih(J5+J6),J6=0,2)
             else
        write(61,142)'ATOM', J1,kinda, 'UNK',J1,(x_dih(J5+J6),J6=0,2)
              endif
           enddo
           write(61,144) "END"
         endif
**************END OF PDB ANIMATION REC********************
         if (recover_sw.eq.1)then
           write(51,*) "Recovered job,start from ",J7,"'gen"
         endif
         write(51,*) ndim2/3," dimension,use seeding,min method=",min
         write(51,*) "Initial time:",dabs(fini-iini)/10.D0," sec"  
         write(52,*) ndim2/3,dabs(fini-iini)/10.D0
         call flush(6)
         close(52)
**************************************************<<<<<<GA>>>>>>*************************************************
         if (MAXGEN.le.1)then
           write(*,*)
           write(*,*) "USE ONLY 1st RUN without GA(2nd RUN)"
         else
           write(*,*) 
           write(*,*) "Start GA!!"
         endif
       endif
       same=0
       J0=0
       STEPMAX=STEPMAX2
       MULMAX=MULMAX2
       sw_ensnum=sw_ensnum2
       sw_parnum=sw_parnum2
       sw_same=sw_same2
       mv1=mv1_2
       mv2=mv2_2
       mv3=mv3_2
       mv4=mv4_2
       mv5=mv5_2 
       mv6=mv6_2
       mv7=mv7_2
       mv8=mv8_2
       mv9=mv9_2
       mv10=mv10_2
       rms_fac=rms_fac2
       if (BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
         if(ERANGE_SW.ne.0)then
           call mpi_bcast(tde_total,ens_num,mpi_real8,0,mpi_comm_world,
     &mpierror)
C         if(myid.eq.0)write(*,*) "TEST1"
         endif
       endif
       do while(J0.lt.MAXGEN.and.same.le.samebc.and.MAXGEN.ne.1)
         J0=J0+1
         if(myid.eq.0)then
           igen=MPI_WTIME()*10.D0
*****************MPIGA START!!***************************
           if (J0.eq.1.or.mod(J0,naccept).eq.0)then
             call mpi_gen(xini,pot_total,ens_num,ndim,np,seed,fit,
     &xw,ex,almin)
             almin_temp(1)=dble(almin)
           endif
         endif
         call mpi_barrier(mpi_comm_world,mpierror)
         call mpi_bcast(pot_total,ens_num,mpi_real8,0,mpi_comm_world,
     &mpierror)
         call mpi_bcast(almin_temp,1,mpi_real8,0,mpi_comm_world,
     &mpierror)
         almin=int(almin_temp(1))
         do J4=1,ens_num
           np_temp(J4)=dble(np(J4))
         enddo
         call mpi_bcast(np_temp,ens_num,mpi_real8,0,mpi_comm_world,
     &mpierror)
         do J4=1,ens_num
           np(J4)=int(np_temp(J4))
         enddo
         call mpi_bcast(fit,ens_num,mpi_real8,0,mpi_comm_world,
     &mpierror)
         if(BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
         call mpi_bcast(beta0_total,ens_num,mpi_real8,0,mpi_comm_world,
     &mpierror)
         call mpi_bcast(s_total,ens_num*ibin,mpi_real8,0,mpi_comm_world,
     &mpierror)
C         call mpi_bcast(g_all_total,ens_num*ibin,mpi_real8,0,
C     &mpi_comm_world,mpierror)
           if(ERANGE_SW.ne.0)then
             call mpi_bcast(EMIN_total,ens_num,mpi_real8,0,
     &mpi_comm_world,mpierror)
           endif
         endif
         call mpi_bcast(temp_total,ens_num,mpi_real8,0,mpi_comm_world,
     &mpierror)
         dtemp=(tempf-tempi)/dble(ens_num-almin)
         do J4=0,node-1
           call startend(J4,node,1,(ens_num-almin),ist7,iend7)
           gdisp4(J4)=(ist7-1)*ndim
           gcount4(J4)=(iend7-ist7+1)*ndim
           gdisp5(J4)=(ist7-1)*2
           gcount5(J4)=(iend7-ist7+1)*2
           gdisp6(J4)=(ist7-1)*ibin
           gcount6(J4)=(iend7-ist7+1)*ibin
           gdisp7(J4)=ist7-1
           gcount7(J4)=iend7-ist7+1
         enddo
         mykount4=gcount4(myid)
         mykount5=gcount5(myid)
         mykount6=gcount6(myid)
         mykount7=gcount7(myid)
         call mpi_barrier(mpi_comm_world,mpierror)
********************************************************
         call mpi_scatterv(xw,gcount4,gdisp4,mpi_real8,x_node(1),
     &mykount4,mpi_real8,0,mpi_comm_world,mpierror)     
         do J4=1,ens_num*2
           ex_temp(J4)=dble(ex(J4))
         enddo
         call mpi_scatterv(ex_temp,gcount5,gdisp5,mpi_real8,
     &ex_node(1),mykount5,mpi_real8,0,mpi_comm_world,mpierror)
         call mpi_barrier(mpi_comm_world,mpierror)
**********minized the potential in mpi****************** 
         J1=0
         J2=0 
         J3=0
         BH_SW=BH_SW2
         do J4=1,gcount7(myid)
           sw_node(J4)=0.D0
           do J5=1,ndim
             xinib(J5)=x_node((J4-1)*ndim+J5)
           enddo
           if (GMAX_SW.eq.3)then
             GMAX=(rms_end+rms_dx2*dble(myid*gcount7(myid)+(J4-1)))*
     &rms_rate
           endif
           call LOCAL_MIN(xinib,pot,ndim,seed,cflag,min_s,min_t,EVAP
     &)
           recursion_sw=.false.
           if (STEPMAX2.eq.1.and.MULMAX2.eq.1)then
             recursion_sw=.true.
             sw_node(J4)=1.D0
           else
             if (temp_ex.eq.2)then
               temp=(tempi+dtemp*dble(myid*gcount7(myid)+(J4-1)))
             else
               temp=temp_total(int(ex_node((J4)*2-1)))
             endif
             if (dexp((pot_total(int(ex_node((J4)*2-1)))-pot)/temp).gt.
     &rnd(seed))then
               if (BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
                 if (ERANGE_SW.ne.0)then
                   TDE=tde_total(int(ex_node(J4*2-1)))
                   EMIN=EMIN_total(int(ex_node(J4*2-1)))
                   EMAX=EMIN+TDE
                 endif
                 BETA0=beta0_total(int(ex_node(J4*2-1)))
                 do J5=1,ibin
C                 beta(J5)=beta_total(((int(ex_node(J4))*2)-2)*ibin+J5)
C                 al(J5)=al_total(((int(ex_node(J4))*2)-2)*ibin+J5)
                   s(J5)=s_total(((int(ex_node(J4))*2)-2)*ibin+J5)
              g_all(J5)=g_all_total(((int(ex_node(J4))*2)-2)*ibin+J5)
C               ken(J5)=int(ken_total(((int(ex_node(J4))*2)-2)*ibin+J5))
                 enddo
               endif
               recursion_sw=.true.
               sw_node(J4)=1.D0
             endif
             if (ex_node(J4*2).ne.0.D0)then
               if (temp_ex.eq.2)then
                 temp2=(tempi+dtemp*dble(myid*gcount7(myid)+(J4-1)))
               else
                 temp2=temp_total(int(ex_node((J4)*2)))
               endif
               if (dexp((pot_total(int(ex_node(J4*2)))-pot)/temp2).gt.
     &rnd(seed))then
                 if (fit((int(ex_node(J4))*2)).gt.
     &fit((int(ex_node(J4))*2-1)))then
                   temp=temp2
                   if (BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
                     if (ERANGE_SW.ne.0)then
                       TDE=tde_total(int(ex_node(J4*2)))
                       EMIN=EMIN_total(int(ex_node(J4*2)))
                       EMAX=EMIN+TDE
                     endif
                     BETA0=beta0_total(int(ex_node(J4*2)))
                     do J5=1,ibin
C                   beta(J5)=beta_total(((int(ex_node(J4))*2)-1)*ibin+J5)
C                   al(J5)=al_total(((int(ex_node(J4))*2)-1)*ibin+J5)
                       s(J5)=s_total(((int(ex_node(J4))*2)-1)*ibin+J5)
               g_all(J5)=g_all_total(((int(ex_node(J4))*2)-1)*ibin+J5)
C                  ken(J5)=int(ken_total(((int(ex_node(J4))*2)-1)*ibin+J5))
                     enddo
                   endif
                   recursion_sw=.true.
                   sw_node(J4)=1.D0
                 endif
               endif
             endif
           endif
           if (recursion_sw)then
             do J5=1,ndim
               xinib_temp(J5)=xinib(J5)
               xinib_temp2(J5)=xinib(J5)
             enddo
             esys=pot
             lmin=pot
             if (BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
               call NEWAB(ibin,pot,ken,beta,al,s,g_all)
               write(*,*)"==> PC=",myid,",BETA0=",BETA0,",TEMP=",
     &TEMP," X ",scalefac
               write(*,*) "EMIN=",EMIN,",EMAX=",EMAX,",TDE=",TDE
             else
               write(*,*)"==> PC=",myid,",MAXSTEP=",STEPMAX,",TEMP=",
     &TEMP," X ",scalefac
             endif
             RATE=0.D0
             do J5=1,MULMAX
               call RECURSION(xinib_temp,xinib_temp2,ndim,esys,lmin,
     &seed,RATE,ibin,beta,al,ken)
               if (lmin.le.pot)then
                 do J6=1,ndim
                   xinib(J6)=xinib_temp2(J6)
                 enddo
                 pot=lmin
               endif
               if (GMAX_SW.eq.1.and.RMS.le.GMAX) then
                 write(*,*)"ENS=",J4,",Run=",J5,",lmin=",esys,
     &",gmin=",pot,",RMS=",RMS
               else
                 write(*,*)"ENS=",J4,",Run=",J5,",lmin=",esys,
     &",gmin=",pot,",GMAX=",GMAX
               endif
               if(BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.
     &or.BH_SW1.eq.4)then
                 call NEWAB(ibin,pot,ken,beta,al,s,g_all)
               endif
             enddo
             call move(RATE)
           endif
           do J5=1,ndim
             J1=J1+1
             x_node(J1)=xinib(J5)
           enddo
           J3=J3+1
           pot_node(J3)=pot
           temp_node(J3)=temp
           if (BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
             do J5=1,ibin
               J2=J2+1
               s_node(J2)=s(J5)
               g_all_node(J2)=g_all(J5)
             enddo
             if (ERANGE_SW.ne.0)then
               EMIN_node(J3)=EMIN
C               tde_node(J3)=TDE
             endif
           endif
           if(RATE.gt.accrat)then
      write(*,*) "ENS=",J4,",",RATE,"(rate) >",accrat,"(acc),STEP="
     &,STEP,",ASTEP=",ASTEP,",IN"
           else if(RATE.eq.accrat)then
      write(*,*) "ENS=",J4,",",RATE,"(rate) =",accrat,"(acc),STEP="
     &,STEP,",ASTEP=",ASTEP,",*1"
           else
      write(*,*) "ENS=",J4,",",RATE,"(rate) <",accrat,"(acc),STEP="
     &,STEP,",ASTEP=",ASTEP,",DE"
           endif
         enddo  
C         write(*,*) "before gatherv",myid,mykount4,mykount7
C         do J3=1,mykount4
C           write(*,*) J3,"xinib_node=",x_node(J3),myid
C         enddo
         call mpi_barrier(mpi_comm_world,mpierror)
         call mpi_gatherv(x_node(1),mykount4,mpi_real8,x_temp,
     &gcount4,gdisp4,mpi_real8,0,mpi_comm_world,mpierror)
         call mpi_gatherv(pot_node(1),mykount7,mpi_real8,pot_temp,
     &gcount7,gdisp7,mpi_real8,0,mpi_comm_world,mpierror)
         call mpi_gatherv(sw_node(1),mykount7,mpi_real8,sw_temp,
     &gcount7,gdisp7,mpi_real8,0,mpi_comm_world,mpierror)
C         call mpi_gatherv(ken_node(1),mykount6,mpi_real8,ken_temp,
C     &gcount6,gdisp6,mpi_real8,0,mpi_comm_world,mpierror)
         if(BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
           call mpi_gatherv(s_node(1),mykount6,mpi_real8,s_temp,
     &gcount6,gdisp6,mpi_real8,0,mpi_comm_world,mpierror)
           call mpi_gatherv(g_all_node(1),mykount6,mpi_real8,
     &g_all_temp,gcount6,gdisp6,mpi_real8,0,mpi_comm_world,mpierror)
C         call mpi_gatherv(beta_node(1),mykount6,mpi_real8,beta_temp,
C     &gcount6,gdisp6,mpi_real8,0,mpi_comm_world,mpierror)
C         call mpi_gatherv(al_node(1),mykount6,mpi_real8,al_temp,
C     &gcount6,gdisp6,mpi_real8,0,mpi_comm_world,mpierror)
           if (ERANGE_SW.ne.0)then
             call mpi_gatherv(EMIN_node(1),mykount7,mpi_real8,
     &EMIN_temp,gcount7,gdisp7,mpi_real8,0,mpi_comm_world,mpierror)
           endif
         endif
         call mpi_gatherv(temp_node(1),mykount7,mpi_real8,temp_temp,
     &gcount7,gdisp7,mpi_real8,0,mpi_comm_world,mpierror)
C         call mpi_allreduce(guide,guide_temp,gbin*gbin*gbin,mpi_integer,
C     &mpi_sum,mpi_comm_world,mpierror)
C         call mpi_allreduce(gbin_total,gbin_all,1,mpi_integer,
C     &mpi_sum,mpi_comm_world,mpierror)
         call mpi_barrier(mpi_comm_world,mpierror)
C         do J5=1,gbin*gbin*gbin
C           guide(J5)=0
C           guide_all(J5)=guide_all(J5)+guide_temp(J5)
C         enddo
C         gbin_total=0
C         write(*,*) "gbin_total=",gbin_all
         if(myid.eq.0)then
           J1=0
           J2=0
           J3=0
           do J6=almin+1,ens_num
             J3=J3+1
             if(sw_temp(J3).eq.1.D0)then
               do J5=1,ndim
                 J1=J1+1
                 xini(np(J6),J5)=x_temp(J1)
                 xinib(J5)=x_temp(J1)
               enddo
               if (BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
                 do J5=1,ibin
                   J2=J2+1
C                 ken_total((np(J6)-1)*ibin+J5)=ken_temp(J2)
                   s_total((np(J6)-1)*ibin+J5)=s_temp(J2)
                   g_all_total((np(J6)-1)*ibin+J5)=g_all_temp(J2)
C                 beta_total((np(J6)-1)*ibin+J5)=beta_temp(J2)
C                 al_total((np(J6)-1)*ibin+J5)=al_temp(J2)
C                 write(*,*) J2,"s_temp=",s_temp(J2)                 
                 enddo
                 if (ERANGE_SW.ne.0)then
                   EMIN_total(np(J6))=EMIN_temp(J3)
                 endif
               endif
               pot_total(np(J6))=pot_temp(J3)
               temp_total(np(J6))=temp_temp(J3)
C              write(*,*) J6,",pot=",pot_total(np(J6)),"=?=",esys2,np(J6)
             else
               J1=J1+ndim
               J2=J2+ibin
             endif
           enddo
           do J6=1,ens_num
             do J5=1,ndim
               xinib(J5)=xini(np(J6),J5)
             enddo
           enddo
******************MPIGA FINISH!!*************************
           fgen=MPI_WTIME()*10.D0
           bgen(J0)=dabs(fgen-igen)/10.D0         
           call ene_sort(pot_total,same,ens_num,np)
           write(*,*)
           if(recover_sw.eq.1)then
             write(*,*) "Generation=",J7,"+",J0,"=",J7+J0," ,same=",same
     &," ,Current temperature=",temp
           else
             write(*,*) "Generation=",J0," ,same=",same
     &," ,Current temperature=",temp
           endif
**************START RECORDING/RECOVERING*************************
           if(recover_rec.eq.1.and.((pot_rec-pot_total(np(1))).ge.
     &recover_ratio).or.J0.ge.MAXGEN)then
             open(54,file=name6//"_recover.rec",status="replace")
             write(54,*) J0+J7
             do J6=1,ens_num
               do J5=1,ndim
                 write(54,137) pot_total(np(J6)),(xini(np(J6),J5))
               enddo
             enddo
             do J6=1,ens_num
               write(54,145) temp_total(np(J6))
             enddo
             if(BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
               do J6=1,ens_num
                 do J5=1,ibin
           write(54,141) beta0_total(np(J6)),s_total((np(J6)-1)*ibin+J5)
                 enddo
               enddo
               if(ERANGE_SW.ne.0)then
                 do J6=1,ens_num
                   write(54,141) EMIN_total(np(J6)),tde_total(np(J6))
                 enddo
               endif
             endif
             close(54)
           endif
           if(rec_sw.eq.0)then
             open(55,file=name6//".rec",status="replace")
      write(*,132) "The min_pot=",pot_total(np(1)),"(ev),",
     &bgen(J0)," sec NREC"
      write(51,135) J0,"'gen",bgen(J0),"sec,",pot_total(np(1))," ",
     &same," NREC"
             write(55,136) ndim2/3,pot_total(np(1)),J0
             write(55,*)
             if (dim3.eq.0) then
               do J5=1,ndim,3
                 if((type.eq.2.or.type.eq.9).and.(J5.gt.(A*3)))then
                   write(55,131) kindb,xini(np(1),J5),xini(np(1),J5+1),
     &xini(np(1),J5+2)
                 else
                   write(55,131) kinda,xini(np(1),J5),xini(np(1),J5+1),
     &xini(np(1),J5+2)
                 endif
               enddo
             else
               do J5=1,ndim
                 write(55,134) kinda,xini(np(1),J5)
               enddo
             endif
             close(55)
           endif
           do J5=1,ens_num
             do J6=1,ndim
               xinib(J6)=xini(np(J5),J6)
             enddo
           enddo
           if(rec_sw.eq.1)then
             if((pot_rec-pot_total(np(1))).ge.rec_ratio.or.J0.
     &ge.MAXGEN)then
      write(*,132) "The min_pot=",pot_total(np(1)),"(ev),",
     &bgen(J0)," sec SREC"
       write(51,135) J0,"'gen",bgen(J0),"sec,",pot_total(np(1))," ",
     &same," SREC"
C   PDB ANIMATION RECORDING!!
               if (pdb_rec.eq.1)then
                 j1=0
                 write(61,143) "HEADER    ",pot_total(np(1)),J0
                 if (dim3.eq.1) then
                   do J5=1,ndim
                     xinib(J5)=xini(np(1),J5)
                   enddo
                   call dih_gen(ndim,xinib,x_dih)
                 endif
                 do J5=1,ndim2,3
                   j1=j1+1
                   if(dim3.eq.0)then
                     x_dih(J5)=xini(np(1),J5)
                     x_dih(J5+1)=xini(np(1),J5+1)
                     x_dih(J5+2)=xini(np(1),J5+2)
                   endif
                   if((type.eq.2.or.type.eq.9).and.(J5.gt.(A*3)))then
         write(61,142)'ATOM', J1,kindb, 'UNK',J1,(x_dih(J5+J6),J6=0,2)
                   else
         write(61,142)'ATOM', J1,kinda, 'UNK',J1,(x_dih(J5+J6),J6=0,2)
                   endif
                 enddo
                 write(61,144) "END"
                 open(62,file=name6//".rec",status="replace")
                 write(62,136) ndim2/3,pot_total(np(1)),J0
                 write(62,*)
                 if(dim3.eq.1)then
                   do J5=1,ndim
                     write(62,134) kinda,xini(np(1),J5)
                   enddo
                 else
                   do J5=1,ndim,3
                     if((type.eq.2.or.type.eq.9).and.(J5.gt.(A*3)))then
                       write(62,131) kindb,xini(np(1),J5),
     &xini(np(1),J5+1),xini(np(1),J5+2)
                     else
                       write(62,131) kinda,xini(np(1),J5),
     &xini(np(1),J5+1),xini(np(1),J5+2)
                     endif
                   enddo
                 endif
                 close(62)
               else  
                 r1=J0/10000
                 r2=J0/1000
                 r2=mod(r2,10)
                 r3=J0/100
                 r3=mod(r3,10)
                 r4=J0/10
                 r4=mod(r4,10)
                 r5=mod(J0,10)
            name8=name6//"_gen"//char(r1+48)//char(r2+48)//char(r3+48)
     &//char(r4+48)//char(r5+48)//".xyz"
                 open(47,file=name8,status="replace")
                 write(47,*) ndim2/3,pot_total(np(1)),J0
                 write(47,*)
                 if (dim3.eq.1)then
                
                 else   
                   do J5=1,ndim,3
                     if((type.eq.2.or.type.eq.9).and.(J5.gt.(A*3)))then
                       write(47,131) kindb,xini(np(1),J5),
     &xini(np(1),J5+1),xini(np(1),J5+2)
                     else
                       write(47,131) kinda,xini(np(1),J5),
     &xini(np(1),J5+1),xini(np(1),J5+2)
                     endif
                   enddo
                 endif
                 close(47)
               endif
             else
      write(*,132) "The min_pot=",pot_total(np(1)),"(ev),",
     &bgen(J0)," sec     "
      write(51,135) J0,"'gen",bgen(J0),"sec,",
     &pot_total(np(1))," ",same,"     "
             endif
             pot_rec=pot_total(np(1))
           endif
           if (rec_sw.eq.2)then
      write(*,132) "The min_pot=",pot_total(np(1)),"(ev),",
     &bgen(J0)," sec     "
      write(51,135) J0,"'gen",bgen(J0),"sec,",
     &pot_total(np(1))," ",same,"     "
**************END OF RECORDING*************************           
           endif  
           gen=J0
           write(*,*)
           call flush(6)
         endif                                            
         call mpi_barrier(mpi_comm_world,mpierror) 
         same_temp(1)=dble(same)
         call mpi_bcast(same_temp,1,mpi_real8,0,mpi_comm_world,
     &mpierror)
         same=int(same_temp(1))
         if (myid.eq.0)then
           if(pot_total(np(1)).lt.upot)then
             upot=pot_total(np(1))
             if(BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
               if(ERANGE_SW.ne.0)then
                 TDE=tde_total(np(1))
                 EMIN=EMIN_total(np(1))
                 EMAX=EMIN+TDE
               endif
               do J5=1,ibin
                 s_out(J5)=s_total((np(1)-1)*ibin+J5)
                 g_all_out(J5)=g_all_total((np(1)-1)*ibin+J5)
               enddo
             endif
             do J5=1,ndim
               outx(J5)=xini(np(1),J5)
             enddo
           endif
         endif
       enddo
C       call mpi_barrier(mpi_comm_world,mpierror)
*************END OF MPI GA**********************
10     if(myid.eq.0)then
C         do J5=1,ndim,3
C           write(*,*) outx(J5),outx(J5+1),outx(J5+2)
C         enddo
C        call CENTRE(outx,ndim/3)
*****************************************
         endtime=MPI_WTIME()*10.D0
         per=dabs(dabs(endtime-initime)-dabs(fini-iini))/
     &(dble(gen)*(10.D0))
         write(53,*) ndim2/3,per
         close(53)
         igen=MPI_WTIME()*10.D0
************RESET RECOVER FOR SAFTY****************8
         STEPMAX=STEPMAX3
         MULMAX=MULMAX3
         BH_SW=BH_SW3
         sw_ensnum=sw_ensnum3
         sw_parnum=sw_parnum3
         sw_same=sw_same3
         mv1=mv1_3
         mv2=mv2_3
         mv3=mv3_3
         mv4=mv4_3
         mv5=mv5_3
         mv6=mv6_3
         mv7=mv7_3
         mv8=mv8_3
         mv9=mv9_3
         mv10=mv10_3
         rms_fac=rms_fac3
         do J5=1,ndim
           xinib_temp(J5)=outx(J5)
           xinib_temp2(J5)=outx(J5)
         enddo
         esys=upot
         lmin=upot
         if(temp_ex.eq.2)then
           temp=tempi
         else
           temp=temp_total(np(1))
         endif
         write(*,*) "------------------------------------------------"
         write(*,*) "Final Run,MULMAX=",MULMAX,",STEPMAX=",STEPMAX
C         call POTENTIAL(ndim,outx,grad,esys2,RMS,.false.,EVAP)
C         write(*,*) "Before final,upot=",esys2
         if(BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
           call NEWAB(ibin,upot,ken,beta_out,al_out,s_out,g_all_out)
           do J5=1,ibin
             ken(J5)=0
C             write(*,*) "BETA=",beta_out(J5),"Al=",al_out(J5)
C     &,"S=",s_out(J5)
           enddo
           write(*,*)"Final ENS=",np(1),",MAXSTEP=",STEPMAX,
     &",BETA0=",BETA0,",TEMP=",TEMP
         else
           write(*,*)"Final ENS=",np(1),",MAXSTEP=",STEPMAX,
     &",TEMP=",TEMP
         endif
         do J5=1,MULMAX
      call RECURSION(xinib_temp,xinib_temp2,ndim,esys,lmin,seed,RATE,
     &ibin,beta_out,al_out,ken)
           if (lmin.le.upot)then
             upot=lmin
             do J6=1,ndim
               outx(J6)=xinib_temp2(J6)
             enddo
           endif
           write(51,140) "FINAL ",J5,dabs((fgen-igen)/10.D0),"sec,",
     &upot," ",same," RUNS"
           call flush(6)
           if(GMAX_SW.eq.1.and.RMS.le.GMAX) then
        write(*,*) "Run=",J5,",lmin=",lmin,",gmin=",upot,",RMS=",
     &RMS
           else
        write(*,*) "Run=",J5,",lmin=",lmin,",gmin=",upot,",GMAX=",
     &GMAX
           endif
           if(BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
             call NEWAB(ibin,upot,ken,beta_out,al_out,s_out,g_all_out)
           endif
         enddo
         call POTENTIAL(ndim,outx,grad,esys2,RMS,.false.,EVAP)
         write(*,*) "Check again upot=",esys2,"=?=",upot
         if(dabs(esys2-upot).gt.0.0000001D0)then
           write(*,*) "ERROR!!! ENERGY NOT CORRECT!!"
           write(*,*) "YOUR COMPONENT MAY BE WRONG!!!"
         endif
         if (BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then         
           open(56,file=name6//"_static.inf",status="replace")
           open(57,file=name6//"_static.emc",status="replace")
      write(56,*) "   BIN  BETA      ALPHA     S       ENERGY     EMC"
           do J5=1,ibin
             EP=EMAX-dble(J5-1)*DE
             EMC=beta_out(J5)*EP+al_out(J5)
             write(56,*) J5,beta_out(J5),al_out(J5),s_out(J5),EP,EMC
             write(57,*) EP,EMC,dexp(-EMC/temp),beta_out(J5),al_out(J5),
     &s_out(J5)
C             write(*,138) J5,beta_out(J5),al_out(J5),s_out(J5),EP,np(1)
           enddo
          write(57,*) "# ----------------------------------------------"
        write(57,*) "# EP      EMC       P       BETA       AL        S"
           write(57,*) "# MIN POT=",upot,",NP=",np(1)
           write(56,*) "Coordinate:"
           if(dim3.eq.1)then
             do J5=1,ndim
               write(56,*) J5,outx(J5)
             enddo
           else
             do J5=1,ndim/3
               J6=J5*3
               write(56,139) J5,outx(J6-2),outx(J6-1),outx(J6)
             enddo
           endif
           write(56,*)
C           call NEWAB(ibin,upot,ken,beta,al,s_out)
C           write(*,*) "AFTER NEWAB"
C           do J5=1,ibin
C             EP=EMAX-dble(J5-1)*DE
C             write(*,138) J5,beta(J5),al(J5),s_out(J5),EP,np(1)
C           enddo 
           write(56,*) "--------------------------------------"         
           write(56,*) "Global Minimum=",upot
           close(56)  
           close(57)
         endif
         if(recover_rec.eq.1.and.(upot.lt.pot_total(np(1))))then
           open(58,file=name6//"_recover.rec",status="replace")
           write(58,*) J0+J7+1
           do J6=1,ens_num
             do J5=1,ndim
               if(J6.eq.1)then
                 write(58,137) upot,outx(J5)
               else
                 write(58,137) pot_total(np(J6-1)),(xini(np(J6-1),J5))
               endif
             enddo
           enddo
           do J6=1,ens_num
             write(58,145) temp_total(np(J6))
           enddo
           if(BH_SW1.eq.1.or.BH_SW2.eq.1.or.BH_SW3.eq.1.or.
     &BH_SW1.eq.4)then
             do J6=1,ens_num
               if(J6.eq.1)then
                 do J5=1,ibin
           write(58,141) beta0_total(np(J6)),s_out(J5)
                 enddo
               else
                 do J5=1,ibin
           write(58,141) beta0_total(np(J6)),s_total((np(J6)-1)*ibin+J5)
                 enddo
               endif
             enddo
             if(ERANGE_SW.ne.0)then
               do J6=1,ens_num
                 if(J6.eq.1)then
                   write(58,141) EMIN,TDE
                 else
                   write(58,141) EMIN_total(np(J6)),tde_total(np(J6))
                 endif
               enddo
             endif
           endif
           close(58)
         endif
         open(59,file=name6//".rec",status="replace")
         write(59,136) ndim2/3,upot,J0+MULMAX
         write(59,*)
         if (dim3.eq.0) then
           do J5=1,ndim,3
             if((type.eq.2.or.type.eq.9).and.(J5.gt.(A*3)))then
               write(59,131) kindb,outx(J5),outx(J5+1),outx(J5+2)
             else
               write(59,131) kinda,outx(J5),outx(J5+1),outx(J5+2)
             endif
           enddo
         else
           do J5=1,ndim
             write(59,134) kinda,outx(J5)
           enddo
         endif
         close(59)
         if (pdb_rec.eq.1)then
           j1=0
           write(61,143) "HEADER    ",pot_total(np(1)),J0+MULMAX
           if (dim3.eq.1) then
             do J5=1,ndim
               xinib(J5)=xini(np(1),J5)
             enddo
             call dih_gen(ndim,xinib,x_dih)
           endif
           do J5=1,ndim2,3
             j1=j1+1
             if (dim3.eq.0)then
               x_dih(J5)=xini(np(1),J5)
               x_dih(J5+1)=xini(np(1),J5+1)
               x_dih(J5+2)=xini(np(1),J5+2)
             endif
             if((type.eq.2.or.type.eq.9).and.(J5.gt.(A*3)))then
         write(61,142)'ATOM', J1,kindb, 'UNK',J1,(x_dih(J5+J6),J6=0,2)
             else
         write(61,142)'ATOM', J1,kinda, 'UNK',J1,(x_dih(J5+J6),J6=0,2)
             endif
           enddo
           write(61,144) "END"
         endif
         endtime=MPI_WTIME()*10.D0
         fgen=MPI_WTIME()*10.D0
         close(61)
        write(*,*)"Total CPU Time= ",dabs(endtime-initime)/10.D0,"(Sec)"
         write(51,*) ndim2/3," dimension",dabs(endtime-initime)/10.D0,
     &"(sec)"
         write(51,*) gen," generations ",per," sec per generation"
         close(51)
       endif
       call mpi_barrier(mpi_comm_world,mpierror)
C       if(myid.eq.0)write(*,*) "ready out of useseed",mpierror
       if (type.eq.1.or.type.eq.2) then
         call gp_io(ndim,upot,outx,0)
       endif
       if (type.eq.3.or.type.eq.4) then
         call protein_io(ndim,upot,outx,0)
       endif
       if (type.eq.5) then
         call eam_io(ndim,upot,outx,0)
C  INPUT MOMENT OUTPUT
       endif
       if (type.eq.7) then
         call chiman_io(ndim,upot,outx,0)
       endif
       if (type.eq.8) then
         call chima3n_io(ndim,upot,outx,0)
       endif
       if (type.eq.9) then
         call qscff_io(ndim,upot,outx,0)
       endif
       call mpi_finalize(mpierror)
       stop
       return
       end
                         
