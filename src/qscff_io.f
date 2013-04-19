        subroutine qscff_io(ndim,upot,outx,io)
        implicit none
        include "mpif.h"
        integer A,A0,DA,ndim,io,j,k,l,t1,t2,t3,t4,atom_num,num,type
        integer nproc,myid,GEN_NUM
        real*8 upot,outx(6000)
        real*8 Dii,ni,mi,aii,ci
        real*8 Djj,nj,mj,ajj,cj
        real*8 Dij0,mij0,nij0,aij0
        character name*6,kinda*2,kindb*2,kindab*2,seed_name*4,A_name*4
        common nproc,myid
        common/name1/ name
        common/engtype/ type
        common/seeding/ num
        common/qscaa/ Dii,ci,mi,ni,aii
        common/qscbb/ Djj,cj,mj,nj,ajj
        common/qscab/ Dij0,mij0,nij0,aij0
        common/kinda/ kinda
        common/kindb/ kindb
        common/kindab/ kindab
        common/alloy1/ A
        if (io.eq.1) then
          open(41,file="atom_num.dat",status="old")
          read(41,*) atom_num,num
          close(41)
          ndim=atom_num*3
          open(10,file="qscff.dat",status="old")
          read(10,*) A,kindab
          read(10,*) kinda,Dii,ci,mi,ni,aii
          read(10,*) kindb,Djj,cj,mj,nj,ajj
          close(10)
          Dij0=DSQRT(Dii*Djj)
          nij0=(ni+nj)/2.D0
          mij0=(mi+mj)/2.D0
          aij0=(aii+ajj)/2.D0
*********************
*    For Alloy      *
*********************
          if (myid.eq.0) then
            write(*,*) "<<Caculating      ",kinda,"+",kindb,">>"
            write(*,*) "@ Particle number",A,"+",atom_num-A
            write(*,*) "@ Total particle number ",atom_num
********************END OF APPLICATION**********************
            t1=(ndim)/3000
            t2=(ndim)/300
            t2=mod(t2,10)
            t3=(ndim)/30
            t3=mod(t3,10)
            t4=mod(((ndim)/3),10)
         name=char(t1+48)//char(t2+48)//char(t3+48)//char(t4+48)//kindab
          endif
        else
          if (myid.eq.0)then
1           format(A18,F25.13,A5)
2           format(A3,I4,F15.10,F15.10,F15.10)
3           format(A15,F25.13,A5)
4           format(A3,F15.10,F15.10,F15.10)
5           format(F15.10,F15.10,F15.10)           
6           format(A6,A8,F25.13)
7           format(a4,I7, 2X,a2,2X,  a3, I6,4X,3(F8.3))
8           format(A3)      
            t1=(ndim)/3000
            t2=(ndim)/300
            t2=mod(t2,10)
            t3=(ndim)/30
            t3=mod(t3,10)
            t4=mod(((ndim)/3),10)
            seed_name=char(t1+48)//char(t2+48)//char(t3+48)//char(t4+48)
            t1=(atom_num-A)/1000
            t2=(atom_num-A)/100
            t2=mod(t2,10)
            t3=(atom_num-A)/10
            t3=mod(t3,10)
            t4=mod((atom_num-A),10)
            A_name=char(t1+48)//char(t2+48)//char(t3+48)//char(t4+48)
            open(14,file=name//A_name//".dat",status="unknown")
            open(15,file=name//A_name//".xyz",status="unknown")
            open(16,file=seed_name,status="unknown")
            open(17,file=name//A_name//".pdb",status="unknown")
            write(14,1) "Minimum Potential=",upot," (ev)"
            write(14,3) "Binding Energy=",upot/((ndim)/3)," (ev)"
            write(*,*) "The coordinate of atoms :"
            write(14,*) "The coordinate of atoms :"
            write(15,*) 
     &((ndim)/3),upot," Q-SC Force Field Potential,A=",A
            write(15,*)"pot=",upot
            write(17,6)'HEADER      ',name,upot
            k=0
            do j=1,ndim,3
              k=k+1
              if (k.le.A)then
                write(15,4) kinda,(outx(j+l),l=0,2)
                write(14,2) kinda,k,(outx(j+l),l=0,2)
                write(*,2) kinda,k,(outx(j+l),l=0,2)
        write(17,7)'ATOM', k,     kinda, 'UNK', k,(outx(j+l),l=0,2)
              else
                write(15,4) kindb,(outx(j+l),l=0,2)
                write(14,2) kindb,k,(outx(j+l),l=0,2)
                write(*,2) kindb,k,(outx(j+l),l=0,2)
        write(17,7)'ATOM', k,     kindb, 'UNK', k,(outx(j+l),l=0,2)
              endif
              write(16,5) (outx(j+l),l=0,2)
            enddo
            write(17,8) "END"
            write(*,*) "-----------------------------------------------"
            write(*,1) "Minimun Potential=",upot," (ev)"
            write(*,3) "Binding Energy=",upot/((ndim)/3)," (ev)"
            close(14)
            close(15)
            close(16)
            close(17)
          endif
        endif
        end
