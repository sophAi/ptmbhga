        subroutine eam_io(ndim,upot,outx,io)
        implicit none
        include "mpif.h"
        integer A,A0,DA,ndim,io,j,k,l,t1,t2,t3,t4,atom_num,num,type
        integer nproc,myid,GEN_NUM
        real*8 upot,outx(6000)
        real*8 eam_a1,eam_epsilon1,eam_C1,eam_n1,eam_m1
        real*8 eam_a2,eam_epsilon2,eam_C2,eam_n2,eam_m2
        real*8 eam_a3,eam_epsilon3,eam_C3,eam_n3,eam_m3
        character name*6,kinda*2,kindb*2,kindab*2,seed_name*4
        common nproc,myid
        common/name1/ name
        common/engtype/ type
        common/seeding/ num
        common/eam_parameter1/ eam_a1,eam_epsilon1,eam_C1,eam_n1,eam_m1
        common/kinda/ kinda
        common/kindb/ kindb
        common/kindab/ kindab
        common/alloy1/ A
        if (io.eq.1) then
          open(41,file="atom_num.dat",status="old")
          read(41,*) atom_num,num
          close(41)
          open(42,file="alloy.dat",status="old")
          read(42,*) A
          close(42)
          ndim=atom_num*3
          open(50,file="eam_aa.dat",status="old")
          read(50,*) kinda,eam_a1,eam_epsilon1,eam_C1,eam_n1,eam_m1
          if(A.ne.atom_num)then
            open(51,file="eam_bb.dat",status="old")
            read(51,*) kindb,eam_a2,eam_epsilon2,eam_C2,eam_n2,eam_m2
            open(52,file="eam_ab.dat",status="old")
            read(52,*) kindab,eam_a3,eam_epsilon3,eam_C3,eam_n3,eam_m3
            close(52)
            close(51)
          endif
          close(50)
*********************
*    For Alloy      *
*********************
          if (myid.eq.0) then
            if (type.eq.5)then
         write(*,*)"<<Caculating EAM,",kinda,"=",atom_num," particles>>"
            endif
            if (A.ne.atom_num)then
              write(*,*) "<<Caculating alloy,",kinda,"+",kindb,">>"
              write(*,*) "@ Particle number",A,"+",atom_num-A
              write(*,*) "@ Total particle number ",atom_num
            endif
********************END OF APPLICATION**********************
            t1=(ndim)/3000
            t2=(ndim)/300
            t2=mod(t2,10)
            t3=(ndim)/30
            t3=mod(t3,10)
            t4=mod(((ndim)/3),10)
            if (type.eq.5)then
          name=char(t1+48)//char(t2+48)//char(t3+48)//char(t4+48)//kinda
            endif
            if (A.ne.atom_num)then
         name=char(t1+48)//char(t2+48)//char(t3+48)//char(t4+48)//kindab
            endif
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
            open(14,file=name//".dat",status="unknown")
            open(15,file=name//".xyz",status="unknown")
            open(16,file=seed_name,status="unknown")
            open(17,file=name//".pdb",status="unknown")
            write(14,1) "Minimum Potential=",upot," (ev)"
            write(14,3) "Binding Energy=",upot/((ndim)/3)," (ev)"
            write(*,*) "The coordinate of atoms :"
            write(14,*) "The coordinate of atoms :"
            write(15,*) ((ndim)/3),upot," EAM Potential,A=",A
            write(15,*)"pot=",upot
            write(17,6)'HEADER      ',name,upot
            k=0
            write(*,*) "NAME=",name,kinda
            do j=1,ndim,3
              k=k+1
              if (type.eq.1.or.type.eq.5.or.(A.ne.atom_num.and.(k.le.A)
     &))then
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
          endif
            close(14)
            close(15)
            close(16)
            close(17)
        endif
        end
