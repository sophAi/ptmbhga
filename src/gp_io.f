        subroutine gp_io(ndim,upot,outx,io)
        implicit none
        include "mpif.h"
        integer A,A0,DA,ndim,io,j,k,l,t1,t2,t3,t4,atom_num,num,type
        integer nproc,myid,GEN_NUM
        real*8 upot,outx(6000)
        real*8 eshlon,eta,p,q,rzero,rzero_max
        real*8 eshlon_a,eta_a,p_a,q_a,rzero_a
        real*8 eshlon_b,eta_b,p_b,q_b,rzero_b
        real*8 eshlon_ab,eta_ab,p_ab,q_ab,rzero_ab
        character name*6,kinda*2,kindb*2,kindab*2,seed_name*4
        common nproc,myid
        common/name1/ name
        common/engtype/ type
        common/seeding/ num
        common/parameter/ eshlon,eta,p,q,rzero
        common/parametera/ eshlon_a,eta_a,p_a,q_a,rzero_a
        common/parameterb/ eshlon_b,eta_b,p_b,q_b,rzero_b
        common/parameterab/ eshlon_ab,eta_ab,p_ab,q_ab,rzero_ab
        common/kinda/ kinda
        common/kindb/ kindb
        common/kindab/ kindab
        common/alloy1/ A
        common/alloy2/ rzero_max
        if (io.eq.1) then
          open(41,file="atom_num.dat",status="old")
          read(41,*) atom_num,num
          close(41)
          ndim=atom_num*3
          open(10,file="parameter1.dat",status="old")
          read(10,*) kinda,eshlon_a,eta_a,p_a,q_a,rzero_a
          close(10)
          eshlon=eshlon_a
          eta=eta_a
          p=p_a
          q=q_a
          rzero=rzero_a
*********************
*    For Alloy      *
*********************
          open(11,file="alloy.dat",status="old")
          read(11,*) A
          close(11)
          open(12,file="parameter2.dat",status="old")
          read(12,*) kindb,eshlon_b,eta_b,p_b,q_b,rzero_b
          close(12)
          open(13,file="parameter3.dat",status="old")
          read(13,*) kindab,eshlon_ab,eta_ab,p_ab,q_ab,rzero_ab
          close(13)
          if (type.eq.1)then
            rzero_max=1.D0
          endif
          if (type.eq.2)then
        if(rzero_a.ge.rzero_b.and.rzero_a.ge.rzero_ab) rzero_max=rzero_a
        if(rzero_b.ge.rzero_a.and.rzero_b.ge.rzero_ab) rzero_max=rzero_b
      if(rzero_ab.ge.rzero_a.and.rzero_ab.ge.rzero_b) rzero_max=rzero_ab
          endif
          if (myid.eq.0) then
            if (type.eq.1)then
             write(*,*)"<<Caculating ",kinda," ",atom_num," particles>>"
            endif
            if (type.eq.2)then
              write(*,*) "<<Caculating      ",kinda,"+",kindb,">>"
              write(*,*) "@ Particle number",A,"+",atom_num-A
              write(*,*) "@ Total particle number ",atom_num
            endif
********************END OF APPLICATION**********************
            t1=int((ndim)/3000)
            t2=int((ndim)/300)
            t2=int(mod(t2,10))
            t3=int((ndim)/30)
            t3=int(mod(t3,10))
            t4=int(mod(((ndim)/3),10))
            if (type.eq.1)then
          name=char(t1+48)//char(t2+48)//char(t3+48)//char(t4+48)//kinda
            endif
            if (type.eq.2)then
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
            t1=int((ndim)/3000)
            t2=int((ndim)/300)
            t2=int(mod(t2,10))
            t3=int((ndim)/30)
            t3=int(mod(t3,10))
            t4=mod(((ndim)/3),10)
            seed_name=char(t1+48)//char(t2+48)//char(t3+48)//char(t4+48)
            open(14,file=name//".dat",status="unknown")
            open(15,file=name//".xyz",status="unknown")
            open(16,file=seed_name,status="unknown")
            open(17,file=name//".pdb",status="unknown")
            write(14,1) "Minimum Potential=",upot," (ev)"
            write(14,3) "Binding Energy=",upot/dble(ndim/3)," (ev)"
            write(*,*) "The coordinate of atoms :"
            write(14,*) "The coordinate of atoms :"
            if (type.eq.1)then
              write(15,*) int(ndim/3),upot," Gupta Potential Pure"
            else
              write(15,*) 
     &int(ndim/3),upot," Gupta Potential Alloy,A=",A
            endif
            write(15,*)"pot=",upot
            write(17,6)'HEADER      ',name,upot
            k=0
            do j=1,ndim,3
              k=k+1
              if (type.eq.1.or.(type.eq.2.and.(k.le.A)))then
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
            write(*,3) "Binding Energy=",upot/dble(ndim/3)," (ev)"
            close(14)
            close(15)
            close(16)
            close(17)
          endif
        endif
        end
