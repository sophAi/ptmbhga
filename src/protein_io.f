        subroutine protein_io(ndim,upot,outx,io)
        include "mpif.h"
        integer chk_dim,I_TYPE(6000),ndim,io,ndim_test,i
        integer type
        integer t1,t2,t3,t4,j1,j2,nproc,myid,num,residues_num(6000)
        real*8 upot,outx(6000),ECOIL
        real*8 ACOIL,BCOIL,CCOIL,DCOIL,THETA0,CUTOFF_ANGLE,CUTOFF
        real*8 R3,A0,A1,A2,A3,A4,A5,B1,B2
        character name*6,load(6000)*3,kindpr*2,load2*2,seed_name*4
        common nproc,myid
        common /seeding/ num
        common /name1/ name
        common /kinda/ kindpr
        common /engtype/ type
        common /INTER_TYPE1/ I_TYPE
        common /LONG_RANGE1/ CUTOFF,ECOIL
        common /BOND1/ THETA0,CUTOFF_ANGLE
        common /DIH1/ ACOIL,BCOIL,CCOIL,DCOIL,R3,A0,A1,A2,A3,A4,A5,B1,B2
        if(io.eq.1)then
          open(14,file="atom_num.dat",status="old")
          read(14,*) atom_num,num
          close(14)
          ndim=atom_num*3
          ndim_test=0
          chk_dim=0
          open(15,file="protein1.dat",status="old")
          read(15,*)
     &CUTOFF,THETA0,CUTOFF_ANGLE,EPS,EPS2,EPS3,EPS4,ACOIL,BCOIL,CCOIL,
     &DCOIL,ECOIL
          close(15)
          R3=DSQRT(3.D0)
          A3=(4.*EPS-4.*EPS1+22.*EPS2)/9.
          A5=4.*EPS2-2.*A3
          A1=(-1.)*A3/8.-5.*EPS2/4.
          A2=2.*EPS4/3.+13.*EPS2/3.-4.*EPS1/3.-3.*A3/2.
          A4=(-2.)*A2
          B1=(-10.)*EPS4/(R3*3.)
          B2=8.*EPS4/(R3*3.)
          A0=(-1.)*(3.*R3*B1/8.+9.*R3*B2/32.+
     &A1/2.+A2/4.+A3/8.+A4/16.+A5/32.)
          open(16,file="protein2.dat",status="old")
11        read(16,*,end=22)
          ndim_test=ndim_test+1
          goto 11
22        rewind(16)
          do i=1,ndim_test
            read(16,*) load(i),residues_num(i)
            do j=1,residues_num(i)
              chk_dim=chk_dim+1
              if (load(i).eq."B")then
                I_TYPE(chk_dim)=1
              endif
              if (load(i).eq."L")then
                I_TYPE(chk_dim)=2
              endif
              if (load(i).eq."N")then
                I_TYPE(chk_dim)=3
              endif
              if (load(i).eq."LB")then
                I_TYPE(chk_dim)=2
                chk_dim=chk_dim+1
                I_TYPE(chk_dim)=1
              endif
              if (load(i).eq."BBL")then
                I_TYPE(chk_dim)=1
                chk_dim=chk_dim+1
                I_TYPE(chk_dim)=1
                chk_dim=chkdim+1
                I_TYPE(chk_dim)=2
              endif
              if (load(i).eq."LBB")then
                I_TYPE(chk_dim)=2
                chk_dim=chk_dim+1
                I_TYPE(chk_dim)=1
                chk_dim=chk_dim+1
                I_TYPE(chk_dim)=1
              endif
              if (load(i).eq."BLB")then
                I_TYPE(chk_dim)=1
                chk_dim=chk_dim+1
                I_TYPE(chk_dim)=2
                chk_dim=chk_dim+1
                I_TYPE(chk_dim)=1
              endif
            enddo
          enddo
          close(16)
          if(myid.eq.0)then
            if(chk_dim.ne.(ndim/3))then
              write(*,*) "Particles and residues are not equal!!"
              write(*,*) "Particle number=",ndim/3
              write(*,*) "Residue number=",chk_dim
              write(*,*) "Please run config again!!,program stop!!"
              return
            endif
          endif
          if(type.eq.3)then
            kindpr="Be"
            write(*,*) "<<Doing Beta Protein Folding!>>"
          endif
          if(type.eq.4)then
            kindpr="Al"
            write(*,*) "<<Doing Alpha Protein Folding!>>"
          endif
          if(type.eq.6.or.type.eq.7.or.type.eq.8)then
            kindpr="Cu"
            write(*,*) "<<Doing Chan's Universal Protein Folding!>>"
            write(*,*) 
     &"@ EPS1:",EPS,",EPS2:",EPS2,",EPS3:",EPS3,",EPS4:",EPS4
          endif
          write(*,*)
     &"@ A:",ACOIL,",B:",BCOIL,",C:",CCOIL,",D:",DCOIL,",E:",ECOIL
       write(*,*) "@ ",("(",load(i),residues_num(i),")",i=1,ndim_test)
          t1=ndim/3000
          t2=ndim/300
          t2=mod(t2,10)
          t3=ndim/30
          t3=mod(t3,10)
          t4=mod((ndim/3),10)
          name=
     &char(t1+48)//char(t2+48)//char(t3+48)//char(t4+48)//kindpr
        else
          if (myid.eq.0)then
1           format(A18,F25.13,A5)  
2           format(A3,I4,F15.10,F15.10,F15.10)
4           format(A3,F15.10,F15.10,F15.10)
3           format(A15,F25.13,A5)
5           format(F15.10,F15.10,F15.10)
            t1=ndim/3000
            t2=ndim/300
            t2=mod(t2,10)
            t3=ndim/30
            t3=mod(t3,10)
            t4=mod((ndim/3),10)
            seed_name=char(t1+48)//char(t2+48)//char(t3+48)//char(t4+48)
            open(18,file=name//".dat",status="unknown")
            open(19,file=name//".xyz",status="unknown")
            open(20,file=seed_name,status="unknown")
            write(*,1) "Minimun Potential=",upot," (ev)"
            write(*,3) "Binding Energy=",upot/(ndim/3)," (ev)"
            write(18,1) "Minimum Potential=",upot," (ev)"
            write(18,3) "Binding Energy=",upot/(ndim/3)," (ev)"
            write(18,*) "The coordinate of atoms :"
            write(19,*) (ndim/3),upot," Protein"
            write(19,*)"pot=",upot
            do j1=1,(ndim/3)
              j2=j1*3
              if (I_TYPE(j1).eq.1)then
                load2="Be"
              endif
              if (I_TYPE(j1).eq.2)then
                load2="Li"
              endif
              if (I_TYPE(j1).eq.3)then
                load2="Ni"
              endif 
              write(18,2) load2,j1,outx(j2-2),outx(j2-1),outx(j2)
              write(*,2) load2,j1,outx(j2-2),outx(j2-1),outx(j2)
              write(19,4) load2,outx(j2-2),outx(j2-1),outx(j2)
              write(20,5) outx(j2-2),outx(j2-1),outx(j2)
            enddo
            close(18)
            close(19)
          endif
        endif      
        end      
