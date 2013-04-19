        subroutine chiman_io(ndim,upot,dih,io)
        include "mpif.h"
        integer atom_num,ndim,io,nproc,myid,t1,t2,t3,t4,j1,j2,j3
        integer ch_st(1000),ch_end(1000),neu_num,num
        real*8 FCOIL,GCOIL,THETA1,dih(1000),upot,outx(6000),EXC_D
        real*8 d11,d12,d13,d21,d22,d23,d_x,d_y,d_z,mag_d,pei
        real*8 o_x(1000),o_y(1000),o_z(1000),a,b,c,d,e,f,g,h
        real*8 h_x(1000),h_y(1000),h_z(1000),dih_step,factor
        parameter(pei=3.1415962D0)
        character name*6,seed_name*4,kinda*2,kindb*2,kindab*2,pdb_name*6
        common nproc,myid
        common /CHIMA1/ FCOIL,GCOIL,THETA1,EXC_D
        common /CHIMA2/ ch_st,ch_end
        common /CHIMA3/ neu_num
        common /CHIMA4/ dih_step,a,b,c,d,e,f,g,h
        common /name1/ name
        common /pdbname1/ pdb_name
        common /pdb2/ o_x,o_y,o_z,h_x,h_y,h_z
        common /seeding/ num
        common /dih1/ factor
        common /kinda/kinda
        common /kindb/kindb
        common /kindab/kindab
        factor=1.D0
        if (io.eq.1)then
          open(67,file="atom_num.dat",status="old")
          read(67,*) atom_num,num
          close(67)
          ndim=atom_num
          open(68,file="chima1.dat",status="old")
          read(68,*) FCOIL,GCOIL,THETA1,EXC_D,dih_step
          read(68,*) a,b,c,d,e,f,g,h
          THETA1=(pei-THETA1)/2.D0
          close(68)
          open(69,file="chima2.dat",status="old")
          neu_num=0
57        read(69,*,end=59)
          neu_num=neu_num+1
          goto 57
59        rewind(69)
          do j2=1,neu_num
            read(69,*) ch_st(j2),ch_end(j2)
          enddo
          close(69)
          if (myid.eq.0)then
      write(*,*)"<<Caculating CHEN AND IMAMURA ",atom_num," particles>>"
            if (ch_st(1).ne.0.or.ch_end(1).ne.0)then
      write(*,*) "@ Heterpolymer!"
      write(*,*) "@ The Neutral Start / End from:"
      write(*,*) "@ ",("(",ch_st(j2),ch_end(j2),")",j2=1,neu_num)
            else 
              write(*,*) "@ Homopolymer!"
            endif
          endif
          t1=ndim/1000
          t2=ndim/100
          t2=mod(t2,10)
          t3=ndim/10
          t3=mod(t3,10)
          t4=mod(ndim,10)
          name=char(t1+48)//char(t2+48)//char(t3+48)//char(t4+48)//"Cu"
          
        else
          if (myid.eq.0)then
1           format(A18,F25.13,A5)
2           format(A3,I4,F15.10,F15.10,F15.10)
3           format(A15,F25.13,A5)
4           format(A3,F15.10,F15.10,F15.10)
5           format(F15.10,F15.10,F15.10)      
            t1=ndim/1000
            t2=ndim/100
            t2=mod(t2,10)
            t3=ndim/10
            t3=mod(t3,10)
            t4=mod(ndim,10)
            seed_name=char(t1+48)//char(t2+48)//char(t3+48)//char(t4+48)
            open(14,file=name//".dat",status="unknown")
            open(15,file=name//".xyz",status="unknown")
            open(16,file=seed_name,status="unknown")
            write(14,1) "Minimum Potential=",upot," (ev)"
            write(14,3) "Binding Energy=",upot/real(ndim)," (ev)"
            write(*,*) "The coordinate of atoms :"
            write(14,*) "The coordinate of atoms :"
            write(15,*) ndim*3+5,upot," CHIMA DIH"
            write(15,*)
            k=0
            kinda="Cu"
            kindb="Cu"
            kindab="Cu"
            call dih_gen(ndim,dih,outx)
            do j1=1,ndim+3
              j3=j1*3
              write(15,4) kinda,outx(j3-2),outx(j3-1),outx(j3)
              write(14,2) kinda,j1,outx(j3-2),outx(j3-1),outx(j3)
              write(*,2) kinda,j1,outx(j3-2),outx(j3-1),outx(j3)
              write(16,5) outx(j3-2),outx(j3-1),outx(j3)
            enddo
            do j1=1,ndim+1
              j3=j1*3
              d11=outx(j3+4)-outx(j3+1)
              d12=outx(j3+5)-outx(j3+2)
              d13=outx(j3+6)-outx(j3+3)
              d21=outx(j3+1)-outx(j3-2)
              d22=outx(j3+2)-outx(j3-1)
              d23=outx(j3+3)-outx(j3)
              d_x=d12*d23-d13*d22
              d_y=d13*d21-d11*d23
              d_z=d11*d22-d12*d21
              mag_d=dsqrt(d_x**2.D0+d_y**2.D0+d_z**2.D0)
              o_x(j1)=(2.D0/3.D0)*outx(j3-2)+(1.D0/3.D0)*outx(j3+1)-
     &FCOIL*d_x/(2.D0*mag_d)
              o_y(j1)=(2.D0/3.D0)*outx(j3-1)+(1.D0/3.D0)*outx(j3+2)-
     &FCOIL*d_y/(2.D0*mag_d)
              o_z(j1)=(2.D0/3.D0)*outx(j3)+(1.D0/3.D0)*outx(j3+3)-
     &FCOIL*d_z/(2.D0*mag_d)
              h_x(j1)=(1.D0/3.D0)*outx(j3-2)+(2.D0/3.D0)*outx(j3+1)+
     &FCOIL*d_x/(2.D0*mag_d)
              h_y(j1)=(1.D0/3.D0)*outx(j3-1)+(2.D0/3.D0)*outx(j3+2)+
     &FCOIL*d_y/(2.D0*mag_d)
              h_z(j1)=(1.D0/3.D0)*outx(j3)+(2.D0/3.D0)*outx(j3+3)+
     &FCOIL*d_z/(2.D0*mag_d)
                                                    
              write(15,4) "H",h_x(j1),h_y(j1),h_z(j1)
              write(15,4) "O",o_x(j1),o_y(j1),o_z(j1)
            enddo
            write(*,*) "----------------------------------------------"
            write(*,1) "Minimun Potential=",upot," (ev)",name
            write(*,3) "Binding Energy=",upot/real(ndim)," (ev)"
            close(14)
            close(15)
            close(16)
            pdb_name=name
            call pdb((ndim*3+9),outx)
          endif
        endif
        end
