        program test_dih
        integer eff_co,i,j,j1,j3,turn,fl_n,neu_num
        real*8 dih(17),x(6000),FCOIL,GCOIL,THETA1,angled,pei,theta0
        real*8 d11,d12,d13,d21,d22,d23,mag_d,o_x(1000),o_y(1000),fluct
        real*8 o_z(1000),h_x(1000),h_y(1000),h_z(1000),d_x,d_y,d_z
        real*8 vt(1000),PROTOH_FUNC,EREAL,a1,a2,a3,a4,a5,a6,a7,a8
        real*8 ch_st(1000),ch_end(1000),EXC_D,f
        character pdb_name*6
        logical EVAP
        common/CHIMA1/ FCOIL,GCOIL,THETA1,EXC_D
        common/CHIMA2/ ch_st,ch_end
        common/CHIMA3/ neu_num
        common /pdbname1/ pdb_name
        common /dih1/ f
        parameter(pei=3.1415926D0)
        pdb_name="testdh"
        THETA1=(pei-(pei*105.D0/180.D0))/2.D0
        FCOIL=1.5D0
        EXC_D=1.D0
        GCOIL=1.D0
        EREAL=0.D0
        f=1.D0
        open(15,file="testdih.dat",status="old")
        read(15,*) eff_co,angled,turn,fluct,fl_n
        read(15,*) a1,a2,a3,a4,a5,a6,a7,a8
        eff_co=eff_co-3
        close(15)
        theta0=angled*pei/180.D0
        write(*,*) "The dih angles:",angled,"<degree>  =",theta0,"<rad>"
        do i=1,eff_co
          if(mod(i,fl_n).eq.0)then
            dih(i)=fluct*pei/180.D0
          else
            dih(i)=theta0
          endif
        enddo
        if(turn.ne.0)then
          dih(turn-3)=a1*pei/180.D0
C       determine beta width(slope) of ladder shape;side 1
          dih(turn-2)=a2*pei/180.D0        
C       determine beta turning point angle:side1
          dih(turn-1)=a3*pei/180.D0
C       determine two side's angle /\
          dih(turn)=a4*pei/180.D0
C       determine beta turning point angle;side2
          dih(turn+1)=a5*pei/180.D0
          dih(turn+2)=a6*pei/180.D0
C       determine beta width(slope) of ladder shape;side2
          dih(turn+3)=a7*pei/180.D0
          dih(turn+4)=a8*pei/180.D0
        endif
        call dih_gen(eff_co,dih,x)
        EREAL=PROTOH_FUNC(eff_co,dih,EVAP)
        open(1,file="test.xyz",status="replace")
        write(1,*) eff_co+3+((eff_co+1)*2),EREAL
        write(1,*) 
        do i=1,eff_co
      write(*,*) "dih",i,"=",dih(i)*180.D0/pei,"<deg>=",
     &dih(i),"<rad>"
        enddo
        do i=1,eff_co+3
          j=i*3
          write(1,*) "Rb",x(j-2),x(j-1),x(j)
          write(*,*) "Rb",i,x(j-2),x(j-1),x(j)
        enddo
        do j1=1,eff_co+1
          j3=j1*3
          d11=x(j3+4)-x(j3+1)
          d12=x(j3+5)-x(j3+2)
          d13=x(j3+6)-x(j3+3)
          d21=x(j3+1)-x(j3-2)
          d22=x(j3+2)-x(j3-1)
          d23=x(j3+3)-x(j3)
          d_x=d12*d23-d13*d22
          d_y=d13*d21-d11*d23
          d_z=d11*d22-d12*d21
          mag_d=dsqrt(d_x**2.D0+d_y**2.D0+d_z**2.D0)
          o_x(j1)=(2.D0/3.D0)*x(j3-2)+(1.D0/3.D0)*x(j3+1)-
     &FCOIL*d_x/(2.D0*mag_d)
          o_y(j1)=(2.D0/3.D0)*x(j3-1)+(1.D0/3.D0)*x(j3+2)-
     &FCOIL*d_y/(2.D0*mag_d)
          o_z(j1)=(2.D0/3.D0)*x(j3)+(1.D0/3.D0)*x(j3+3)-
     &FCOIL*d_z/(2.D0*mag_d)
          h_x(j1)=(1.D0/3.D0)*x(j3-2)+(2.D0/3.D0)*x(j3+1)+
     &FCOIL*d_x/(2.D0*mag_d)
          h_y(j1)=(1.D0/3.D0)*x(j3-1)+(2.D0/3.D0)*x(j3+2)+
     &FCOIL*d_y/(2.D0*mag_d)
          h_z(j1)=(1.D0/3.D0)*x(j3)+(2.D0/3.D0)*x(j3+3)+
     &FCOIL*d_z/(2.D0*mag_d)
          write(1,*) "O",o_x(j1),o_y(j1),o_z(j1)
          write(1,*) "H",h_x(j1),h_y(j1),h_z(j1)
        enddo  
        write(*,*) "ENERGY=",EREAL,",Bond energy=",EREAL/eff_co
        close(1)
        call pdb((eff_co*3+9),x)
        stop
        end
