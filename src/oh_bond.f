*
*  P.J.Hsu's dihedral angle protein  generator
*  This subroutine will transfer dihedral angles into x-y-z coordinates
*  The dimension redused from 3n to n-9;For positive values of dih 
*  ,the torsion is right hand
        subroutine oh_bond(eff_co,dih,enr,evap)
        integer eff_co,j1,j2,j3,j4,j5,neu_num,ch_st(1000),ch_end(1000)
        integer jmax
        real*8 dih(eff_co),x(eff_co*3+9),vt(2000),vtr(2000),enr,vmax
        real*8 FCOIL,GCOIL,THETA1,THETA2,pei,dist,enr_exvol,enr_dih
        real*8 mag_b,mag_c,mag_d,d11,d12,d13,d21,d22,d23,c_x,c_y,c_z
        real*8 d_x,d_y,d_z,f_x,f_y,f_z,b_x,b_y,b_z,x_3,y_3,z_3
        real*8 rzero,o_x(eff_co),o_y(eff_co),o_z(eff_co)
        real*8 h_x(eff_co),h_y(eff_co),h_z(eff_co),f,EXC_D
        logical evap
        parameter(pei=3.14159D0)
        common /CHIMA1/ FCOIL,GCOIL,THETA1,EXC_D
        common /CHIMA2/ ch_st,ch_end
        common /CHIMA3/ neu_num
        common /dih1/ f
        common /vt1/ vt,vtr
        common /vt2/ jmax
        evap=.false.
        rzero=(2.D0**(1.D0/6.D0))*FCOIL
        THETA2=THETA1*2.D0
        vmax=-1.0D6
        enr=0.D0
        enr_dih=0.D0
        x(1)=0.D0
        x(2)=0.D0
        x(3)=0.D0
        x(4)=dcos(THETA2)*FCOIL
        x(5)=dsin(THETA2)*FCOIL
        x(6)=0.D0
        x(7)=(1.D0+dcos(THETA2))*FCOIL
        x(8)=dsin(THETA2)
        x(9)=0.D0
        do j1=1,eff_co+1
C          if (abs(dih(j1)).gt.pei)then
C            if(J1.NE.eff_co+1)then
C              evap=.true.      
C              enr=9999999.D0
C            endif
C          endif
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
          mag_d=dsqrt(d_x**2+d_y**2+d_z**2)
          if(j1.eq.eff_co+1) goto 42
          c_x=x(j3+4)-x(j3+1)
          c_y=x(j3+5)-x(j3+2)
          c_z=x(j3+6)-x(j3+3)
          mag_c=dsqrt(c_x**2+c_y**2+c_z**2)
          f_x=(d_y*c_z-d_z*c_y)/(mag_c*mag_d)
          f_y=(d_z*c_x-d_x*c_z)/(mag_c*mag_d)
          f_z=(d_x*d_y-d_y*c_x)/(mag_c*mag_d)
          x_3=x(j3+1)+x(j3+4)-x(j3-2)
          y_3=x(j3+2)+x(j3+5)-x(j3-1)
          z_3=x(j3+3)+x(j3+6)-x(j3)
          b_x=x_3-x(j3+4)-dcos(THETA2)*(x(j3+4)-x(j3+1))
          b_y=y_3-x(j3+5)-dcos(THETA2)*(x(j3+5)-x(j3+2))
          b_z=z_3-x(j3+6)-dcos(THETA2)*(x(j3+6)-x(j3+3))
          mag_b=dsqrt(b_x**2+b_y**2+b_z**2)
      x(j3+7)=x_3+mag_b*dsin(dih(j1)/f)*d_x/mag_d-(1.D0-
     &dcos(dih(j1)/f))*b_x
      x(j3+8)=y_3+mag_b*dsin(dih(j1)/f)*d_y/mag_d-(1.D0-
     &dcos(dih(j1)/f))*b_y
      x(j3+9)=z_3+mag_b*dsin(dih(j1)/f)*d_z/mag_d-(1.D0-
     &dcos(dih(j1)/f))*b_z
42       o_x(j1)=(2.D0/3.D0)*x(j3-2)+(1.D0/3.D0)*x(j3+1)-
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
C         write(*,*) "dih",j1,"=",dih(j1)         
       enddo
       do j1=1,eff_co+1
         if (j1.le.eff_co.and.abs(dih(j1)).gt.pei)then
           enr_dih=enr_dih+20*GCOIL*(dih(j1)-pei)**2.
         endif
         do j2=1,eff_co+1
           if(j1.eq.j2) goto 14
           
           do j5=1,neu_num
             if(j1.ge.ch_st(j5).and.j1.le.ch_end(j5))goto 14
             if(j2.ge.ch_st(j5).and.j2.le.ch_end(j5))goto 14
           enddo
           dist=(o_x(j1)-h_x(j2))**2+(o_y(j1)-h_y(j2))**2+
     &(o_z(j1)-h_z(j2))**2
           dist=(FCOIL/(dsqrt(dist)+rzero))**6
           vt(j1)=vt(j1)+4.D0*GCOIL*(dist*(dist-1.D0))
           if (vt(j1).gt.vmax)then
             vmax=vt(j1)
             jmax=j1
           endif
           enr=enr+4.D0*GCOIL*(dist*(dist-1.D0))
14       enddo
C           write(*,*) j1,"vt=",vt(j1),enr
12     enddo
       do j1=1,eff_co+3
         j3=j1*3
         do j2=j1+3,eff_co+3
           j4=j2*3
           dist=(x(j3-2)-x(j4-2))**2+(x(j3-1)-x(j4-1))**2+
     &(x(j3)-x(j4))**2
C          dist=1.D0/(dist**3.D0)
C          enr_exvol=enr_exvol+(dist-1)*dist
            if (dsqrt(dist).le.EXC_D)then
              enr_exvol=enr_exvol+(EXC_D/dsqrt(dist))**12
C             write(*,*) "stuck!!",FCOIL,dsqrt(dist)
            endif
         enddo
       enddo
        enr=enr+enr_exvol+enr_dih
        vtr(jmax)=vtr(jmax)+1.D0
       return
       end
            
