        subroutine dih_gen(eff_co,dih,x)
        integer eff_co,j1,j3
C        real*8 dih(eff_co/3),x(eff_co+9)
        real*8 dih(12),x(45)
        real*8 FCOIL,GCOIL,THETA1,THETA2,pei
        real*8 mag_b,mag_c,mag_d,d11,d12,d13,d21,d22,d23,c_x,c_y,c_z
        real*8 d_x,d_y,d_z,f_x,f_y,f_z,b_x,b_y,b_z,x_3,y_3,z_3
        parameter(pei=3.14159D0)
        common /CHIMA1/ FCOIL,GCOIL,THETA1
        THETA2=pei-THETA1
        x(1)=0.D0
        x(2)=0.D0
        x(3)=0.D0
        x(4)=FCOIL*dcos(THETA2)
        x(5)=FCOIL*dsin(THETA2)
        x(6)=0.D0
        x(7)=FCOIL*(1.D0+dcos(THETA2))
        x(8)=FCOIL*dsin(THETA2)
        x(9)=0.D0
        do j1=1,eff_co/3
          j3=j1*3
          write(*,*) dih(j1)
          c_x=x(j3+4)-x(j3+1)
          c_y=x(j3+5)-x(j3+2)
          c_z=x(j3+6)-x(j3+3)
          mag_c=dsqrt(c_x**2.D0+c_y**2.D0+c_z**2.D0)
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
          f_x=(d_y*c_z-d_z*c_y)/(mag_c*mag_d)
          f_y=(d_z*c_x-d_x*c_z)/(mag_c*mag_d)
          f_z=(d_x*d_y-d_y*c_x)/(mag_c*mag_d)
C         x_3=dcos(THETA2)*c_x+dsin(THETA2)*f_x/mag_d+x(j3+4)
C         y_3=dcos(THETA2)*c_y+dsin(THETA2)*f_y/mag_d+x(j3+5)
C         z_3=dcos(THETA2)*c_z+dsin(THETA2)*f_z/mag_d+x(j3+6)
          x_3=x(j3+1)+x(j3+4)-x(j3-2)
          y_3=x(j3+2)+x(j3+5)-x(j3-1)
          z_3=x(j3+3)+x(j3+6)-x(j3)
          b_x=x_3-x(j3+4)-dcos(THETA2)*(x(j3+4)-x(j3+1))
          b_y=y_3-x(j3+5)-dcos(THETA2)*(x(j3+5)-x(j3+2))
          b_z=z_3-x(j3+6)-dcos(THETA2)*(x(j3+6)-x(j3+3))
          mag_b=dsqrt(b_x**2.D0+b_y**2.D0+b_z**2.D0)
      x(j3+7)=x_3+mag_b*dsin(dih(j1))*d_x/mag_d-(1.D0-dcos(dih(j1)))*b_x
      x(j3+8)=y_3+mag_b*dsin(dih(j1))*d_y/mag_d-(1.D0-dcos(dih(j1)))*b_y
      x(j3+9)=z_3+mag_b*dsin(dih(j1))*d_z/mag_d-(1.D0-dcos(dih(j1)))*b_z
       enddo
       return
       end
            
