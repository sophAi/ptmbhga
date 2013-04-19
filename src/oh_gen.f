        SUBROUTINE OH_GEN(NATOMS,X,O_X,O_Y,O_Z,H_X,H_Y,H_Z)
        INTEGER NATOMS,J1,J3
        REAL*8 X(30000),O_X(NATOMS-2),O_Y(NATOMS-2),O_Z(NATOMS-2)
        REAL*8 H_X(NATOMS-2),H_Y(NATOMS-2),H_Z(NATOMS-2),MAG_CRO
        REAL*8 D11,D12,D13,D21,D22,D23
        REAL*8 MAG_CRO_AB,SQR_ABC,CROSS_ABCX,CROSS_ABCY,CROSS_ABCZ
        DO J1=1,NATOMS-2
          J3=J1*3
          D11=X(J3+1)-X(J3-2)
          D12=X(J3+2)-X(J3-1)
          D13=X(J3+3)-X(J3)
          D21=X(J3+4)-X(J3+1)
          D22=X(J3+5)-X(J3+2)
          D23=X(J3+6)-X(J3+3)
          CROSS_ABCX=D12*D23-D13*D22
          CROSS_ABCY=D13*D21-D11*D23
          CROSS_ABCZ=D11*D22-D12*D21
          MAG_CRO_AB=CROSS_ABCX**2.+CROSS_ABCY**2.+CROSS_ABCZ**2.
          SQR_ABC=DSQRT(MAG_CRO_AB)
          O_X(J1)=(2.D0/3.D0)*X(J3-2)+(1.D0/3.D0)*X(J3+1)+
     &FCOIL*CROSS_ABCX/(2.D0*SQR_ABC)
          O_Y(J1)=(2.D0/3.D0)*X(J3-1)+(1.D0/3.D0)*X(J3+2)+
     &FCOIL*CROSS_ABCY/(2.D0*SQR_ABC)
          O_Z(J1)=(2.D0/3.D0)*X(J3)+(1.D0/3.D0)*X(J3+3)+
     &FCOIL*CROSS_ABCZ/(2.D0*SQR_ABC)
          H_X(J1)=(1.D0/3.D0)*X(J3-2)+(2.D0/3.D0)*X(J3+1)-
     &FCOIL*CROSS_ABCX/(2.D0*SQR_ABC)
          H_Y(J1)=(1.D0/3.D0)*X(J3-1)+(2.D0/3.D0)*X(J3+2)-
     &FCOIL*CROSS_ABCY/(2.D0*SQR_ABC)
          H_Z(J1)=(1.D0/3.D0)*X(J3)+(2.D0/3.D0)*X(J3+3)-
     &FCOIL*CROSS_ABCZ/(2.D0*SQR_ABC)
          WRITE(*,*) H_X(J1),H_Y(J2),H_Z(J3)
        ENDDO
        RETURN
        END

          
