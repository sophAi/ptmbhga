        SUBROUTINE CHIMA(NATOMS,X,ENR,EVAP)
        INTEGER NATOMS,J1,J2,J3,J4,J5,CH_ST(1000),CH_END(1000),NEU_NUM
        INTEGER JMAX
        REAL*8 X(NATOMS*3),ENR_CHIMA,ENR_BOND,ENR,DIST,DIST_BOND,CUTOFF
        REAL*8 RADIUS,NONE,FCOIL,GCOIL,RZERO,THETA1
        REAL*8 D11,D12,D13,D21,D22,D23,CROSS_ABCX,CROSS_ABCY,CROSS_ABCZ
        REAL*8 O_X(NATOMS-2),O_Y(NATOMS-2),O_Z(NATOMS-2),MAG_CRO_AB
        REAL*8 H_X(NATOMS-2),H_Y(NATOMS-2),H_Z(NATOMS-2),SQR_ABC
        REAL*8 DIST1,DIST2,DIST3,ENR_DIST,ENR_EXVOL,GRAD(NATOMS*3)
        REAL*8 GRAD_BOND,GRAD_DIST1,GRAD_DIST2
        REAL*8 VMIN,VMAX,DMAX,VT(2000),VTR(2000),EXC_D
        LOGICAL EVAP
        COMMON /CHIMA1/ FCOIL,GCOIL,THETA1,EXC_D
        COMMON /CHIMA2/ CH_ST,CH_END
        COMMON /CHIMA3/ NEU_NUM
        COMMON /accuracy3/ RADIUS,NONE
        COMMON /vt1/ VT,VTR
        COMMON /vt2/ JMAX
        COMMON /vt3/ VMIN,VMAX,DMAX
C******************8
C       LO=FCOIL,epsilon=GCOIL
        EVAP=.FALSE.
        VMAX=-1.0D6
        VMIN=1.0D6
        DMAX=-1.0D0
        DO J3=1,NATOMS*3
          GRAD(J3)=0.D0
        ENDDO
        RZERO=(2.D0**(1./6.))*FCOIL
        CUTOFF=2.D0*DCOS(THETA1)*FCOIL
        ENR=0.D0
        ENR_CHIMA=0.D0
        ENR_BOND=0.D0
        ENR_DIST=0.D0
        ENR_EXVOL=0.D0
C        write(*,*) GCOIL,FCOIL,RZERO,CUTOFF
        DO J1=1,NATOMS-2
C          DO J5=1,NEU_NUM
Ci           IF (J1.GE.CH_ST(J5).AND.J1.LE.CH_END(J5)) GOTO 15
C          ENDDO
          VT(J1)=0.D0
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
C     BOND ANGLE
15        DIST_BOND=(X(J3+4)-X(J3-2))**2.+(X(J3+5)-X(J3-1))**2.+
     &(X(J3+6)-X(J3))**2.
          DIST_BOND=DSQRT(DIST_BOND)
          ENR_BOND=ENR_BOND+20.D0*GCOIL*(DIST_BOND-CUTOFF)**2.
C GRADIENT OF BOND ANGLE
          GRAD_BOND=40.D0*GCOIL*(DST_BOND-CUTOFF)/DIST_BOND
          GRAD(J3-2)=GRAD_BOND*(X(J3-2)-X(J3+4))
          GRAD(J3-1)=GRAD_BOND*(X(J3-1)-X(J3+5))
          GRAD(J3)=GRAD_BOND*(X(J3)-X(J3+6))
          GRAD(J3+4)=GRAD(J3+4)+GRAD_BOND*(X(J3+4)-X(J3-2))
          GRAD(J3+5)=GRAD(J3+5)+GRAD_BOND*(X(J3+5)-X(J3-1))
          GRAD(J3+6)=GRAD(J3+6)+GRAD_BOND*(X(J3+6)-X(J3))
C END OF BOND-ANGLE ANGLE
          DIST1=(X(J3+1)-X(J3-2))**2.+(X(J3+2)-X(J3-1))**2.+
     &(X(J3+3)-X(J3))**2.
          ENR_DIST=ENR_DIST+20.D0*(DSQRT(DIST1)-1.D0)**2.
          GRAD_DIST1=40.D0*GCOIL*(DSQRT(DIST1)-1.D0)/DSQRT(DIST1)
          GRAD(J3-2)=GRAD(J3-2)+GRAD_DIST1*(X(J3-2)-X(J3+1))
          GRAD(J3-1)=GRAD(J3-1)+GRAD_DIST1*(X(J3-1)-X(J3+2))
          GRAD(J3)=GRAD(J3)+GRAD_DIST1*(X(J3)-X(J3+3))
          GRAD(J3+1)=GRAD(J3+1)+GRAD_DIST1*(X(J3+1)-X(J3-2))
          GRAD(J3+2)=GRAD(J3+2)+GRAD_DIST1*(X(J3+2)-X(J3-1))
          GRAD(J3+3)=GRAD(J3+3)+GRAD_DIST1*(X(J3+3)-X(J3))
          IF(J1.EQ.NATOMS-2)THEN
            DIST2=(X(J3+4)-X(J3+1))**2.+(X(J3+5)-X(J3+2))**2.+
     &(X(J3+6)-X(J3+3))**2.
            ENR_DIST=ENR_DIST+20.D0*(DSQRT(DIST2)-1.D0)**2.
            GRAD_DIST2=40.D0*GCOIL*(DSQRT(DIST2)-1.D0)/DSQRT(DIST2)
            GRAD(J3+1)=GRAD(J3+1)+GRAD_DIST2*(X(J3+1)-X(J3+4))
            GRAD(J3+2)=GRAD(J3+2)+GRAD_DIST2*(X(J3+2)-X(J3+5))
            GRAD(J3+3)=GRAD(J3+3)+GRAD_DIST2*(X(J3+3)-X(J3+6))
            GRAD(J3+4)=GRAD(J3+4)+GRAD_DIST2*(X(J3+4)-X(J3+1))
            GRAD(J3+5)=GRAD(J3+5)+GRAD_DIST2*(X(J3+5)-X(J3+2))
            GRAD(J3+6)=GRAD(J3+6)+GRAD_DIST2*(X(J3+6)-X(J3+3))
          ENDIF
        ENDDO
        DO J1=1,NATOMS-2
C          DO J5=1,NEU_NUM
C            IF(J1.GE.CH_ST(J5).AND.J1.LE.CH_END(J5))then
C              VT(J1)=0.D0
C              GOTO 12
C            ENDIF
C          ENDDO
          DO J2=1,NATOMS-2
C         DO J2=J1,NATOMS-2
            DO J5=1,NEU_NUM
              IF(J1.GE.CH_ST(J5).AND.J1.LE.CH_END(J5))GOTO 14
              IF(J1.EQ.J2)GOTO 14
              IF(J2.GE.CH_ST(J5).AND.J2.LE.CH_END(J5))GOTO 14
            ENDDO
            DIST=(O_X(J1)-H_X(J2))**2.+(O_Y(J1)-H_Y(J2))**2.+
     &(O_Z(J1)-H_Z(J2))**2.
            DIST=(FCOIL/(DSQRT(DIST)+RZERO))**6.
            VT(J1)=VT(J1)+4.D0*GCOIL*(DIST*(DIST-1.D0))
C            ENR_CHIMA=ENR_CHIMA+4.D0*GCOIL*(DIST*(DIST-1.D0))
14        ENDDO
          ENR_CHIMA=ENR_CHIMA+VT(J1)
          IF (VT(J1).GT.VMAX)THEN
            VMAX=VT(J1)
            JMAX=J1
          ENDIF
          IF (VT(J1).LT.VMIN) VMIN=VT(J1)                     
12      ENDDO
        VTR(JMAX)=VTR(JMAX)+1.D0
        DO J1=1,NATOMS
          J3=J1*3
          DO J2=J1+2,NATOMS
            J4=J2*3
            DIST3=(X(J3-2)-X(J4-2))**2.+(X(J3-1)-X(J4-1))**2.+
     &(X(J3)-X(J4))**2.
C           ENR_EXVOL=ENR_EXVOL+((FCOIL**2.)/DIST3)**6.-GCOIL
             IF (DSQRT(DIST3).LE.EXC_D) THEN
               ENR_EXVOL=ENR_EXVOL+(EXC_D/DSQRT(DIST3))**12.
             ENDIF
C           WRITE(*,*) "DIST3=",DSQRT(DIST3),FCOIL,ENR_EXVOL
          ENDDO
        ENDDO
C        write(*,*) ENR_BOND,ENR_CHIMA,ENR_EXVOL
        ENR=ENR_BOND+ENR_CHIMA+ENR_EXVOL+ENR_DIST
        RETURN
        END


          
