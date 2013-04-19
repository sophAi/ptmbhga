        SUBROUTINE BOND(NATOMS,X,ENR)
        INTEGER NATOMS,J1,J3
        REAL*8 ENR,COS_THEB,THEB,THETA0,CUTOFF_ANGLE
        REAL*8 X(NATOMS*3),MAG_BA,MAG_BC,MAG_ABC
        REAL*8 ENR_TEMP,ENR_TEMP2,DIST1,DIST2
        COMMON /BOND1/ THETA0,CUTOFF_ANGLE
        THETA0=1.8326D0
        CUTOFF_ANGLE=1.D0
        ENR_TEMP=0.D0
        ENR_TEMP2=0.D0
        ENR=0.D0
        DO J1=1,NATOMS-2
          J3=J1*3
          COS_THEB=(X(J3)-X(J3+3))*(X(J3+6)-X(J3+3))
     &         +(X(J3-1)-X(J3+2))*(X(J3+5)-X(J3+2))
     &         +(X(J3-2)-X(J3+1))*(X(J3+4)-X(J3+1))
          MAG_BA=(X(J3)-X(J3+3))**2.+(X(J3-1)-X(J3+2))**2.+
     &           (X(J3-2)-X(J3+1))**2.
          MAG_BC=(X(J3+6)-X(J3+3))**2.+(X(J3+5)-X(J3+2))**2.
     &          +(X(J3+4)-X(J3+1))**2.
          MAG_ABC=DSQRT(MAG_BA*MAG_BC)
          THEB=DACOS(COS_THEB/MAG_ABC)
C         WRITE(*,*) THEB,COS_THEB
C         IF (DABS(THEB-THETA0).LE.CUTOFF_ANGLE)THEN
          ENR_TMEP=ENR_TEMP+20.D0*(THEB-THETA0)**2.
C         ELSE
C           ENR_TEMP=0.D0
C         ENDIF
C*****BOND DISTANCE AND ORDER CONSTRAIN*************
          DIST1=(X(J3+1)-X(J3-2))**2.+(X(J3+2)-X(J3-1))**2.+
     &(X(J3+3)-X(J3))**2.
          ENR_TEMP2=ENR_TEMP2+20.D0*(DSQRT(DIST1)-CUTOFF_ANGLE)**2.
          IF(J1.EQ.NATOMS-2)THEN
            DIST2=(X(J3+4)-X(J3+1))**2.+(X(J3+5)-X(J3+2))**2.+
     &(X(J3+6)-X(J3+3))**2.
            ENR_TEMP2=ENR_TEMP2+20.D0*(DSQRT(DIST2)-CUTOFF_ANGLE)**2.
          ENDIF
        ENDDO
        ENR=ENR_TEMP+ENR_TEMP2
        RETURN
        END
