**********************************************************
*  Energy and for Many-Body Alloy Potential By W.T.Lin.  *
**********************************************************
      SUBROUTINE ALLOY(NATOMS,X,ENR,VT,EVAP)
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4,J5,J6,NATOMS,A,SW
      REAL*8 ZETA,EPSILON,P,Q,RZERO
      REAL*8 ZETA1,EPSILON1,P1,Q1,RZERO1
      REAL*8 ZETA2,EPSILON2,P2,Q2,RZERO2
      REAL*8 ZETA3,EPSILON3,P3,Q3,RZERO3
      REAL*8 ATT(NATOMS),REP(NATOMS),ATT_DUMMY,REP_DUMMY
      REAL*8 ENR_REP,ENR_ATT,VT(NATOMS),RADIUS
      REAL*8 ENR_REPI,ENR_ATTI,NONE
      REAL*8 X(NATOMS*3), DIST,ENR
      LOGICAL EVAP
      COMMON /parametera/ EPSILON1,ZETA1,P1,Q1,RZERO1
      COMMON /parameterb/ EPSILON2,ZETA2,P2,Q2,RZERO2
      COMMON /parameterab/ EPSILON3,ZETA3,P3,Q3,RZERO3
      COMMON /accuracy3/ RADIUS,NONE
      COMMON /alloy1/ A
      COMMON /cgsw/ SW
      EVAP=.FALSE.
      ENR=0.0D0
      ENR_REPI=0.0D0
      ENR_ATTI=0.0D0
      DO J1=1,NATOMS
        ATT(J1)=0.0D0
        REP(J1)=0.0D0
        VT(J1)=0.0D0
      ENDDO
      DO J1=1,NATOMS
        J3=3*J1
        DIST=X(J3-2)**2.+X(J3-1)**2.+X(J3)**2.
        DIST=DSQRT(DIST)
         
C       IF (DIST.GT.RADIUS) THEN
C         EVAP=.TRUE.
C         IF(SW.EQ.11)RETURN
C         ENR=ENR+(DIST-RADIUS)**2.
C       ENDIF
        ENR_REP=0.0D0
        ENR_ATT=0.0D0
        REP(J1)=0.0D0
        ATT(J1)=0.0D0
        DO J2=1,NATOMS
          J4=3*J2
          IF(J2.NE.J1)THEN
            DIST=(X(J3-2)-X(J4-2))**2.+(X(J3-1)-X(J4-1))**2.+
     &(X(J3)-X(J4))**2.
            IF ((J1.LE.A) .AND. (J2.LE.A) ) THEN
               P=P1
               Q=Q1
               EPSILON=EPSILON1
               ZETA=ZETA1
            ELSE IF ( (J1.GT.A ) .AND .(J2.GT.A ) ) THEN
               P=P2
               Q=Q2
               EPSILON=EPSILON2
               ZETA=ZETA2
            ELSE
               P=P3
               Q=Q3
               EPSILON=EPSILON3
               ZETA=ZETA3
            ENDIF 
            ATT_DUMMY=DEXP(2.*Q*(1.-DSQRT(DIST)))
            REP_DUMMY=DEXP(P*(1.-DSQRT(DIST)))
            ATT(J1)=ATT(J1)+ZETA**2*ATT_DUMMY
            REP(J1)=REP(J1)+EPSILON*REP_DUMMY
            ENR_ATT=ENR_ATT+ZETA**2*ATT_DUMMY
            ENR_REP=ENR_REP+EPSILON*REP_DUMMY
          ENDIF
        ENDDO
        ENR_REPI=ENR_REPI+ENR_REP
        ENR_ATTI=ENR_ATTI+DSQRT(ENR_ATT)
        VT(J1)=REP(J1)-DSQRT(ATT(J1))
      ENDDO
      ENR=(ENR_REPI)-(ENR_ATTI)
      RETURN
      END
