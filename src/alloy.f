*************************************************************************
*Energy and Gradient for Many-Body Alloy Potential By W.T.Lin And K.L.Wu*
*************************************************************************
      SUBROUTINE ALLOY(NDIM,X,ENR,EVAP)
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4,J5,J6,NATOMS,A,SW,NDIM,JMAX
      DOUBLE PRECISIONZETA,EPSILON,P,Q,RZERO
     &,ZETA1,EPSILON1,P1,Q1,RZERO1
     &,ZETA2,EPSILON2,P2,Q2,RZERO2
     &,ZETA3,EPSILON3,P3,Q3,RZERO3
     &,ATT(NDIM),REP(NDIM),ATT_DUMMY,REP_DUMMY
     &,ENR_REP,ENR_ATT,VT(2000),RADIUS,VTR(2000)
     &,ENR_REPI,ENR_ATTI,NONE,VMIN,VMAX,DMAX
     &,X(NDIM), DIST,ENR,GRAD(6000)
     &,GRAD_ATTM(2000,2000),GRAD_REPM(2000,2000)
     &,DUMMYX,DUMMYY,DUMMYZ,XMUL2,XMUL3,GRAD_ATT(NDIM)
      LOGICAL EVAP,GRADT
      COMMON /parametera/ EPSILON1,ZETA1,P1,Q1,RZERO1
      COMMON /parameterb/ EPSILON2,ZETA2,P2,Q2,RZERO2
      COMMON /parameterab/ EPSILON3,ZETA3,P3,Q3,RZERO3
      COMMON /accuracy3/ RADIUS,NONE
      COMMON /alloy1/ A
      COMMON /cgsw/ SW
      COMMON /GRAD1/ GRAD
      COMMON /GRAD2/ GRADT
      COMMON /vt1/ VT,VTR
      COMMON /vt2/ JMAX
      COMMON /vt3/ VMIN,VMAX,DMAX
      NATOMS=NDIM/3
      EVAP=.FALSE.
      VMAX=-1.D06
      VMIN=1.0D6
      DMAX=-1.0D0
      ENR=0.0D0
      ENR_REPI=0.0D0
      ENR_ATTI=0.0D0
      DO J1=1,NATOMS
        J3=3*J1
        DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
        IF (DSQRT(DIST).GT.DMAX) DMAX=DSQRT(DIST)
        IF (DIST.GT.RADIUS) THEN
          EVAP=.TRUE.
C          WRITE(*,*) "EVAP",DIST,RADIUS
          IF(SW.EQ.11)RETURN
          ENR=ENR+((DIST-RADIUS)/RADIUS)**2
        ENDIF
        ENR_REP=0.0D0
        ENR_ATT=0.0D0
        REP(J1)=0.0D0
        ATT(J1)=0.0D0
        GRAD_ATT(J1)=0.D0
        GRAD_REPM(J1,J1)=0.D0
        GRAD_ATTM(J1,J1)=0.D0
        DO J2=1,NATOMS
          J4=3*J2
          IF(J2.NE.J1)THEN
            DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+
     &(X(J3)-X(J4))**2
            DIST=DSQRT(DIST)
            IF ((J1.LE.A) .AND. (J2.LE.A) ) THEN
               P=P1
               Q=Q1
               EPSILON=EPSILON1
               ZETA=ZETA1
               RZERO=RZERO1
            ELSE IF ( (J1.GT.A ) .AND .(J2.GT.A ) ) THEN
               P=P2
               Q=Q2
               EPSILON=EPSILON2
               ZETA=ZETA2
               RZERO=RZERO2
            ELSE
               P=P3
               Q=Q3
               EPSILON=EPSILON3
               ZETA=ZETA3
               RZERO=RZERO3
            ENDIF
            ATT_DUMMY=DEXP(2.D0*Q*(1.D0-(DIST/RZERO)))
            REP_DUMMY=DEXP(P*(1.D0-(DIST/RZERO)))
            ATT(J1)=ATT(J1)+(ZETA**2)*ATT_DUMMY
            REP(J1)=REP(J1)+EPSILON*REP_DUMMY
            ENR_ATT=ENR_ATT+(ZETA**2)*ATT_DUMMY
            ENR_REP=ENR_REP+EPSILON*REP_DUMMY
            GRAD_REPM(J1,J2)=EPSILON*(-1.D0)*P*REP_DUMMY/(DIST*RZERO)
            GRAD_ATTM(J1,J2)=(ZETA**2.D0)*
     &(-1.D0)*Q*ATT_DUMMY/(DIST*RZERO)
          ENDIF
        ENDDO
        GRAD_ATT(J1)=DSQRT(ENR_ATT)
        ENR_REPI=ENR_REPI+ENR_REP
        ENR_ATTI=ENR_ATTI+GRAD_ATT(J1)
        VT(J1)=REP(J1)-DSQRT(ATT(J1))
        IF (VT(J1).GT.VMAX)THEN
          VMAX=VT(J1)
          JMAX=J1
        ENDIF
        IF (VT(J1).LT.VMIN) VMIN=VT(J1)
      ENDDO
      ENR=(ENR_REPI)-(ENR_ATTI)
      VTR(JMAX)=VTR(JMAX)+1.D0
      DO J1=1,NATOMS
        J3=J1*3
        DUMMYX=0.0D0
        DUMMYY=0.0D0
        DUMMYZ=0.0D0
        DO J2=1,NATOMS
          J4=J2*3
          XMUL2=GRAD_REPM(J1,J2)-(GRAD_ATTM(J1,J2)/GRAD_ATT(J1))
          XMUL3=GRAD_REPM(J2,J1)-(GRAD_ATTM(J2,J1)/GRAD_ATT(J2))
          DUMMYX=DUMMYX+(XMUL2+XMUL3)*(X(J3-2)-X(J4-2))
          DUMMYY=DUMMYY+(XMUL2+XMUL3)*(X(J3-1)-X(J4-1))
          DUMMYZ=DUMMYZ+(XMUL2+XMUL3)*(X(J3)  -X(J4))
        ENDDO   
        DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
        IF (DIST.GT.RADIUS) THEN
          DUMMYX=DUMMYX+((DIST-RADIUS)/RADIUS)*X(J3-2)
          DUMMYY=DUMMYY+((DIST-RADIUS)/RADIUS)*X(J3-1) 
          DUMMYZ=DUMMYZ+((DIST-RADIUS)/RADIUS)*X(J3)
        ENDIF  
        GRAD(J3-2)=DUMMYX
        GRAD(J3-1)=DUMMYY
        GRAD(J3)=DUMMYZ
      ENDDO
      RETURN
      END
