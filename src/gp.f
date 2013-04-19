**************************************************
*  Energy and Gradient for Many-Body Potential.  *
**************************************************
      SUBROUTINE GP(NDIM,X,ELJ,EVAP)
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4,J5,J6,NATOMS,SW,NDIM,JMAX
      DOUBLE PRECISION ZETA,EPSILON,P,Q,RZERO,DX,DIST_DX,GRAD(6000)
     &,ATT(NDIM),REP(NDIM),ATT_DUMMY,REP_DUMMY
     &,ELJ_REP,ELJ_ATT,VT(2000),RADIUS,VTR(2000)
     &,ELJ_REPI,ELJ_ATTI,NONE,GRAD_ATTM(2000,2000)
     &,X(NDIM), DIST,ELJ,GRAD_REPM(2000,2000)
     &,DUMMYX,DUMMYY,DUMMYZ,XMUL2,XMUL3,GRAD_ATT(NDIM)
     &,VMIN,VMAX,DMAX
      LOGICAL EVAP,GRADT
      COMMON /parameter/ EPSILON,ZETA,P,Q,RZERO
      COMMON /accuracy3/ RADIUS,NONE
      COMMON /cgsw/ SW
      COMMON /GRAD1/ GRAD
      COMMON /GRAD2/ GRADT
      COMMON /vt1/ VT,VTR
      COMMON /vt2/ JMAX
      COMMON /vt3/ VMIN,VMAX,DMAX
      NATOMS=NDIM/3
      EVAP=.FALSE.
      VMAX=-1.0D6
      VMIN=1.0D6
      DMAX=-1.0D0
      ELJ=0.0D0
      ELJ_REPI=0.0D0
      ELJ_ATTI=0.0D0
      DO J1=1,NATOMS
        ATT(J1)=0.0D0
        REP(J1)=0.0D0
        VT(J1)=0.0D0
        VTR(J1)=0.0D0
      ENDDO
      DO J1=1,NATOMS
        J3=3*J1
        DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
        IF (DSQRT(DIST).GT.DMAX) DMAX=DSQRT(DIST)
        IF (DIST.GT.RADIUS) THEN
          EVAP=.TRUE.
C          WRITE(*,*) "EVAP",RADIUS,">",DIST
          IF (SW.EQ.11)RETURN
            ELJ=ELJ+(DIST-RADIUS)**2
        ENDIF
        ELJ_REP=0.0D0
        ELJ_ATT=0.0D0
        REP(J1)=0.0D0
        ATT(J1)=0.0D0
        GRAD_REPM(J1,J1)=0.0D0
        GRAD_ATTM(J1,J1)=0.0D0
        DO J2=1,NATOMS
          J4=3*J2
          IF(J2.NE.J1)THEN
            DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+
     &(X(J3)-X(J4))**2
            DIST=DSQRT(DIST)
            ATT_DUMMY=DEXP(2.D0*Q*(1.D0-DIST))
            REP_DUMMY=DEXP(P*(1.D0-DIST))
            ATT(J1)=ATT(J1)+ATT_DUMMY
            REP(J1)=REP(J1)+REP_DUMMY
            ELJ_ATT=ELJ_ATT+ATT_DUMMY
            ELJ_REP=ELJ_REP+REP_DUMMY
            GRAD_REPM(J1,J2)=EPSILON*(-1.D0)*P*REP_DUMMY/DIST
            GRAD_ATTM(J1,J2)=ZETA*(-2.D0)*Q*ATT_DUMMY/DIST
          ENDIF
        ENDDO
        GRAD_ATT(J1)=DSQRT(ELJ_ATT)
        ELJ_REPI=ELJ_REPI+ELJ_REP
        ELJ_ATTI=ELJ_ATTI+GRAD_ATT(J1)
        VT(J1)=(EPSILON*REP(J1))-(ZETA*DSQRT(ATT(J1)))
        IF (VT(J1).GT.VMAX)THEN
          VMAX=VT(J1)
          JMAX=J1
        ENDIF
        IF (VT(J1).LT.VMIN) VMIN=VT(J1)
      ENDDO
      ELJ=(EPSILON*ELJ_REPI)-(ZETA*ELJ_ATTI)
      VTR(JMAX)=VTR(JMAX)+1.D0
      DO J1=1,NATOMS
        J3=J1*3
        DUMMYX=0.0D0
        DUMMYY=0.0D0
        DUMMYZ=0.0D0
        DO J2=1,NATOMS
          J4=J2*3
          XMUL2=GRAD_REPM(J1,J2)-(GRAD_ATTM(J1,J2)/(2.D0*GRAD_ATT(J1)))
          XMUL3=GRAD_REPM(J2,J1)-(GRAD_ATTM(J2,J1)/(2.D0*GRAD_ATT(J2)))
          DUMMYX=DUMMYX+(XMUL2+XMUL3)*(X(J3-2)-X(J4-2))
          DUMMYY=DUMMYY+(XMUL2+XMUL3)*(X(J3-1)-X(J4-1))
          DUMMYZ=DUMMYZ+(XMUL2+XMUL3)*(X(J3)  -X(J4))
        ENDDO
        DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
        IF (DIST.GT.RADIUS) THEN
          DUMMYX=DUMMYX+(DIST-RADIUS)*X(J3-2)
          DUMMYY=DUMMYY+(DIST-RADIUS)*X(J3-1)
          DUMMYZ=DUMMYZ+(DIST-RADIUS)*X(J3)
        ENDIF
        GRAD(J3-2)=DUMMYX
        GRAD(J3-1)=DUMMYY
        GRAD(J3)=DUMMYZ
      ENDDO
      RETURN
      END
