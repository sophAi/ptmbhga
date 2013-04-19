*******************************************************************************
*Energy and Gradient for Q-SC Many Body Force-field Potential By Po-Jen Hsu   *
*******************************************************************************
      SUBROUTINE QSCFF(NDIM,X,ENR,EVAP)
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4,J5,J6,NATOMS,A,SW,NDIM,JMAX
      REAL*8 Dii,ni,mi,aii,ci,Dij2
      REAL*8 Djj,nj,mj,ajj,cj
      REAL*8 Dij,nij,mij,aij,cij,Dij0,mij0,nij0,aij0
      REAL*8 ATT(NDIM),REP(NDIM),ATT_DUMMY,REP_DUMMY
      REAL*8 ENR_REP,ENR_ATT,VT(2000),RADIUS,VTR(2000)
      REAL*8 ENR_REPI,ENR_ATTI,NONE,VMIN,VMAX,DMAX
      REAL*8 X(NDIM), DIST,ENR,GRAD(6000)
      REAL*8 GRAD_ATTM(2000,2000),GRAD_REPM(2000,2000)
      REAL*8 DUMMYX,DUMMYY,DUMMYZ,XMUL2,XMUL3,GRAD_ATT(NDIM)
      LOGICAL EVAP,GRADT
      COMMON /qscaa/ Dii,ci,mi,ni,aii
      COMMON /qscbb/ Djj,cj,mj,nj,ajj
      COMMON /qscab/ Dij0,mij0,nij0,aij0
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
      IF (A.EQ.NATOMS)THEN
        RADIUS=2.D0*aii
      ELSE
        RADIUS=2.D0*aij0
      ENDIF
      DO J1=1,NATOMS
        J3=3*J1
        IF (J1.LE.A)THEN
          Cij=Ci
          Dij2=Dii
        ELSE
          Cij=Cj
          Dij2=Djj
        ENDIF
C        DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
C        IF (DSQRT(DIST).GT.DMAX) DMAX=DSQRT(DIST)
C        IF (DIST.GT.RADIUS) THEN
C          EVAP=.TRUE.
C          WRITE(*,*) "EVAP",DIST,RADIUS
C          IF(SW.EQ.11)RETURN
C          ENR=ENR+((DIST-RADIUS)/RADIUS)**2
C        ENDIF
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
              Dij=Dii
              nij=ni
              mij=mi
              aij=aii
            ELSE IF ( (J1.GT.A ) .AND .(J2.GT.A ) ) THEN
              Dij=Djj
              nij=nj
              mij=mj
              aij=ajj
            ELSE
              Dij=Dij0
              nij=nij0
              mij=mij0
              aij=aij0
            ENDIF
C            write(*,*) Dij,cij,mij,nij,aij,DIST
            ATT_DUMMY=(aij/DIST)**mij
            REP_DUMMY=(aij/DIST)**nij
            ATT(J1)=ATT(J1)+ATT_DUMMY
            REP(J1)=REP(J1)+Dij*REP_DUMMY/2.D0
            ENR_ATT=ENR_ATT+ATT_DUMMY
            ENR_REP=ENR_REP+Dij*REP_DUMMY/2.D0
            GRAD_REPM(J1,J2)=(-1.D0/2.D0)*nij*Dij*
     &(REP_DUMMY)/(DIST**2.D0)
            GRAD_ATTM(J1,J2)=(1.D0/2.D0)*mij*Cij*Dij2
     &*ATT_DUMMY/(DIST**2.D0)
          ENDIF
        ENDDO
        GRAD_ATT(J1)=DSQRT(ENR_ATT)
        ENR_REPI=ENR_REPI+ENR_REP
        ENR_ATTI=ENR_ATTI+GRAD_ATT(J1)*Cij*Dij2
        VT(J1)=REP(J1)-Cij*Dij2*DSQRT(ATT(J1))
        IF (VT(J1).GT.VMAX)THEN
          VMAX=VT(J1)
          JMAX=J1
        ENDIF
        IF (VT(J1).LT.VMIN) VMIN=VT(J1)
      ENDDO
      ENR=(ENR_REPI)-(ENR_ATTI)
C      write(*,*) "ENR=",ENR,ENR_REPI,ENR_ATTI
C      pause
      VTR(JMAX)=VTR(JMAX)+1.D0
      DO J1=1,NATOMS
        J3=J1*3
        DUMMYX=0.0D0
        DUMMYY=0.0D0
        DUMMYZ=0.0D0
        DO J2=1,NATOMS
          J4=J2*3
          XMUL2=GRAD_REPM(J1,J2)+(GRAD_ATTM(J1,J2)/GRAD_ATT(J1))
          XMUL3=GRAD_REPM(J2,J1)+(GRAD_ATTM(J2,J1)/GRAD_ATT(J2))
          DUMMYX=DUMMYX+(XMUL2+XMUL3)*(X(J3-2)-X(J4-2))
          DUMMYY=DUMMYY+(XMUL2+XMUL3)*(X(J3-1)-X(J4-1))
          DUMMYZ=DUMMYZ+(XMUL2+XMUL3)*(X(J3)  -X(J4))
        ENDDO   
C        DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
C        IF (DIST.GT.RADIUS) THEN
C          DUMMYX=DUMMYX+((DIST-RADIUS)/RADIUS)*X(J3-2)
C          DUMMYY=DUMMYY+((DIST-RADIUS)/RADIUS)*X(J3-1) 
C          DUMMYZ=DUMMYZ+((DIST-RADIUS)/RADIUS)*X(J3)
C        ENDIF  
        GRAD(J3-2)=DUMMYX
        GRAD(J3-1)=DUMMYY
        GRAD(J3)=DUMMYZ
      ENDDO
      RETURN
      END
