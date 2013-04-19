C
C  Energy and gradient for Many-Body Potential.
C
      SUBROUTINE ADFUNCGP(NDIM,X,GRAD)
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4,NATOMS,NDIM
      REAL*8 ZETA,EPSILON,P,Q,RZERO
      REAL*8 X(NDIM),GRAD(NDIM)
      REAL*8 GRAD_ATT(NDIM/3),ELJ_ATT,XMUL3
      REAL*8 DIST
      REAL*8 GRAD_ATTM(6000,6000)
      REAL*8 GRAD_REPM(6000,6000)
      REAL*8 DUMMYX, DUMMYY, DUMMYZ, XMUL2
      COMMON /parameter/ EPSILON,ZETA,P,Q,RZERO
C     COMMON /accuracy3/ RADIUS,NONE
      NATOMS=NDIM/3
      DO J1=1,NATOMS
         J3=3*J1
         DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
C        IF (DIST.GT.RADIUS) THEN
C           EVAP=.TRUE.
C           ELJ=ELJ+(DIST-RADIUS)**2.
C           WRITE(*,*) "DIST=",DIST,RADIUS
C           PAUSE "EVP!!"
C        ENDIF
         GRAD_REPM(J1,J1)=0.0D0
         GRAD_ATTM(J1,J1)=0.0D0
         ELJ_ATT=0.0D0
         DO J2=1,NATOMS
            J4=3*J2
C   J3=I,J4=J
            IF(J1.NE.J2)THEN
              DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2
     &+(X(J3)-X(J4))**2
              DIST=DSQRT(DIST)
C             WRITE(*,*) DIST
            GRAD_REPM(J1,J2)=EPSILON*(-1.D0)*P*DEXP(P*(1.D0-DIST))/DIST
           GRAD_ATTM(J1,J2)=ZETA*(-2.D0)*Q*DEXP(2.D0*Q*(1.D0-DIST))/DIST
            ELJ_ATT=ELJ_ATT+GRAD_ATTM(J1,J2)*DIST/(-2.D0*Q*ZETA)
            ENDIF
         ENDDO
         GRAD_ATT(J1)=DSQRT(ELJ_ATT)
      ENDDO

      DO J1=1,NATOMS
         J3=3*J1
         DUMMYX=0.0D0
         DUMMYY=0.0D0
         DUMMYZ=0.0D0
         DO J2=1,NATOMS
           J4=3*J2
         XMUL2=GRAD_REPM(J1,J2)-(GRAD_ATTM(J1,J2)/(2.D0*GRAD_ATT(J1)))
         XMUL3=GRAD_REPM(J2,J1)-(GRAD_ATTM(J2,J1)/(2.D0*GRAD_ATT(J2)))
           DUMMYX=DUMMYX+(XMUL2+XMUL3)*(X(J3-2)-X(J4-2))
           DUMMYY=DUMMYY+(XMUL2+XMUL3)*(X(J3-1)-X(J4-1))
           DUMMYZ=DUMMYZ+(XMUL2+XMUL3)*(X(J3)  -X(J4))
         ENDDO
C     FINISH VIJ
         DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
C     PASS
C        IF (DIST.GT.RADIUS) THEN
C          DUMMYX=DUMMYX+(DIST-RADIUS)*X(J3-2)
C          DUMMYY=DUMMYY+(DIST-RADIUS)*X(J3-1)
C          DUMMYZ=DUMMYZ+(DIST-RADIUS)*X(J3)
C          WRITE(*,*) "WARNING!!RADIUS=",RADIUS,",DIST=",DIST
C        ENDIF
C     ENDPASS
C        FAC=DSQRT(DUMMYX**2+DUMMYY**2+DUMMYZ**2)
         GRAD(J3-2)=DUMMYX
         GRAD(J3-1)=DUMMYY
         GRAD(J3)=DUMMYZ
C TEST OUTPUT
C        WRITE(*,*) "OUTPUT TEST"
C        WRITE(*,*) V(J3-2),V(J3-1),V(J3)
C        PAUSE
      ENDDO
      RETURN
      END
