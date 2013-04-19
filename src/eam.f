************************************************************
*  Energy and Gradient for EAM and alloy Potential.        *
************************************************************
      SUBROUTINE EAM(NDIM,X,ELJ,EVAP)
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4,J5,J6,NATOMS,SW,NDIM,JMAX,AN,P1,P3,Q1,Q3
      REAL*8 A,EPSILON,C,N,M,SUMP,SUMQ
      REAL*8 A1,EPSILON1,C1,N1,M1,GRAD(6000)
      REAL*8 A2,EPSILON2,C2,N2,M2,A3,EPSILON3,C3,N3,M3
      REAL*8 PHI(NDIM),F(NDIM),PHI_DUMMY,F_DUMMY
      REAL*8 ELJ_PHI,ELJ_F,VT(2000),RADIUS,VTR(2000)
      REAL*8 ELJ_PHII,ELJ_FI,NONE
      REAL*8 X(NDIM),DIST,ELJ,GRAD_PHI(2000,2000)
      REAL*8 GRAD_F(2000,2000),DISTP,DISTQ
      REAL*8 P_DUMMY,Q_DUMMY 
      REAL*8 DUMMYX,DUMMYY,DUMMYZ,XMUL2,XMUL3,F_DUMMYI
      REAL*8 VMIN,VMAX,DMAX
      LOGICAL EVAP,GRADT
      COMMON /alloy3/ AN
      COMMON /accuracy3/ RADIUS,NONE
      COMMON /cgsw/ SW
      COMMON /GRAD1/ GRAD
      COMMON /GRAD2/ GRADT
      COMMON /vt1/ VT,VTR
      COMMON /vt2/ JMAX
      COMMON /vt3/ VMIN,VMAX,DMAX
      COMMON /eam_parameter1/ A,EPSILON,C,N,M
c      write(*,*) a1,a2,a3,c1,c2,c3,n1,n2,n3,m1,m2,m3,epsilon1,epsilon2
c     &           ,epsilon3
      NATOMS=NDIM/3
      EVAP=.FALSE.
      VMAX=-1.0D6
      VMIN=1.0D6
      DMAX=-1.0D0
      ELJ=0.0D0
      ELJ_PHII=0.0D0
      ELJ_FI=0.0D0
      DO J1=1,NATOMS
        PHI(J1)=0.0D0
        F(J1)=0.0D0
        VT(J1)=0.0D0
        VTR(J1)=0.0D0
      END DO
      DO J1=1,NATOMS
        J3=3*J1
        DIST=X(J3-2)**2.D0+X(J3-1)**2.D0+X(J3)**2.D0
        IF (DSQRT(DIST).GT.DMAX) DMAX=DSQRT(DIST)
        IF (DIST.GT.RADIUS) THEN
          EVAP=.TRUE.
          IF (SW.EQ.11)RETURN
C            ELJ=ELJ+(DIST-RADIUS)**2
        ENDIF
        ELJ_PHI=0.0D0
        ELJ_F=0.0D0
        PHI(J1)=0.0D0
        F(J1)=0.0D0
        GRAD_PHI(J1,J1)=0.0D0
        GRAD_F(J1,J1)=0.0D0
        DO J2=1,NATOMS
          J4=3*J2
          IF(J2.NE.J1)THEN
            DIST=(X(J3-2)-X(J4-2))**2.D0+(X(J3-1)-X(J4-1))**2.D0+
     &(X(J3)-X(J4))**2.D0
            DIST=DSQRT(DIST)
c            write(*,*) "DIST=",DIST
            PHI_DUMMY = ((A/DIST))**N
            F_DUMMY = ((A/DIST))**M  
            PHI(J1)= PHI(J1)+ PHI_DUMMY
            F(J1)= F(J1)+ F_DUMMY
            ELJ_PHI= ELJ_PHI+PHI_DUMMY
            ELJ_F= ELJ_F+F_DUMMY
          ENDIF
        ENDDO
        F_DUMMYI=DSQRT(ELJ_F)
        ELJ_PHII=ELJ_PHII+ELJ_PHI
        ELJ_FI=ELJ_FI+F_DUMMYI        
        VT(J1)=((1.D0/2.D0)*EPSILON*(PHI(J1)))-(EPSILON*C*DSQRT(F(J1)))
        IF (VT(J1).GT.VMAX)THEN
          VMAX=VT(J1)
          JMAX=J1
        ENDIF
        IF (VT(J1).LT.VMIN) VMIN=VT(J1)
      ENDDO
      ELJ=((1.D0/2.D0)*EPSILON*ELJ_PHII)-(EPSILON*C*ELJ_FI)
c      write(*,*) "ELJ=",ELJ,ELJ_PHII*EPSILON*0.5D0,ELJ_FI*epsilon*c
      VTR(JMAX)=VTR(JMAX)+1.D0
      DO J1=1,NATOMS
        J3=3*J1
        DIST=X(J3-2)**2.D0+X(J3-1)**2.D0+X(J3)**2.D0
        IF (DSQRT(DIST).GT.DMAX) DMAX=DSQRT(DIST)
        IF (DIST.GT.RADIUS) THEN
          EVAP=.TRUE.
          IF (SW.EQ.11)RETURN
C            ELJ=ELJ+(DIST-RADIUS)**2
        ENDIF
        DUMMYX=0.0D0
        DUMMYY=0.0D0
        DUMMYZ=0.0D0
        DO J2=1,NATOMS
          J4=3*J2
          IF (J2.NE.J1)THEN
            DIST=(X(J3-2)-X(J4-2))**2.D0+(X(J3-1)-X(J4-1))**2.D0+
     &(X(J3)-X(J4))**2.D0
            DIST=DSQRT(DIST)
C            write(*,*) "DIST=",DIST
            GRAD_PHI(J1,J2)=(-1.D0*EPSILON)*N*(A**N)*
     &((DIST)**(-1.D0*N-2.D0))
            SUMP=0.0D0
            DO P1=1,NATOMS
              P3=3*P1
              IF(P1.NE.J1)THEN
                DISTP=(X(J3-2)-X(P3-2))**2.D0+(X(J3-1)-X(P3-1))**2.D0+
     &(X(J3)-X(P3))**2.D0
                DISTP=DSQRT(DISTP)
C                write(*,*) "DISTP=",DISTP
                P_DUMMY=(1.D0/DISTP)**M
                SUMP=SUMP+P_DUMMY
              ENDIF
            ENDDO
            SUMQ=0.0D0
            DO Q1=1,NATOMS
              Q3=3*Q1
              IF (Q1.NE.J2)THEN
                DISTQ=(X(J4-2)-X(Q3-2))**2.D0+(X(J4-1)-X(Q3-1))**2.D0+
     &(X(J4)-X(Q3))**2.D0
                DISTQ=DSQRT(DISTQ)
C                write(*,*) "DISTQ=",DISTQ
                Q_DUMMY=(1.D0/DISTQ)**M
                SUMQ=SUMQ+Q_DUMMY
              ENDIF
            ENDDO
C            write(*,*) "SUMQ=",SUMQ,",SUMP=",SUMP
            GRAD_F(J1,J2)= (1.D0/2.D0)*EPSILON*C*(A**(M/2.D0))*M*
     &(DIST)**(-1.D0*M-2.D0)*((SUMP)**(-1.D0/2.D0)+(SUMQ)**
     &(-1.D0/2.D0))
        DUMMYX=DUMMYX+(GRAD_PHI(J1,J2)+GRAD_F(J1,J2))*(X(J3-2)-X(J4-2))
        DUMMYY=DUMMYY+(GRAD_PHI(J1,J2)+GRAD_F(J1,J2))*(X(J3-1)-X(J4-1))
        DUMMYZ=DUMMYZ+(GRAD_PHI(J1,J2)+GRAD_F(J1,J2))*(X(J3)-X(J4))
          ENDIF
        ENDDO
C        DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
C        IF (DIST.GT.RADIUS) THEN
C          DUMMYX=DUMMYX+(DIST-RADIUS)*X(J3-2)
C          DUMMYY=DUMMYY+(DIST-RADIUS)*X(J3-1)
C          DUMMYZ=DUMMYZ+(DIST-RADIUS)*X(J3)
C        ENDIF
        GRAD(J3-2)=DUMMYX
        GRAD(J3-1)=DUMMYY
        GRAD(J3)=DUMMYZ
C        write(*,*) "grad",grad(j3-2),grad(j3-1),grad(j3)
      ENDDO
C      pause
      RETURN
      END

