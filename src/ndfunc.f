C
C  Numerical derivative of an input Potential
C
      SUBROUTINE NDFUNC(NDIM,X,GRAD,ENR,FUNCT)
      IMPLICIT NONE
      INTEGER J1,J2,NDIM,SW,type 
      REAL*8 X_DX(NDIM),DX,ENR_DX
      REAL*8 X(NDIM), DIST, GRAD(NDIM) 
      REAL*8 ENR, DUMMYX, DUMMYY, DUMMYZ
      REAL*8 FUNCT,NONE1,BC,RADIUS_TEMP,RADIUS
      LOGICAL EVAP
      EXTERNAL FUNCT
      COMMON /accuracy3/ RADIUS,NONE1
      COMMON /accuracy4/ DX,BC
      COMMON /cgsw/ SW
      COMMON /engtype/ type
      EVAP=.FALSE.
      ENR_DX=0.0D0
      IF(SW.EQ.11)THEN    
        RADIUS_TEMP=RADIUS
      ELSE
        RADIUS_TEMP=RADIUS
        RADIUS=10000000.D0
      ENDIF
      DO J1=1,NDIM
        X_DX(J1)=X(J1)
      ENDDO
      DO J1=1,NDIM
        X_DX(J1)=X_DX(J1)+DX
****************************************
*    INSERT YOUR POTENTIAL HERE        *
****************************************
        ENR_DX=FUNCT(NDIM,X_DX,EVAP)
        GRAD(J1)=(ENR_DX-ENR)/(DX*BC)
        X_DX(J1)=X(J1)
C        WRITE(*,*) "ENR_FX,ENR,GRAD IN NDFUNC=",ENR_DX,ENR,GRAD(J1)
      ENDDO
      RADIUS=RADIUS_TEMP
      IF(SW.EQ.11.or.type.eq.7) GOTO 13
      DO J1=1,NDIM/3
        J2=3*J1
        DUMMYX=GRAD(J2-2)
        DUMMYY=GRAD(J2-1)
        DUMMYZ=GRAD(J2)
        DIST=X(J2-2)**2+X(J2-1)**2+X(J2)**2
        IF (DIST.GT.RADIUS) THEN
          DUMMYX=DUMMYX+(DIST-RADIUS)*X(J2-2)/BC
          DUMMYY=DUMMYY+(DIST-RADIUS)*X(J2-1)/BC
          DUMMYZ=DUMMYZ+(DIST-RADIUS)*X(J2)/BC
        ENDIF
        GRAD(J2-2)=DUMMYX
        GRAD(J2-1)=DUMMYY
        GRAD(J2)=DUMMYZ
C       WRITE(*,*) "GRAD IN DFUNC=",GRAD(J2-2),GRAD(J2-1),GRAD(J2)
      ENDDO
C     PAUSE
13    RETURN
      END
