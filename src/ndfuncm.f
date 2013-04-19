C
C  Numerical derivative of an input Potential
C
      SUBROUTINE NDFUNCM(NDIM,X,GRAD,ENR,FUNCT)
      IMPLICIT NONE
      INTEGER J1,J2,NDIM,SW 
      DOUBLE PRECISION X_DX(NDIM),DX,ENR_DX
      DOUBLE PRECISION X(NDIM), DIST, GRAD(NDIM) 
      DOUBLE PRECISION ENR, DUMMYX, DUMMYY, DUMMYZ
      DOUBLE PRECISION FUNCT,NONE1,BC,RADIUS_TEMP,RADIUS
      LOGICAL EVAP
      EXTERNAL FUNCT
      COMMON /accuracy3/ RADIUS,NONE1
      COMMON /accuracy4/ DX,BC
      COMMON /cgsw/ SW
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
        ENR_DX=FUNCT(NDIM,X_DX)
        GRAD(J1)=(ENR_DX-ENR)/DX
        X_DX(J1)=X(J1)
      ENDDO
      RADIUS=RADIUS_TEMP
      IF(SW.EQ.11) GOTO 13
      DO J1=1,NDIM/3
        J2=3*J1
        DUMMYX=GRAD(J2-2)
        DUMMYY=GRAD(J2-1)
        DUMMYZ=GRAD(J2)
        DIST=X(J2-2)**2+X(J2-1)**2+X(J2)**2
        IF (DIST.GT.RADIUS) THEN
          DUMMYX=DUMMYX+(DIST-RADIUS)*X(J2-2)
          DUMMYY=DUMMYY+(DIST-RADIUS)*X(J2-1)
          DUMMYZ=DUMMYZ+(DIST-RADIUS)*X(J2)
        ENDIF
        GRAD(J2-2)=DUMMYX/BC
        GRAD(J2-1)=DUMMYY/BC
        GRAD(J2)=DUMMYZ/BC
C       WRITE(*,*) "GRAD IN DFUNC=",GRAD(J2-2),GRAD(J2-1),GRAD(J2)
      ENDDO
C     PAUSE
13    RETURN
      END
