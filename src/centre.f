C
C**************************************************************************
C
C  Subroutine CENTRE moves the centre of mass to the origin.
C
C*********************************************************************
C
C
C
C
      SUBROUTINE CENTRE(X,NATOMS)
      IMPLICIT NONE
      INTEGER I,NATOMS
      REAL*8 XMASS, YMASS, ZMASS, X(NATOMS*3)
      REAL*8 RADIUS,DX
      COMMON /accuracy4/ DX,RADIUS
      XMASS=0.0D0
      YMASS=0.0D0
      ZMASS=0.0D0
      DO 10 I=1,NATOMS
         XMASS=XMASS+X(3*(I-1)+1)
         YMASS=YMASS+X(3*(I-1)+2)
         ZMASS=ZMASS+X(3*(I-1)+3)
10    CONTINUE
      XMASS=XMASS/NATOMS
      YMASS=YMASS/NATOMS
      ZMASS=ZMASS/NATOMS
      DO 20 I=1,NATOMS
         X(3*(I-1)+1)=X(3*(I-1)+1)-XMASS
         X(3*(I-1)+2)=X(3*(I-1)+2)-YMASS
         X(3*(I-1)+3)=X(3*(I-1)+3)-ZMASS
20    CONTINUE
C
C  Check that all the atoms are in the container. If not then rescale.
C
C      DISTMAX=0.0D0
C      DO J1=1,NATOMS
C         DIST=X(3*J1-2)**2+X(3*J1-1)**2+X(3*J1)**2
C         IF (DIST.GT.DISTMAX) DISTMAX=DIST
C      ENDDO
C      IF (DISTMAX.GT.RADIUS) THEN
C         DISTMAX=DSQRT(DISTMAX/RADIUS)*0.99D0
C         DO J1=1,NATOMS
C            X(3*J1-2)=X(3*J1-2)*DISTMAX
C            X(3*J1-1)=X(3*J1-1)*DISTMAX
C            X(3*J1)  =X(3*J1)  *DISTMAX
C         ENDDO
C      ENDIF
      RETURN
      END
