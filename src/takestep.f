      SUBROUTINE TAKESTEP
      IMPLICIT NONE
      
      INCLUDE 'params.h'
      INCLUDE 'commons.h'

      REAL*8 DPRAND, RANDOM, XMASS, YMASS, ZMASS, 
     1                 DIST(3*MXATMS), DMAX, VMAX, VMIN,
     2                 THETA, PHI, PI, DUMMY
      PARAMETER (PI=3.141592654D0)
      INTEGER J1, J2, JMAX
      LOGICAL CENT
      COMMON /cent/CENT
C     INTEGER NQTOT
C     COMMON /TOT/ NQTOT

C
C  This call can be used to keep the random numbers the same for parallel
C  runs - only for testing!
C
C     CALL SDPRND(NQTOT)

      DO J1=1,3*(NATOMS-NSEED)
         COORDSO(J1)=COORDS(J1)
      ENDDO
      DO J1=1,NATOMS
         VATO(J1)=VAT(J1)
      ENDDO

      DMAX=-1.0D0
      VMAX=-1.0D6
      VMIN=1.0D6
      DO J1=1,NATOMS
         J2=3*J1
         DIST(J1)=DSQRT(COORDS(J2-2)**2+COORDS(J2-1)**2+COORDS(J2)**2)
         IF (DIST(J1).GT.DMAX) DMAX=DIST(J1)
         IF (VAT(J1).GT.VMAX) THEN
            VMAX=VAT(J1)
            JMAX=J1
         ENDIF
         IF (VAT(J1).LT.VMIN) VMIN=VAT(J1)
      ENDDO
C     WRITE(*,*) "ASTEP=",ASTEP
C     WRITE(*,*) "STEP=",STEP
C     PAUSE"XANADU"
      DO J1=1,NATOMS-NSEED
10       J2=3*J1
         IF ((((VAT(J1).GT.ASTEP*VMIN).AND.(J1.EQ.JMAX)).OR.
     &(NATOMS-NSEED.EQ.1))) THEN
           THETA=DPRAND()*PI
           PHI=DPRAND()*PI*2.0D0
           COORDS(J2-2)=DMAX*DSIN(THETA)*DCOS(PHI)
           COORDS(J2-1)=DMAX*DSIN(THETA)*DSIN(PHI)
           COORDS(J2)=  DMAX*DCOS(THETA)
         ELSE
           RANDOM=(DPRAND()-0.5D0)*2.0D0
           COORDS(J2-2)=COORDS(J2-2)+STEP*RANDOM
           RANDOM=(DPRAND()-0.5D0)*2.0D0
           COORDS(J2-1)=COORDS(J2-1)+STEP*RANDOM
           RANDOM=(DPRAND()-0.5D0)*2.0D0
           COORDS(J2)=COORDS(J2)+STEP*RANDOM
C
C Stop atoms leaving the container in this step
C
           DUMMY=COORDS(J2-2)**2+COORDS(J2-1)**2+COORDS(J2)**2
           IF (DUMMY.GT.RADIUS) THEN
              COORDS(J2-2)=COORDS(J2-2)*DSQRT(RADIUS/DUMMY)
              COORDS(J2-1)=COORDS(J2-1)*DSQRT(RADIUS/DUMMY)
              COORDS(J2)=COORDS(J2)*DSQRT(RADIUS/DUMMY)
           ENDIF
         ENDIF
      ENDDO
C
C  Preserve centre of mass if required.
C
      IF (CENT.AND.(.NOT.SEEDT)) THEN
         XMASS=0.0D0
         YMASS=0.0D0
         ZMASS=0.0D0
         DO J1=1,NATOMS
            J2=3*J1
            XMASS=XMASS+COORDS(J2-2)
            YMASS=YMASS+COORDS(J2-1)
            ZMASS=ZMASS+COORDS(J2)
         ENDDO
         XMASS=XMASS/(NATOMS)
         YMASS=YMASS/(NATOMS)
         ZMASS=ZMASS/(NATOMS)
         DO J1=1,NATOMS
            J2=3*J1
            COORDS(J2-2)=COORDS(J2-2)-XMASS
            COORDS(J2-1)=COORDS(J2-1)-YMASS
            COORDS(J2)=  COORDS(J2)-ZMASS
         ENDDO
      ENDIF

      RETURN
      END
