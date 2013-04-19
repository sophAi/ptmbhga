        SUBROUTINE LJ(NATOMS,X,ENR,EVAP)
        IMPLICIT NONE
        INTEGER NATOMS,J1,J2,J3,J4,SW,JMAX
        INTEGER I_TYPE(10000),TYPE
        REAL*8 X(NATOMS*3),ENR,VT(2000),DIST,DIST1,DIST2,DIST3,CUTOFF
        REAL*8 FACTOR1,FACTOR2,RADIUS,NONE,ECOIL,VTR(2000)
        LOGICAL EVAP
        COMMON /LONG_RANGE1/ CUTOFF,ECOIL
        COMMON /accuracy3/ RADIUS,NONE
        COMMON /cgsw/ SW
        COMMON /engtype/ TYPE
        COMMON /INTER_TYPE1/ I_TYPE
        COMMON /vt1/ VT,VTR
        COMMON /vt2/ JMAX
C       I_TYPE=1  B
C       I_TYPE=2  L
C       I_TYPE=3  N
C       RADIUS=500.
C       DO J1=1,NATOMS
C         WRITE(*,*) J1,I_TYPE(J1)
C       ENDDO
        ENR=0.D0
        EVAP=.FALSE.
        DO J1=1,NATOMS
          J3=J1*3
          VT(J1)=0.0D0
          DO J2=J1+1,NATOMS
            J4=J2*3
            DIST1=(X(J3)-X(J4))**2.
            DIST2=(X(J3-1)-X(J4-1))**2.
            DIST3=(X(J3-2)-X(J4-2))**2.
        IF (DIST1.LT.CUTOFF.AND.DIST2.LT.CUTOFF.AND.DIST3.LT.CUTOFF)THEN
               DIST=DIST1+DIST2+DIST3
              IF (DSQRT(DIST).GT.RADIUS)THEN
                EVAP=.TRUE.
                IF (SW.EQ.11)RETURN
                ENR=ENR+(DIST-RADIUS)**2.
              ENDIF
              IF (I_TYPE(J1).EQ.1.AND.I_TYPE(J2).EQ.1)THEN
C  BB           
                FACTOR1=1.D0
                FACTOR2=1.D0
              ELSE IF (I_TYPE(J1).EQ.2.OR.I_TYPE(J2).EQ.2)THEN
C  LB,BL,LL     
                IF (TYPE.EQ.3)THEN      
                  FACTOR1=2.D0/3.D0
                  FACTOR2=(-2.D0)/3.D0
                ENDIF
                IF (TYPE.EQ.4)THEN
                  FACTOR1=1.D0
                  FACTOR2=0.D0
                ENDIF
              ELSE
C  NB,BN,NL,LN,NN
                IF (TYPE.EQ.3)THEN
                  FACTOR1=1.D0
                  FACTOR2=0.D0
                ENDIF
                IF (TYPE.EQ.4)THEN
                  FACTOR1=1.D0
                  FACTOR2=0.D0
                ENDIF
              ENDIF
              ENR=ENR+4.D0*ECOIL*((FACTOR1*DIST**(-6.))-
     &(FACTOR2*DIST**(-3.)))
              VT(J1)=VT(J1)+4.D0*ECOIL*((FACTOR1*DIST**(-6.))-
     &(FACTOR2*DIST**(-3.)))
C             write(*,*) "ENR_LJ=",ENR,ECOIL,DIST
            ENDIF
          ENDDO
        ENDDO
        RETURN
        END
                
              
