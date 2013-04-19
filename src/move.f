        SUBROUTINE MOVE(RATE)
        IMPLICIT NONE
        REAL*8 DGUESS,GTOL,RADIUS,NONE
        REAL*8 STEP,ASTEP,TEMP,RATE,ACCRAT,SCALEFAC,TEMP_TEMP
        REAL*8 MOVE_FAC1,MOVE_FAC2,TEMP_MIN,TEMP_MAX,STEP_TEMP
        REAL*8 STEP_MIN,STEP_MAX,ASTEP_MIN,ASTEP_MAX,ASTEP_TEMP
        LOGICAL FIXBOTH,FIXTEMP,FIXSTEP
        COMMON /bmin1/ DGUESS,ASTEP,STEP,GTOL
        COMMON /STEP1/ FIXTEMP,FIXSTEP
        COMMON /TEMP1/ SCALEFAC,TEMP
        COMMON /accept2/ ACCRAT
        COMMON /accept3/ MOVE_FAC1,MOVE_FAC2
        COMMON /accept4/ STEP_MIN,STEP_MAX,ASTEP_MIN,ASTEP_MAX
        COMMON /accept5/ TEMP_MIN,TEMP_MAX
        COMMON /accept6/ STEP_TEMP,ASTEP_TEMP,TEMP_TEMP
        COMMON /accuracy3/ RADIUS,NONE
**********MUST  BETA0<TEMP<BETAF****************
        IF(RATE.GT.ACCRAT)THEN
          IF (TEMP.GT.TEMP_MIN)THEN
            IF (.NOT.FIXTEMP) TEMP=TEMP/MOVE_FAC1
          ELSE
            WRITE(*,*) "#### Temperature is too SAMLL!temp=",TEMP
            ACCRAT=ACCRAT+MOVE_FAC2
            TEMP=TEMP_TEMP
            WRITE(*,*) "#### Acceptance ratio is INCREASING,now=",ACCRAT
          ENDIF
          IF (STEP.LE.STEP_MAX.AND.ASTEP.LE.ASTEP_MAX)THEN
            IF (.NOT.FIXSTEP) THEN
              STEP=STEP*MOVE_FAC1
              ASTEP=ASTEP*MOVE_FAC1
C            WRITE(*,*) "INCREASE"
            ENDIF
          ELSE
            WRITE(*,*) "#### STEP OR ASTEP is too LARGE!radius=",
     &DSQRT(RADIUS)/2.D0
            WRITE(*,*) "#### STEP=",STEP,",ASTEP=",ASTEP
            ACCRAT=ACCRAT+MOVE_FAC2
            STEP=STEP_TEMP
            ASTEP=ASTEP_TEMP
            WRITE(*,*) "#### Acceptance ratio is INCREASING,now=",ACCRAT
          ENDIF
        ELSE IF(RATE.EQ.ACCRAT)THEN
        ELSE
          IF (TEMP.LE.TEMP_MAX)THEN
            IF (.NOT.FIXTEMP) TEMP=TEMP*MOVE_FAC1
          ELSE
            WRITE(*,*) "#### Temperature is too LARGE!temp=",TEMP
            ACCRAT=ACCRAT-MOVE_FAC2
            TEMP=TEMP_TEMP
            WRITE(*,*) "#### Acceptance ratio is DECREASING now=",ACCRAT
          ENDIF
          IF (STEP.GT.STEP_MIN.AND.ASTEP.GT.ASTEP_MIN)THEN
            IF (.NOT.FIXSTEP) THEN
              STEP=STEP/MOVE_FAC1
              ASTEP=ASTEP/MOVE_FAC1
C              WRITE(*,*) "DECREASE!"
            ENDIF
          ELSE
            WRITE(*,*) "#### STEP OR ASTEP is too SMALL!radius=",
     &DSQRT(RADIUS)/2.D0
            WRITE(*,*) "#### STEP=",STEP,",ASTEP=",ASTEP
            ACCRAT=ACCRAT-MOVE_FAC2
            STEP=STEP_TEMP
            ASTEP=ASTEP_TEMP
            WRITE(*,*) "#### Acceptance ratio is DECREASING,now=",ACCRAT
          ENDIF
        ENDIF
C       IF(.NOT.FIXTEMP) WRITE(*,'(A,F12.6)') ' Temperature is now:',TEMP
C        IF (.NOT.FIXSTEP) THEN
C          WRITE(*,'(A,F12.6,A,F12.6)')
C     &' Maximum centre-of-mass and angular steps arenow:',STEP,' and '
C     &,ASTEP,",RATE=",RATE
C       ENDIF
       TEMP=TEMP*SCALEFAC
       RETURN
       END
