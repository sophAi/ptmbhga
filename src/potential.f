      SUBROUTINE POTENTIAL(NDIM,X,GRAD,EREAL,RMS,GRADT,EVAP)
      IMPLICIT NONE 
      LOGICAL GRADT, EVAP ,GRADT2
      INTEGER J1, J2, NDIM,TYPE,SW,ASW
      REAL*8 EREAL, GRAD(NDIM), X(NDIM)
      REAL*8 RMS,XG,YG,ZG,DUMMY2,GRAD_TEMP(6000)
      REAL*8 GPFUNC,ALLOYFUNC,PROTFUNC,EAMFUNC,PROTOH_FUNC,
     &CHIMA_FUNC,QSCFFFUNC
      EXTERNAL GPFUNC,ALLOYFUNC,PROTFUNC,EAMFUNC,PROTOH_FUNC,
     &CHIMA_FUNC,QSCFFFUNC
      COMMON /engtype/ TYPE
      COMMON /cgsw/ SW
      COMMON /anatsw/ ASW 
      COMMON /GRAD1/ GRAD_TEMP
      COMMON /GRAD2/ GRADT2
****************************************
*          Select energy type          *
****************************************
      GRADT2=GRADT
      If (TYPE.EQ.1)THEN
        EREAL=GPFUNC(NDIM,X,EVAP)
        DO J1=1,NDIM
          GRAD(J1)=GRAD_TEMP(J1)
        ENDDO
        IF(SW.EQ.11.AND.EVAP) RETURN
        IF (GRADT) THEN
          IF (ASW.NE.1)THEN
            CALL NDFUNC(NDIM,X,GRAD,EREAL,GPFUNC)
          ENDIF
        ENDIF
      ENDIF
      IF (TYPE.EQ.2)THEN
        EREAL=ALLOYFUNC(NDIM,X,EVAP)
        DO J1=1,NDIM
          GRAD(J1)=GRAD_TEMP(J1)
        ENDDO
        IF (SW.EQ.11.AND.EVAP) RETURN
        IF (GRADT) THEN
          IF (ASW.NE.1)THEN
            CALL NDFUNC(NDIM,X,GRAD,EREAL,ALLOYFUNC)
          ENDIF 
        ENDIF
      ENDIF 
      IF (TYPE.EQ.3.OR.TYPE.eq.4.OR.TYPE.EQ.6)THEN
        EREAL=PROTFUNC(NDIM,X,EVAP)
        IF (SW.EQ.11.AND.EVAP) RETURN
        IF (GRADT) THEN
          IF (ASW.EQ.1)THEN
C           CALL ADFUNCAL...
          ELSE
            CALL NDFUNC(NDIM,X,GRAD,EREAL,PROTFUNC)
          ENDIF
        ENDIF
      ENDIF
      IF (TYPE.EQ.5)THEN
        EREAL=EAMFUNC(NDIM,X,EVAP)
        IF(ASW.EQ.1)THEN
          DO J1=1,NDIM
            GRAD(J1)=GRAD_TEMP(J1)
          ENDDO
        ENDIF
        IF(SW.EQ.11.AND.EVAP) RETURN
        IF (GRADT) THEN
         IF (ASW.NE.1)THEN
            CALL NDFUNC(NDIM,X,GRAD,EREAL,EAMFUNC)
          ENDIF
       ENDIF
      ENDIF
      IF (TYPE.EQ.7)THEN

        EREAL=PROTOH_FUNC(NDIM,X,EVAP)
C         WRITE(*,*) "POT_EVAP=",EVAP
C         WRITE(*,*) "EREAL=",EREAL
C         IF (EVAP.AND.SW.EQ.11)RETURN
          DO J1=1,NDIM
C             WRITE(*,*) "POT,",J1," X=",X(J1)
          ENDDO
C        ENDIF
        EVAP=.false.
        IF (GRADT) THEN
          IF (ASW.EQ.1)THEN
            
          ELSE
            CALL NDFUNC(NDIM,X,GRAD,EREAL,PROTOH_FUNC)
C            DO J1=1,NDIM
C              WRITE(*,*) "GRAD IN POTENTIAL=",GRAD(J1),X(J1),EREAL
C            ENDDO
C            PAUSE
          ENDIF
        ENDIF
      ENDIF
      IF (TYPE.EQ.8)THEN
        EREAL=CHIMA_FUNC(NDIM,X,EVAP)
        IF (GRADT) THEN
          IF (ASW.EQ.1)THEN
          ELSE
            CALL NDFUNC(NDIM,X,GRAD,EREAL,CHIMA_FUNC)
          ENDIF
        ENDIF
      ENDIF
      IF (TYPE.EQ.9)THEN
        EREAL=QSCFFFUNC(NDIM,X,EVAP)
        DO J1=1,NDIM
          GRAD(J1)=GRAD_TEMP(J1)
        ENDDO
        IF (SW.EQ.11.AND.EVAP) RETURN
        IF (GRADT) THEN
          IF (ASW.NE.1)THEN
            CALL NDFUNC(NDIM,X,GRAD,EREAL,QSCFFFUNC)
          ENDIF
        ENDIF
      ENDIF
********End of select energy type*********
      
C
C  Preserve centre of mass if required.
C
      IF(GRADT)THEN
        DUMMY2=0.0D0
        DO J1=1,NDIM
          DUMMY2=DUMMY2+GRAD(J1)**2
        ENDDO
        RMS=DSQRT(DUMMY2/(NDIM))
      ENDIF
      RETURN
      END
