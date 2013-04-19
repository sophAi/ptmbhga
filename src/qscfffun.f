       FUNCTION QSCFFFUNC(NDIM,X,EVAP)
       implicit none
       LOGICAL EVAP
       INTEGER NDIM,SW
       REAL*8 QSCFFFUNC,X(NDIM),EREAL
       COMMON/cgsw/ SW
*************************************
*  Put your energy subroutine here  *
*  This function is for derivative  *
*************************************
       CALL QSCFF(NDIM,X,EREAL,EVAP)
*************************************
       IF (SW.EQ.11.AND.EVAP) RETURN
       QSCFFFUNC=EREAL
       END
