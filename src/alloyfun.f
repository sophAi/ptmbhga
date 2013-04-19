       FUNCTION ALLOYFUNC(NDIM,X,EVAP)
       implicit none
       LOGICAL EVAP
       INTEGER NDIM,SW
       REAL*8 ALLOYFUNC,X(NDIM),EREAL
       COMMON/cgsw/ SW
*************************************
*  Put your energy subroutine here  *
*  This function is for derivative  *
*************************************
       CALL ALLOY(NDIM,X,EREAL,EVAP)
*************************************
       IF (SW.EQ.11.AND.EVAP) RETURN
       ALLOYFUNC=EREAL
       END
