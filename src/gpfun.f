       FUNCTION GPFUNC(NDIM,X,EVAP)
       implicit none
       LOGICAL EVAP
       INTEGER NDIM,SW
       DOUBLE PRECISION GPFUNC,X(NDIM),EREAL
       COMMON /cgsw/ SW
*************************************
*  Put your energy subroutine here  *
*  This function is for derivative  *
*************************************
       CALL GP(NDIM,X,EREAL,EVAP)
C      WRITE(*,*) "EREAL=",EREAL
       IF(SW.EQ.11.AND.EVAP) RETURN
C        WRITE(*,*) "EVAP EREAL IN GPFUNC=",EREAL
       GPFUNC=EREAL
       RETURN
       END
