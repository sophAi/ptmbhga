       FUNCTION MOMENTFUNC(NDIM,X)
       implicit double precision (A-H,O-Z)
       INTEGER NDIM
       REAL*8 MOMENTFUNC,X(NDIM),BeGamma
c       COMMON /cgsw/ SW
*************************************
*  Put your energy subroutine here  *
*  This function is for derivative  *
*************************************
c       CALL MOMENT(NDIM,X,EREAL)
       MOMENTFUNC=BeGamma(X,NDIM)
C      WRITE(*,*) "EREAL=",EREAL
C        WRITE(*,*) "EVAP EREAL IN GPFUNC=",EREAL
c       MOMENTFUNC=EREAL
       END
