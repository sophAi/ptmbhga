      FUNCTION PROTFUNC(NDIM,X,EVAP)
      LOGICAL EVAP
      INTEGER SW,NATOMS,NDIM
      REAL*8 PROTFUNC,X(NDIM)
      REAL*8 EREAL,EREAL1,EREAL2,EREAL3
      COMMON /cgsw/ SW
*************************************
*   Put your energy subroutine here *
*   This function is for derivative *
*************************************
      NATOMS=NDIM/3
      CALL LJ(NATOMS,X,EREAL1,EVAP)
      IF(SW.EQ.11.AND.EVAP) RETURN
      CALL DIH(NATOMS,X,EREAL2)
      CALL BOND(NATOMS,X,EREAL3)
      EREAL=EREAL1+EREAL2+EREAL3
      PROTFUNC=EREAL
C     WRITE(*,*) "EREAL=",EVAP,PROTFUNC,EREAL1,EREAL2,EREAL3,EREAL
C     PAUSE
      END
