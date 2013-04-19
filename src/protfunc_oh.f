      FUNCTION PROTOH_FUNC(NDIM,X,EVAP)
      INTEGER NATOMS,NDIM,SW
      REAL*8 PROTOH_FUNC,X(NDIM),radius,none1
      REAL*8 EREAL1,EREAL2,EREAL0
      COMMON /cgsw/SW
      COMMON /accuracy3/radius,none1
      LOGICAL EVAP
*************************************
*   Put your energy subroutine here *
*   This function is for derivative *
*************************************
C     CALL LJ(NATOMS,X,EREAL0,VT,EVAP)
C     CALL BOND(NATOMS,X,EREAL2)
      CALL oh_bond(ndim,X,EREAL1,EVAP)
      PROTOH_FUNC=EREAL1
C     WRITE(*,*) "EREAL=",EVAP,PROTFUNC,EREAL1,EREAL2,EREAL3,EREAL
C     WRITE(*,*) "FUNC_EREAL=",EREAL1,",EVAP=",EVAP
      END
