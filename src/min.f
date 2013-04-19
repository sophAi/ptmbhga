        SUBROUTINE LOCAL_MIN(XINIB,MIN_POT,NDIM,SEED,CFLAG,MIN_S,
     &MIN_T,EVAP)
        IMPLICIT NONE
        INTEGER pdb_rec,dim3,J1,J2
        INTEGER MIN,NDIM,SEED,MIN_S,MIN_T
        REAL*8 XINIB(NDIM),MIN_POT,RADIUS,NONE,DG,DIST
        COMMON/accuracy3/ RADIUS,NONE
        COMMON/guide1/ DG
        COMMON/pdbrec/ pdb_rec,dim3
        COMMON/cent/ CENT2
        COMMON/min1/ min
        LOGICAL CFLAG,EVAP,CENT2
C        DO J1=1,NDIM/3
C          J2=J1*3
C          WRITE(*,*) "BEFORE=",XINIB(J2-2),XINIB(J2-1),XINIB(J2)
C        ENDDO
        IF (MIN.EQ.1)THEN
          call simplex_min(XINIB,MIN_POT,NDIM,SEED,CFLAG,EVAP)
        ELSE IF(min.eq.2)THEN
          call conj_min(XINIB,MIN_POT,NDIM,SEED,CFLAG,EVAP)
        ELSE IF(min.eq.3)THEN
          call conj_min(XINIB,MIN_POT,NDIM,SEED,CFLAG,EVAP)
          call simplex_min(XINIB,MIN_POT,NDIM,SEED,CFLAG,EVAP)
        ELSE IF(min.eq.4)THEN
          call simplex_min(XINIB,MIN_POT,NDIM,SEED,CFLAG,EVAP)
          call conj_min(XINIB,MIN_POT,NDIM,SEED,CFLAG,EVAP)
        ELSE
          CALL BMIN(XINIB,MIN_POT,NDIM,SEED,CFLAG,EVAP)
C          IF(MIN_POT.gt.-99999999.D0.and.MIN_POT.lt.9999999999.D0)then
C          ELSE
C          pause "NAN!!!"
C          ENDIF
        ENDIF
C        WRITE(*,*) CFLAG
        IF(CENT2.and.dim3.ne.1)THEN
          CALL CENTRE(XINIB,NDIM/3)
        ENDIF
        IF(CFLAG) THEN
          MIN_S=MIN_S+1
          MIN_T=MIN_T+1
          IF (dim3.eq.1) THEN
C            DO J1=1,NDIM
C              GBIN1=INT((XINIB(J1)+DSQRT(RADIUS))/DG)+1
C              GUIDE(GBIN)=GUIDE(GBIN)+1
C              GBIN_TOTAL=GBIN_TOTAL+1
C            ENDDO
          ELSE
            DO J1=1,NDIM/3
              J2=J1*3
              DIST=XINIB(J2-2)**2+XINIB(J2-1)**2+XINIB(J2)**2
              IF (DIST.LE.RADIUS)THEN
C                GBIN1=INT((XINIB(J2-2)+DSQRT(RADIUS))/DG)
C                GBIN2=INT((XINIB(J2-1)+DSQRT(RADIUS))/DG)
C                GBIN3=INT((XINIB(J2)+DSQRT(RADIUS))/DG)
C                GBIN_N=(1+GBIN1+GBIN2*GBIN+GBIN3*GBIN*GBIN)
CC                GUIDE(GBIN_N)=GUIDE(GBIN_N)+1
C                GBIN_TOTAL=GBIN_TOTAL+1
C          WRITE(*,*) ,J1,GBIN_N,GBIN1,GBIN2,GBIN3,XINIB(J2-2),
C     &XINIB(J2-1),XINIB(J2)
              ENDIF
            ENDDO
          ENDIF
        ELSE
          MIN_T=MIN_T+1
C          WRITE(*,*) "Local minimum searching did not converge!"
        ENDIF
C        PAUSE
C         write(*,*) "MIN_POT IN min.f=",MIN_POT
        RETURN
        END
