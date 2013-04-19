*
*       Conjugate Gradient Driver
*
        SUBROUTINE BMIN(XINIB,POT,NDIM,SEED,CFLAG,EVAP)
        IMPLICIT NONE
        INTEGER MAXI,IFLAG,ITER,NDIM,NMUPMAX,SEED,RESETSEED
        INTEGER MUPDATE,IPRINT(2),MP,LP,I,FALSEIT,TAKESTEP
        INTEGER SW,INFO,GMAX_SW,RMS_TRAN_SW,J1   
        PARAMETER(NMUPMAX=4)
        REAL*8 GRAD(NDIM),XINIB(NDIM)
        REAL*8 POT,RMS,EPS,POT2
        REAL*8 WORK(NDIM*(2*NMUPMAX+1)+2*NMUPMAX)
        REAL*8 XTOL,GMAX,GTOL,STPMIN,STPMAX
        REAL*8 DIAG(NDIM),dx,RMS_FAC
        REAL*8 RADIUS,ASTEP,STEP,GTOL_INPUT
        REAL*8 XINIB_ORG(NDIM),DGUESS,RND,eng_temp,eng
        COMMON /accuracy4/ dx,RADIUS
        COMMON /bmin1/ DGUESS,ASTEP,STEP,GTOL_INPUT
        COMMON /bmin2/ MAXI,MUPDATE,TAKESTEP,RESETSEED,FALSEIT
        COMMON /bmin3/ GMAX,GMAX_SW
        COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
        COMMON /cgsw/ SW
        COMMON /DEBUG1/ INFO
        COMMON /RMS1/ RMS
        COMMON /RMS3/ RMS_FAC
        COMMON /RMS4/ RMS_TRAN_SW
        LOGICAL CFLAG,EVAP
        EXTERNAL LB2
        SW=0
        IFLAG=0
        XTOL=2.2204460492503D-16
        IPRINT(1)=-1
        IPRINT(2)=0 
        MP=6
        LP=6
        STPMIN=1.D-20
        STPMAX=1.D+20 
        CALL POTENTIAL(NDIM,XINIB,GRAD,POT,RMS,.TRUE.,EVAP)
C        IF(POT.gt.-99999999.D0.and.POT.lt.999999999.D0)then
C        else
C          write(*,*) "BEFORE ITER,NAN IN potential!!EREAL=",POT,RMS
C          do J1=1,ndim
C            write(*,*) J1,",XINI=",XINIB(J1),",GRAD=",GRAD(J1)
C          ENDDO
C          PAUSE
C        ENDIF
        DO ITER=1,MAXI
          IF (GMAX_SW.EQ.1) THEN
            IF (RMS.LE.GMAX) THEN
               EPS=RMS
             ELSE
               EPS=GMAX
             ENDIF
          ELSE
            EPS=GMAX
          ENDIF
          IF(RMS_TRAN_SW.eq.1)then
C            EPS=EPS+(EPS*RMS_FAC)*(-1.D0**(INT((RND(SEED)*10.D0))))
            EPS=EPS-RMS_FAC*EPS
            RMS_TRAN_SW=0
            RMS=EPS
            GMAX=EPS
C      write(*,*) ITER,",RMS TRAN",EPS,",before=",eng_temp,",after=",pot
          ENDIF
          eng_temp=pot
          CALL LBFGS(NDIM,MUPDATE,XINIB,POT,GRAD,
     &.FALSE.,DIAG,IPRINT,EPS,XTOL,WORK,IFLAG,DGUESS)
          CALL POTENTIAL(NDIM,XINIB,GRAD,POT,RMS,.TRUE.,EVAP)
C         IF(POT.gt.-99999999.D0.and.POT.lt.999999999.D0)then
C          else
C            write(*,*) ITER,"ITER,NAN IN potential!!EREAL=",POT,RMS
C            do J1=1,ndim
C              write(*,*) J1,",XINI=",XINIB(J1),",GRAD=",GRAD(J1)
C            ENDDO
C          PAUSE
C          ENDIF
          IF (IFLAG.EQ.0) THEN
            CFLAG=.TRUE.
            RETURN
          ENDIF
          IF (IFLAG.LT.0) THEN
C            WRITE(*,*) ITER,",WARNING - after LBFGS IFLAG=",IFLAG
            FALSEIT=FALSEIT+1
            IF(RESETSEED.EQ.1) THEN
              CALL dim_gen(XINIB,NDIM,SEED)
            ENDIF
            IF(TAKESTEP.EQ.1) THEN
              CALL tran(XINIB_ORG,XINIB,NDIM,SEED)
            ELSE IF(RESETSEED.NE.1.AND.TAKESTEP.NE.1)THEN
              RETURN
            ENDIF
            CFLAG=.FALSE.
            RETURN
          ENDIF
C         CALL FLUSH(6)
        ENDDO
        CFLAG=.FALSE. 
        RETURN
        END
