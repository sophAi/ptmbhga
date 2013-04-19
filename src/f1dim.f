      FUNCTION f1dim(ncom,pcom,xicom,x,EVAP)
      implicit none
      REAL*8 f1dim,x,RMS
CU    USES func
      INTEGER j,ncom
      REAL*8 pcom(ncom),xicom(ncom),xt(ncom),GRAD(ncom)
      REAL*8 EREAL
      LOGICAL EVAP
C     COMMON /f1com/ pcom,xicom
      do j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
      enddo
      call POTENTIAL(ncom,xt,GRAD,EREAL,RMS,.false.,EVAP)
      f1dim=EREAL
      return
      END
