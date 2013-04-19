      FUNCTION df1dim(ncom,pcom,xicom,x)
      implicit none
      REAL*8 x
CU    USES dfunc
      INTEGER j,ncom
      REAL*8 df(ncom),pcom(ncom),xicom(ncom),xt(ncom),df1dim
      REAL*8 EREAL,RMS
      LOGICAL EVAP
C      COMMON /f1com/ pcom,xicom    
      do j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
      enddo
      call POTENTIAL(ncom,xt,df,EREAL,RMS,.true.,EVAP)
      df1dim=0.D0
      do j=1,ncom
        df1dim=df1dim+df(j)*xicom(j)
      enddo
      return
      END
