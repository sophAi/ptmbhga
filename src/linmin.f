      SUBROUTINE linmin(p,xi,n,fret,EVAP)
      IMPLICIT NONE
      INTEGER j,n
CU    USES brent,f1dim,mnbrak
      REAL*8 fret,p(n),xi(n),TOL,f1dim,df1dim
      PARAMETER (TOL=1.D-4)
      REAL*8 ax,bx,fa,fb,fx,xmin,xx,pcom(n),xicom(n),dbrent
C     COMMON /f1com/ pcom,xicom
      EXTERNAL f1dim,df1dim
      LOGICAL EVAP
      do 11 j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
11    continue
158   ax=0.D0
      xx=1.D0
      call mnbrak(n,pcom,xicom,ax,xx,bx,fa,fx,fb,f1dim,EVAP)
      if(EVAP)return
      fret=dbrent(n,pcom,xicom,ax,xx,bx,f1dim,df1dim,TOL,xmin)
      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
12    continue
      return
      END
