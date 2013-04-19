      FUNCTION dbrent(n,pcom,xicom,ax,bx,cx,f,df,tol,xmin)
      implicit none
      INTEGER ITMAX,n,iter
      REAL*8 dbrent,ax,bx,cx,tol,xmin,df,f,ZEPS,pcom(n),xicom(n)
      EXTERNAL df,f
      PARAMETER (ITMAX=20000,ZEPS=1.0e-10)
      REAL*8 a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2,u,u1,
     *u2,v,w,x,xm
      LOGICAL ok1,ok2,EVAP
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
*********************************
C     write(*,*) "min a=",a,",max b=",b,",mid point c=",c
      w=v
      x=v
      e=0.
      fx=f(n,pcom,xicom,x,EVAP)
      if(EVAP)return
      fv=fx
      fw=fx
      dx=df(n,pcom,xicom,x)
      dv=dx
      dw=dx
      do 11 iter=1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*dabs(x)+ZEPS
        tol2=2.*tol1
        if(dabs(x-xm).le.(tol2-.5*(b-a))) goto 3
        if(dabs(e).gt.tol1) then
          d1=2.*(b-a)
          d2=d1
          if(dw.ne.dx) d1=(w-x)*dx/(dx-dw)
          if(dv.ne.dx) d2=(v-x)*dx/(dx-dv)
          u1=x+d1
          u2=x+d2
          ok1=((a-u1)*(u1-b).gt.0.).and.(dx*d1.le.0.)
          ok2=((a-u2)*(u2-b).gt.0.).and.(dx*d2.le.0.)
          olde=e
          e=d
          if(.not.(ok1.or.ok2))then
            goto 1
          else if (ok1.and.ok2)then
            if(dabs(d1).lt.dabs(d2))then
              d=d1
            else
              d=d2
            endif
          else if (ok1)then
            d=d1
          else
            d=d2
          endif
          if(dabs(d).gt.dabs(0.5*olde))goto 1
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(dx.ge.0.) then
          e=a-x
        else
          e=b-x
        endif
        d=0.5*e
2       if(dabs(d).ge.tol1) then
          u=x+d
          fu=f(n,pcom,xicom,u,EVAP)
          if(EVAP)return
        else
          u=x+sign(tol1,d)
          fu=f(n,pcom,xicom,u,EVAP)
          if(EVAP)return
          if(fu.gt.fx)goto 3
        endif
        du=df(n,pcom,xicom,u)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          dv=dw
          w=x
          fw=fx
          dw=dx
          x=u
          fx=fu
          dx=du
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            dv=dw
            w=u
            fw=fu
            dw=du
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
            dv=du
          endif
        endif
11    continue
3     xmin=x
      dbrent=fx
      return
      END

