      SUBROUTINE mnbrak(n,pcom,xicom,ax,bx,cx,fa,fb,fc,f,EVAP)
      IMPLICIT NONE
      INTEGER n
      REAL*8 ax,bx,cx,fa,fb,fc,f,GOLD,GLIMIT,TINY
      EXTERNAL f
      PARAMETER (GOLD=1.618034D0,GLIMIT=100.D0, TINY=1.e-20)
      REAL*8 dum,fu,q,r,u,ulim,pcom(n),xicom(n)
      LOGICAL EVAP
      fa=f(n,pcom,xicom,ax,EVAP)
      if(EVAP)return
      fb=f(n,pcom,xicom,bx,EVAP)
      if(EVAP)return
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      fc=f(n,pcom,xicom,cx,EVAP)
      if(EVAP)return
C      write(*,*) "fa=",fa,",fb=",fb,",fc=",fc
C      write(*,*) "ax=",ax,",bx=",bx,",cx=",cx
C      pause "mnbark pause"
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))

        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
          fu=f(n,pcom,xicom,u,EVAP)
          if(EVAP)return
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
          fu=f(n,pcom,xicom,u,EVAP)
          if(EVAP)return
        else if((cx-u)*(u-ulim).gt.0.)then
          fu=f(n,pcom,xicom,u,EVAP)
          if(EVAP)return
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)

            fb=fc
            fc=fu
            fu=f(n,pcom,xicom,u,EVAP)
           if(EVAP)return
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
          fu=f(n,pcom,xicom,u,EVAP)
          if(EVAP)return
        else
          u=cx+GOLD*(cx-bx)
          fu=f(n,pcom,xicom,u,EVAP)
          if(EVAP)return
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
C     write(*,*) "inverse statment fa.gt.fb"
      return
      END
