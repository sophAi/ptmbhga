      SUBROUTINE frprmn(p,n,ftol,fret,EVAP,cflag)
      implicit none
      INTEGER iter,n,ITMAX
      REAL*8 fret,ftol,p(n),EPS,EREAL
      PARAMETER (ITMAX=200,EPS=1.e-10)
      INTEGER its,j,k
      REAL*8 dgg,fp,gam,gg,g(n),h(n),xi(n)
      LOGICAL EVAP,cflag
      call POTENTIAL(n,p,xi,EREAL,.true.,EVAP)
      fp=EREAL
      if(EVAP)then
        return
      endif
      do 11 j=1,n
        g(j)=-xi(j)
        h(j)=g(j)
        xi(j)=h(j)
11    continue
      do 14 its=1,ITMAX
        iter=its
        call linmin(p,xi,n,fret,EVAP)
        if(EVAP)then
          return
        endif
        if(2.*abs(fret-fp).le.ftol*(abs(fret)+abs(fp)+EPS))then
          cflag=.true.
          return
        endif
        call POTENTIAL(n,p,xi,EREAL,.true.,EVAP)
        fp=EREAL
C        do k=1,n
C          write(*,*) xi(k),p(k),n
C        enddo
C        pause
        if (EVAP)then
          return
        endif
        gg=0.
        dgg=0.
        do 12 j=1,n
          gg=gg+g(j)**2
C         dgg=dgg+xi(j)**2
          dgg=dgg+(xi(j)+g(j))*xi(j)
12      continue
        if(gg.eq.0.)return
          gam=dgg/gg
        do 13 j=1,n
          g(j)=-xi(j)
          h(j)=g(j)+gam*h(j)
          xi(j)=h(j)
13      continue
14    continue
      return
      end
