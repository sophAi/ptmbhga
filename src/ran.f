       function ran(seed)
       implicit none
       integer seed,IA,IM,IQ,IR,NTAB,NDIV
       real*8 ran,AM,EPS,RNMX
       parameter (IA=16807,IM=2147483647,AM=1.D0/IM,IQ=127773)
       parameter (IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2D-7)
       parameter (RNMX=1.D0-EPS)
       integer j,k,iv(NTAB),iy
       save iv,iy
       if(seed.le.0.or.iy.eq.0) then
         seed=max(abs(seed),1)
         do j=NTAB+8,1,-1
           k=seed/IQ
           seed=IA*(seed-k*IQ)-IR*k
           if(seed.lt.0) seed=seed+IM
           if(j.le.NTAB) iv(j)=seed
         enddo
         iy=iv(1)
       endif
       k=seed/IQ
       seed=IA*(seed-k*IQ)-IR*k
       if(seed.lt.0) seed=seed+IM
       j=1+iy/NDIV
       iy=iv(j)
       iv(j)=seed
       ran=dmin1(AM*iy,RNMX)
       return
       end
       function rndi(seed)
       implicit none
       integer seed,rndi
       real*8 ran
       rndi=1+int(1000000.D0*ran(seed))
       end
