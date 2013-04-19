       subroutine amoeba(p,y,ftol,eff_co,cflag)
       implicit none
       integer iter,eff_co
       integer ITMAX,i,ihi,ilo,inhi,j,m,n   
       real*8 p(eff_co+1,eff_co),y(eff_co+1),ftol,EREAL
       real*8 GRAD(eff_co),RMS
C      parameter (ITMAX=2**30-2)
       parameter (ITMAX=5000)
       real*8 rtol,sum,swap,ysave,ytry,psum(eff_co),amotry
       logical EVAP,cflag
       iter=0
1      do n=1,eff_co
         sum=0.D0
         do m=1,eff_co+1
           sum=sum+p(m,n)
         enddo
         psum(n)=sum
       enddo 
2      ilo=1
       if(y(1).gt.y(2)) then
         ihi=1
         inhi=2
       else
         ihi=2
         inhi=1
       endif
       do i=1,eff_co+1
         if(y(i).le.y(ilo)) ilo=i 
         if(y(i).gt.y(ihi)) then
           inhi=ihi
           ihi=i
         else if(y(i).gt.y(inhi)) then 
           if(i.ne.ihi) inhi=i 
         endif
       enddo  
       rtol=2.D0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
       if(rtol.lt.ftol) then
         swap=y(1)
         y(1)=y(ilo)
         y(ilo)=swap
         do n=1,eff_co
           swap=p(1,n)
           p(1,n)=p(ilo,n)
           p(ilo,n)=swap
         enddo
         cflag=.true.
         return
       endif
C       if(iter.ge.ITMAX) pause 'ITMAX exceeded amoba'
       iter=iter+2
       ytry=amotry(p,y,psum,ihi,-1.D0,eff_co)
       if(ytry.le.y(ilo)) then
         ytry=amotry(p,y,psum,ihi,2.D0,eff_co)   
       else if(ytry.ge.y(inhi)) then
         ysave=y(ihi)
         ytry=amotry(p,y,psum,ihi,0.5D0,eff_co)
         if(ytry.ge.ysave) then
           do i=1,eff_co+1
             if(i.ne.ilo) then
               do j=1,eff_co
                 psum(j)=0.5*(p(i,j)+p(ilo,j))
                 p(i,j)=psum(j)
               enddo
               call POTENTIAL(eff_co,psum,GRAD,EREAL,RMS,.false.
     &,EVAP)
               y(i)=EREAL
C               write(*,*) "EREAL IN SIMPLEX=",EREAL,eff_co
             endif
           enddo
           iter=iter+eff_co
           goto 1
         endif
       else
         iter=iter-1
       endif
       goto 2
       end

       function amotry(p,y,psum,ihi,fac,eff_co)
       implicit none 
       integer ihi,eff_co,j
       real*8 p(eff_co+1,eff_co),y(eff_co+1),psum(eff_co),fac
       real*8 fac1,fac2,ytry,ptry(eff_co),GRAD(eff_co),amotry
       real*8 EREAL,vt(eff_co),RMS
       logical EVAP
       fac1=(1.D0-fac)/eff_co
       fac2=fac1-fac
       do j=1,eff_co
         ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
       enddo
       call POTENTIAL(eff_co,ptry,GRAD,EREAL,RMS,.false.,EVAP)
       ytry=EREAL
       if(ytry.lt.y(ihi)) then   
         y(ihi)=ytry
         do j=1,eff_co
           psum(j)=psum(j)-p(ihi,j)+ptry(j)
           p(ihi,j)=ptry(j)
         enddo
       endif
       amotry=ytry
       return
       end
 
