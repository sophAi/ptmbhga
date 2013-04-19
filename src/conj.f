       subroutine conj_min(xinib,min_pot,ndim,seed,r,s)
       implicit none
       integer i,j,k,o,ndim,seed,s,sw
       real*8 xinib(ndim),xini(ndim),min_pot
       real*8 y,ltemp
       real*8 FTOL,r,rnd
       logical EVAP
       common/accuracy6/ FTOL 
       common/cgsw/ sw
       EVAP=.false.
       sw=11
152    call frprmn(xinib,ndim,FTOL,y,EVAP)
       if(EVAP)then
         if(s.eq.1)then
           return
         endif
         call eulerangles(xinib,ndim,seed,r)
         go to 152
       endif
*************************************
       min_pot=y
       end
       
        
       subroutine ene_sort(pot,same,ens_num,np)
       implicit none
       integer same,ens_num,iwork 
       integer np(ens_num)
       integer a,b         
       real*8 pot(ens_num),L
       common/accuracy2/ L
       do a=1,ens_num
         np(a)=a
       enddo
       do a=1,ens_num
         do b=a+1,ens_num
           if(pot(np(a)).gt.pot(np(b))) then
             iwork=np(a)
             np(a)=np(b)
             np(b)=iwork 
           endif
         enddo
       enddo
       same=0
       do a=1,ens_num
         if(abs(pot(np(1))-pot(np(a))).le.L) then
           same=same+1
         endif
       enddo
       end
       
       subroutine fitness(pot,par_num,ens_num,fit,fitg,np)
       integer i,ens_num   
       integer par_num,np(ens_num)
       real*8 pot(ens_num),fit(ens_num)
       real*8 fitg(0:500)
       real*8 max_pot,min_pot,di_pot,sum_fit
       sum_fit=0.D0
       fitg(0)=0.D0
       min_pot=pot(np(1))
       max_pot=pot(np(par_num))
       di_pot=max_pot-min_pot
       do i=1,par_num
         fit(i)=(max_pot-pot(np(i)))/di_pot
       enddo
       do i=1,par_num
         sum_fit=sum_fit+fit(i)
       enddo
       do i=1,par_num
         fit(i)=fit(i)/sum_fit
       enddo
       do i=0,par_num-1
         fitg(i+1)=fitg(i)+fit(i+1)
       enddo
       end
       
       function selec_par(par_num,ens_num,fitg,seed)
       implicit none
       integer seed,par_num,selec_par,ens_num
       real*8 fitg(0:500)
       real*8 p,rnd
       integer i
       p=rnd(seed)
       do i=1,par_num
         if(p.lt.fitg(i)) then
           selec_par=i
           goto 20
         endif
       enddo
       write(*,*) "ERROR_0"
       selec_par=1       
20     end
        
       subroutine inversion(xoya,xko,eff_co,seed)
       implicit none
       integer eff_co,seed
       real*8 xoya(eff_co),xko(eff_co)
       integer i,r,s,rndi,srnd
       r=mod(rndi(seed),eff_co)
       s=mod(rndi(seed),eff_co)
       do while(r.eq.s) 
         s=mod(rndi(seed),eff_co) 
       enddo
       do i=1,eff_co
         xko(i)=xoya(i)
       enddo
       xko(r)=xoya(s)
       xko(s)=xoya(r)
       end
      
       subroutine arithmetic(xoya1,xoya2,xko,eff_co)
       implicit none
       integer eff_co
       real*8 xoya1(eff_co),xoya2(eff_co),xko(eff_co)
       integer i
       do i=1,eff_co
         xko(i)=0.5D0*(xoya1(i)+xoya2(i))
       enddo
       end
       
       subroutine geometic(xoya1,xoya2,xko,eff_co)
       implicit none
       integer eff_co
       real*8 xoya1(eff_co),xoya2(eff_co),xko(eff_co)
       real*8 temp
       integer i 
       do i=1,eff_co
         temp=abs(xoya1(i)*xoya2(i))
         xko(i)=sqrt(temp)
       enddo
       end 

       subroutine crossing(xoya1,xoya2,xko1,xko2,eff_co,seed)
       implicit none
       integer eff_co
       real*8 xoya1(eff_co),xoya2(eff_co),xko1(eff_co),xko2(eff_co)
       real*8 p,rnd
       integer i,seed
       do i=1,eff_co
         p=rnd(seed)
         if(p.lt.0.5D0) then
           xko1(i)=xoya1(i)
           xko2(i)=xoya2(i)
         else
           xko1(i)=xoya2(i)
           xko2(i)=xoya1(i)
         endif
       enddo 
       end

       subroutine twopoint(xoya1,xoya2,xko1,xko2,eff_co,seed)
       implicit none
       integer i,seed,s,rndi,j,eff_co
       real*8 xoya1(eff_co),xoya2(eff_co),xko1(eff_co)
       real*8 xko2(eff_co),temp(2*eff_co)
       s=mod(rndi(seed),eff_co)
       do i=1,eff_co
         temp(i)=xoya1(i)
         temp(i+eff_co)=xoya2(i)
       enddo
       j=1
       do i=1,eff_co
         xko1(i)=temp(s+i)
         if ((s+i+eff_co).le.(eff_co*2))then
           xko2(i)=temp(s+i+eff_co)
         else
           xko2(i)=temp(j)
           j=j+1
         endif
       enddo
       end                                           
