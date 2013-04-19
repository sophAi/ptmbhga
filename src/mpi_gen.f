       subroutine mpi_gen(xini,pot,ens_num,ndim,np,seed,fit,xw,ex
     &,counter)
       implicit none
       integer ens_num,seed,ndim,np(ens_num),type
       real*8 xini(ens_num,ndim),pot(ens_num),fit(ens_num)
       integer counter,par_num
       integer selec_par,no,no1,no2,i,j,k,ex(ens_num*2),sw_eular1
       real*8 xoya(ndim),xko(ndim),xoya1(ndim),xoya2(ndim)
       real*8 xko1(ndim),xko2(ndim),eular_fac1
       real*8 p,rnd,fitg(ens_num+1),exp1
       real*8 aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii,xw(ens_num*ndim)
       logical evap
       common/engtype/ type
       common/mpigen1/ aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii
       common/mpigen2/par_num
       common/eularangle1/ sw_eular1
       common/eularangle3/ eular_fac1
       counter=ens_num
       call fitness(pot,par_num,ens_num,fit,fitg,np)
       k=0
C       write(*,*) gbin_all,gbin
C       do i=1,gbin*gbin*gbin
C        if (dble(guide_all(i))/dble(gbin_all).ne.0)then
C        write(*,*) i,",guide/gbin=",dble(guide_all(i))/dble(gbin_all),
C     &gbin_all
C        exp1=-dble(guide_all(i))/(dble(gbin_all)*0.001)
C        exp1=dexp(exp1)
C        write(*,*) "exp=",exp1
C        write(*,*) dble(gbin_all)*rnd(seed)
C        endif
C       enddo
       do while(counter.ge.par_num)
         p=rnd(seed)
****************Inversion***************************
         if(p.lt.aaa) then
           no=np(selec_par(par_num,ens_num,fitg,np,seed))
           do i=1,ndim
             xoya(i)=xini(no,i)
           enddo
           if(sw_eular1.eq.1)then
             call eularangles_tran(xoya,ndim,seed,eular_fac1)
           endif
           call inversion(xoya,xko,ndim,seed)
           do i=1,ndim
             xw(k*ndim+i)=xko(i)
           enddo
           k=k+1
           ex(k*2-1)=no
           ex(k*2)=0
           counter=counter-1
***************Arithmetic**************************
         else if(p.lt.(aaa+bbb)) then
           no1=np(selec_par(par_num,ens_num,fitg,np,seed))
           no2=np(selec_par(par_num,ens_num,fitg,np,seed))
           do i=1,ndim
             xoya1(i)=xini(no1,i)
             xoya2(i)=xini(no2,i)
           enddo
           if (sw_eular1.eq.1)then
             call eularangles_tran(xoya1,ndim,seed,eular_fac1)
             call eularangles_tran(xoya2,ndim,seed,eular_fac1)
           endif
           call arithmetic(xoya1,xoya2,xko,ndim)
           do i=1,ndim
             xw(k*ndim+i)=xko(i)
           enddo
           k=k+1
           ex(k*2-1)=no1
           ex(k*2)=no2
           counter=counter-1
***************Geometic****************************
         else if(p.lt.(aaa+bbb+ccc)) then
           no1=np(selec_par(par_num,ens_num,fitg,np,seed))
           no2=np(selec_par(par_num,ens_num,fitg,np,seed))
           do i=1,ndim
             xoya1(i)=xini(no1,i)
             xoya2(i)=xini(no2,i)
           enddo
           if (sw_eular1.eq.1)then
             call eularangles_tran(xoya1,ndim,seed,eular_fac1)
             call eularangles_tran(xoya2,ndim,seed,eular_fac1)
           endif
           call geometic(xoya1,xoya2,xko,ndim)
           do i=1,ndim
             xw(k*ndim+i)=xko(i)
           enddo
           k=k+1
           ex(k*2-1)=no1
           ex(k*2)=no2
           counter=counter-1
***************Crossing****************************
         else if(p.lt.(aaa+bbb+ccc+ddd)) then
           no1=np(selec_par(par_num,ens_num,fitg,np,seed))
           no2=np(selec_par(par_num,ens_num,fitg,np,seed))
           do i=1,ndim
             xoya1(i)=xini(no1,i)
             xoya2(i)=xini(no2,i)
           enddo
           if (sw_eular1.eq.1)then
             call eularangles_tran(xoya1,ndim,seed,eular_fac1)
             call eularangles_tran(xoya2,ndim,seed,eular_fac1)
           endif
           call crossing(xoya1,xoya2,xko1,xko2,ndim,seed)
           do i=1,ndim
             xw(k*ndim+i)=xko1(i)
           enddo
           k=k+1
           ex(k*2-1)=no1
           ex(k*2)=no2
           do i=1,ndim
             xw(k*ndim+i)=xko2(i)
           enddo
           k=k+1
           ex(k*2-1)=no1
           ex(k*2)=no2
           counter=counter-2
***************2-Point Crossover********************
         else if (p.lt.(aaa+bbb+ccc+ddd+eee)) then
           no1=np(selec_par(par_num,ens_num,fitg,np,seed))
           no2=np(selec_par(par_num,ens_num,fitg,np,seed))
           do i=1,ndim
             xoya1(i)=xini(no1,i)
             xoya2(i)=xini(no2,i)
           enddo
           if (sw_eular1.eq.1)then
             call eularangles_tran(xoya1,ndim,seed,eular_fac1)
             call eularangles_tran(xoya2,ndim,seed,eular_fac1)
           endif
           call twopoint(xoya1,xoya2,xko1,xko2,ndim,seed)
           do i=1,ndim
             xw(k*ndim+i)=xko1(i)
           enddo
           k=k+1
           ex(k*2-1)=no1
           ex(k*2)=no2
           do i=1,ndim
             xw(k*ndim+i)=xko2(i)
           enddo
           k=k+1
           ex(k*2-1)=no1
           ex(k*2)=no2  
           counter=counter-2
***************mutation***************************
         else if (p.lt.(aaa+bbb+ccc+ddd+eee+fff)) then
           no=np(selec_par(par_num,ens_num,fitg,np,seed))
           do i=1,ndim
             xko(i)=xini(no,i)
           enddo
           call mutation(xko,xko2,ndim,seed)
           do i=1,ndim
             xw(k*ndim+i)=xko2(i)
           enddo
           k=k+1
           ex(k*2-1)=no
           ex(k*2)=0
           counter=counter-1
         else if (p.lt.(aaa+bbb+ccc+ddd+eee+fff+ggg)) then
           no=np(selec_par(par_num,ens_num,fitg,np,seed))
           do i=1,ndim
             xoya(i)=xini(no,i)
           enddo
           call dec_poi(xoya,xko,ndim,seed)
           do i=1,ndim
             xw(k*ndim+i)=xko(i)
           enddo
           k=k+1
           ex(k*2-1)=no
           ex(k*2)=0
           counter=counter-1
         else if (p.lt.(aaa+bbb+ccc+ddd+eee+fff+ggg+hhh)) then
           no=np(selec_par(par_num,ens_num,fitg,np,seed))
           do i=1,ndim
             xoya(i)=xini(no,i)
           enddo
           call sma_poi(xoya,xko,ndim,seed)
           do i=1,ndim
             xw(k*ndim+i)=xko(i)
           enddo
           k=k+1
           ex(k*2-1)=no
           ex(k*2)=0
           counter=counter-1
         else
           no=np(selec_par(par_num,ens_num,fitg,np,seed))
           do i=1,ndim
             xko(i)=xini(no,i)
           enddo
           call alpha_beta(xko,ndim,seed)
           do i=1,ndim
             xw(k*ndim+i)=xko(i)
           enddo
           k=k+1
           ex(k*2-1)=no
           ex(k*2)=0
           counter=counter-1  
         endif
       enddo
       end
