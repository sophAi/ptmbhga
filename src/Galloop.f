       program Galloop
       implicit none
       include "mpif.h"
       integer mpierror
       integer ens_num,seed,type,ndim,cent2,nproc,myid,ibin,node
       real*8 outx(6000),upot
       logical cent
       common nproc,myid
       common/engtype/ type
       common/cent/cent
*===========================MPI INITIAL===========================
       call mpi_init(mpierror)
       call mpi_comm_size(mpi_comm_world,nproc,mpierror)
       call mpi_comm_rank(mpi_comm_world,myid,mpierror)
       if(myid.eq.0)then
         write(*,*) "      =MPI-PTMBHPGA version 2.8="
         write(*,*) "        Last update:2004/9/10  "
       endif
       open(54,file="type.dat",status="old")
       read(54,*) type,cent2
       if(cent2.eq.1)then
         cent=.True.
       else
         cent=.False.
       endif
       close(54)
       if (type.eq.1.or.type.eq.2) then
         call gp_io(ndim,upot,outx,1)
       endif
       if (type.eq.3.or.type.eq.4) then
         call protein_io(ndim,upot,outx,1)
       endif
       if (type.eq.5) then
         call eam_io(ndim,upot,outx,1)
       endif
       if (type.eq.7) then
         call chiman_io(ndim,upot,outx,1)
       endif
       if (type.eq.8) then
         call chima3n_io(ndim,upot,outx,1)
       endif
       if (type.eq.9) then
         call qscff_io(ndim,upot,outx,1)
       endif
       call file_io(ibin,ens_num)
*==================================================================
       seed=INT(MPI_WTIME()*1000000.D0*2.D0+1.D0)+myid*13
       node=nproc
       call useseed(node,ndim,ibin,ens_num,seed,outx,upot)
       stop
       call mpi_barrier(mpi_comm_world,mpierror)
       write(*,*)"mpierror in gallop",mpierror,myid
       if(myid.eq.0)write(*,*) "out of seed",mpierror    
       if (type.eq.1.or.type.eq.2) then
         call gp_io(ndim,upot,outx,0)
       endif
       if (type.eq.3.or.type.eq.4) then
         call protein_io(ndim,upot,outx,0)
       endif
       if (type.eq.5) then
         call eam_io(ndim,upot,outx,0)
C  INPUT MOMENT OUTPUT
       endif
       if (type.eq.7) then
         call chiman_io(ndim,upot,outx,0)
       endif
       if (type.eq.8) then
         call chima3n_io(ndim,upot,outx,0)
       endif
       if (type.eq.9) then
         call qscff_io(ndim,upot,outx,0)
       endif
       call mpi_finalize(mpierror)
       stop
       end
   
       subroutine conj_min(xinib,min_pot,ndim,seed,cflag,EVAP)
       implicit none
       integer ndim,seed,s,sw
       real*8 xinib(ndim),min_pot
       real*8 FTOL,r,y
       logical EVAP,cflag
       common/accuracy6/ FTOL 
       common/cgsw/ sw
       EVAP=.false.
       cflag=.false.
       sw=11
152    call frprmn(xinib,ndim,FTOL,y,EVAP,cflag)
       if(EVAP)then
C         write(*,*) "EVAP!RETURN BACK"
         call dim_gen(xinib,ndim,seed)
         go to 152
       endif
*************************************
       min_pot=y
       return
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
         if(dabs(pot(np(1))-pot(np(a))).le.L) then
           same=same+1
         endif
       enddo
       return
       end
      
       subroutine sw_ene_sort(pot,same,ens_num,np)
       implicit none
       integer same,ens_num,iwork
       integer np(100)
       integer a,b
       real*8 pot(100),L
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
         if(dabs(pot(np(1))-pot(np(a))).le.L) then
           same=same+1
         endif
       enddo
       return
       end

 
       subroutine fitness(pot,par_num,ens_num,fit,fitg,np)
       implicit none
       integer i,ens_num   
       integer par_num,np(ens_num)
       real*8 pot(ens_num),fit(ens_num),fitg(ens_num+1)
       real*8 max_pot,min_pot,di_pot,sum_fit
       sum_fit=0.D0
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
       fitg(1)=fit(1)
       do i=1,par_num
         fitg(i+1)=fitg(i)+fit(i+1)
       enddo
       return
       end
    
       subroutine sw_fitness(pot,par_num,ens_num,fit,fitg,np)
       implicit none
       integer i,ens_num
       integer par_num,np(100)
       real*8 pot(100),fit(100),fitg(101)
       real*8 max_pot,min_pot,di_pot,sum_fit
       sum_fit=0.D0
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
       fitg(1)=fit(1)
       do i=1,par_num
         fitg(i+1)=fitg(i)+fit(i+1)
       enddo
       return
       end

       
       function selec_par(par_num,ens_num,fitg,np,seed)
       implicit none
       integer seed,par_num,selec_par,ens_num
       integer np(ens_num)
       real*8 fitg(ens_num+1)
       real*8 p,rnd
       integer i
       p=rnd(seed)
       do i=1,par_num
         if(p.lt.fitg(i)) then
           selec_par=np(i)
           goto 20
         endif
       enddo
       write(*,*) "ERROR_0,P=",P
C       do i=1,par_num
C         write(*,*) i,",fitg=",fitg(i),",p=",p
C       enddo
       selec_par=1       
20     return
       end

       function sw_selec_par(par_num,ens_num,fitg,np,seed)
       implicit none
       integer seed,par_num,sw_selec_par,ens_num
       integer np(100)
       real*8 fitg(101)
       real*8 p,rnd
       integer i
       p=rnd(seed)
       do i=1,par_num
         if(p.lt.fitg(i)) then
           sw_selec_par=np(i)
           goto 20
         endif
       enddo
       write(*,*) "ERROR_0,P=",P
       sw_selec_par=1
20     return
       end


C      1 is for low occupation number/2 is for high occupation number
C       function selec_com(i,ndim,gbin,guide_alli,gbin_all,io,seed)
C       implicit none
C       integer select_com,gbin,guide_all(gbin*gbin*gbin),io,gbin_all
C       integer FIX_GFAC
C       real*8 GFAC,SCALEFAC,TEMP,P,sel,no_sel
C       common/TEMP1/ SCALEFAC,TEMP
C       common/guide2/ GFAC,FIX_GFAC
C       p=rnd(seed)
C       if (io.eq.1) then
           
           
       subroutine inversion(xoya,xko,eff_co,seed)
       implicit none
       integer eff_co,seed
       real*8 xoya(eff_co),xko(eff_co),rnd,RADIUS,NONE,DG,GFAC
       real*8 SCALEFAC,TEMP,p
       integer i,r,s,dim3,pdb_rec,FIX_GFAC,DG_SW
       logical sel
       common/accuracy3/ RADIUS,NONE
       common/pdbrec/ pdb_rec,dim3
       common/guide1/ DG,DG_SW
       common/guide2/GFAC,FIX_GFAC
       common/TEMP1/ SCALEFAC,TEMP
C       r=mod(rndi(seed),eff_co)
C       s=mod(rndi(seed),eff_co)
       if(DG_SW.eq.1)then
15         r=int(rnd(seed)*dble(eff_co))+1
           s=int(rnd(seed)*dble(eff_co))+1
           if (r.gt.eff_co.or.r.le.0.or.s.gt.eff_co.or.s.le.0) goto 15
         do while(r.eq.s)
           s=int(rnd(seed)*dble(eff_co))+1
         enddo
         do i=1,eff_co
           xko(i)=xoya(i)
         enddo
         xko(r)=xoya(s)
         xko(s)=xoya(r)
         return
       endif
       p=rnd(seed)
       r=1
       s=2
       if(FIX_GFAC.eq.1) GFAC=TEMP
       do while(r.eq.s)
         if (dim3.ne.1)then
           r=int(rnd(seed)*dble(eff_co))+1
           s=int(rnd(seed)*dble(eff_co))+1
         else
           sel=.false.
           do while(.not.sel)
             r=int(rnd(seed)*dble(eff_co))+1
           enddo
16         sel=.false.
           do while(.not.sel)
17           s=int(rnd(seed)*dble(eff_co))+1
           enddo
         endif
       enddo
       do i=1,eff_co
         xko(i)=xoya(i)
       enddo
       xko(r)=xoya(s)
       xko(s)=xoya(r)
       return
       end
      
       subroutine arithmetic(xoya1,xoya2,xko,eff_co)
       implicit none
       integer eff_co
       real*8 xoya1(eff_co),xoya2(eff_co),xko(eff_co)
       integer i
       do i=1,eff_co
         xko(i)=0.5D0*(xoya1(i)+xoya2(i))
       enddo
       return
       end
       
       subroutine geometic(xoya1,xoya2,xko,eff_co)
       implicit none
       integer eff_co
       real*8 xoya1(eff_co),xoya2(eff_co),xko(eff_co)
       real*8 temp
       integer i 
       do i=1,eff_co
         temp=dabs(xoya1(i)*xoya2(i))
         xko(i)=dsqrt(temp)
       enddo
       return
       end 

       subroutine crossing(xoya1,xoya2,xko1,xko2,eff_co,seed)
       implicit none
       integer pdb_rec,dim3,FIX_GFAC,DG_SW
       integer eff_co
       real*8 xoya1(eff_co),xoya2(eff_co),xko1(eff_co),xko2(eff_co)
       real*8 p,rnd,SCALEFAC,TEMP,RADIUS,NONE,DG,GFAC
       integer i,seed
       common /accuracy3/ RADIUS,NONE
       common /TEMP1/ SCALEFAC,TEMP
       common /pdbrec/ pdb_rec,dim3
       common /guide1/ DG,DG_SW
       common /guide2/ GFAC,FIX_GFAC
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
       return
       end

       subroutine twopoint(xoya1,xoya2,xko1,xko2,eff_co,seed)
       implicit none
       integer pdb_rec,dim3,FIX_GFAC,DG_SW
       integer i,seed,s,rndi,j,k,eff_co
       real*8 xoya1(eff_co),xoya2(eff_co),xko1(eff_co)
       real*8 xko2(eff_co),temp_x(eff_co+eff_co),rnd
       real*8 SCALEFAC,TEMP,RADIUS,NONE,DG,GFAC
       common /accuracy3/ RADIUS,NONE
       common /TEMP1/ SCALEFAC,TEMP
       common /guide1/ DG,DG_SW
       common /guide2/GFAC,FIX_GFAC
       common /pdbrec/ pdb_rec,dim3
C       s=mod(rndi(seed),eff_co)
15     s=int(rnd(seed)*dble(eff_co))+1
       if(s.gt.eff_co.or.s.le.0)goto 15
       j=s-1
       k=eff_co-j          !3n-(s-1)
       do i=1,j
         xko1(k+i)=xoya1(i)
         xko2(k+i)=xoya2(i)
       enddo
       do i=1,k
         xko1(i)=xoya2(j+i)
         xko2(i)=xoya1(j+i)
       enddo
       return
       end                                           

C       subroutine crossover(xoya1,xoya2,xko1,xko2,eff_co)
C       implicit none
C       integer eff_co.seed
C       real*8 xoya1(eff_co),xoya2(eff_co),xko1(eff_co),
C     &xko2(eff_co),rnd
C       integer i,s
C15     s=int(rnd(seed)*dble(eff_co))+1
C       if(s.gt.eff_co.or.s.le.0)goto 15 
C       do i=1,eff_co
C         xko1(i)=0.5D0*(xoya1(i)+xoya2(i))
C       enddo
C       return
C       end

 
       subroutine regenerate(xoya,xko,eff_co,seed)
       implicit none
       integer pdb_rec,dim3,FIX_GFAC,j1,j2,j3,j4
       integer eff_co,seed,s,rndi,DG_SW
       real*8 xoya(eff_co),xko(eff_co),rnd,SCALEFAC,TEMP,RADIUS,NONE,
     &GFAC,DG,p,r,exp1
       logical sel
       common/accuracy1/ r
       common/accuracy3/ RADIUS,NONE
       common/TEMP1/ SCALEFAC,TEMP
       common/pdbrec/ pdb_rec,dim3
       common/guide1/ DG,DG_SW
       common/guide2/ GFAC,FIX_GFAC 
       p=rnd(seed)
       if(FIX_GFAC.eq.1) GFAC=TEMP
       if(dim3.eq.0) then        !for  3N
15       s=int(rnd(seed)*dble(eff_co/3))+1
         if (s.gt.(eff_co/3).or.s.le.0) goto 15
         do j1=1,eff_co/3
           j2=j1*3
           if(j1.le.s)then
             sel=.false.
             j4=j3
           else
             xko(j2-2)=xoya(j2-2)
             xko(j2-1)=xoya(j2-1)
             xko(j2)=xoya(j2)
           endif
         enddo
       else if(dim3.eq.1)then      ! for N
18       s=int(rnd(seed)*dble(eff_co))+1
         if (s.gt.eff_co.or.s.le.0)goto 18
         do j1=1,eff_co
           if(j1.ge.s)then
             sel=.false.
             j4=j3
           else
             xko(j1)=xoya(j1)
           endif        
         enddo
       endif
       return
       end

       subroutine dec_poi(xoya,xko,eff_co,seed)
       implicit double precision (a-h,o-z)
       integer eff_co
       dimension xoya(eff_co),xko(eff_co),temp(eff_co)
       do k=1,eff_co
         temp(k)=xoya(k)
       enddo
       igo_en=eff_co-1
10     igo=0
       do i=1,eff_co-1
         if (dabs(temp(i+1)).gt.dabs(temp(i))) then
           if (dabs(temp(i+1)).eq.temp(i+1)) then
             xko(i+1)=dabs(temp(i))
           else
             xko(i+1)=-dabs(temp(i))
           endif
           if (dabs(temp(i)).eq.temp(i)) then
             xko(i)=dabs(temp(i+1))
           else
             xko(i)=-dabs(temp(i+1))
           endif
         else
           igo=igo+1
         endif
       enddo
       do j=1,eff_co
         temp(i)=xko(i)
       enddo
       if (igo.lt.igo_en) goto 10
       return
       end

       subroutine sma_poi(xoya,xko,eff_co,seed)
       implicit double precision (a-h,o-z)
       integer eff_co
       dimension xoya(eff_co),xko(eff_co),temp(eff_co)
       do k=1,eff_co
         temp(k)=xoya(k)
       enddo
       do i=1,eff_co 
10       if (dabs(temp(i)).gt.1.D0) then
         temp(i)=temp(i)/10.D0
         goto 10
         endif
       enddo
       do j=1,eff_co
         xko(i)=temp(i)
       enddo
       return
       end


       subroutine dih_tran(xold,xnew,ndim,seed)
       implicit none
       integer pdb_rec,dim3,FIX_GFAC,DG_SW
       integer ndim,seed,n,i,j,k,jmax
       real*8 xnew(ndim),vt(2000),vtr(2000),b_x,b_y,b_z,xold(ndim)
       real*8 d11,d12,d13,d21,d22,d23,vmax,vmin,d_x,d_y,d_z,rnd
       real*8 mag_cro_ab,b_mag,theta,pei,dguess,gtol,gtol_input
       real*8 step,astep,sqr_abc,SCALEFAC,TEMP,RADIUS,NONE,DG,GFAC
       parameter(pei=3.14159D0,astep=0.4D0,step=0.5D0)
       common /vt1/ vt,vtr
       common /vt2/ jmax
       common /accuracy3/ RADIUS,NONE
       common /TEMP1/ SCALEFAC,TEMP
       common /pdbrec/ pdb_rec,dim3
       common /guide1/ DG,DG_SW
       common /guide2/ GFAC,FIX_GFAC
C       common /bmin1/ dguess,astep,step,gtol_input
       vmax=-1.0D6
       vmin=1.0D6
       do j=1,ndim
         xnew(j)=xold(j)
       enddo
       do j=1,ndim/3-2
         if (vt(j).gt.vmax) then
           vmax=vt(j)
           jmax=j
         endif
         if (vt(j).lt.vmin) vmin=vt(j)
       enddo
       do i=1,ndim/3-2
         k=3*i
C         write(*,*) ,i,VT(i)
C         if (dexp(vmin).lt.dexp(vt(i))) then
C        if ((vt(i).gt.astep*vmin).and.(i.eq.jmax)) then
C        if (i.eq.jmax)then
C           theta=pei*rnd(seed)/2.D0
          if (VT(i).gt.0.36D0*vmin.and.rnd(seed).gt.dexp(vt(i)))then
            theta=pei*rnd(seed)*((-1.D0)**(int(rnd(seed))*100))
C           theta=pei*((-1.D0)**(int(rnd(seed))*100))/4.D0
C            write(*,*) "MOVE!!",i,vmin*rnd(seed),vt(i)
            d11=xnew(k+1)-xnew(k-2)
            d12=xnew(k+2)-xnew(k-1)
            d13=xnew(k+3)-xnew(k)
            d21=xnew(k+4)-xnew(k+1)
            d22=xnew(k+5)-xnew(k+2)
            d23=xnew(k+6)-xnew(k+3)
            d_x=d21*d13-d23*d12
            d_y=d23*d11-d21*d13
            d_z=d21*d12-d22*d11
            mag_cro_ab=d_x**2+d_y**2+d_z**2
            sqr_abc=dsqrt(mag_cro_ab)
            b_x=xnew(k+1)-xnew(k+4)/2.D0-xnew(k-2)/2.D0
            b_y=xnew(k+2)-xnew(k+5)/2.D0-xnew(k-1)/2.D0
            b_z=xnew(k+3)-xnew(k+6)/2.D0-xnew(k)/2.D0
            b_mag=dsqrt(b_x**2.+b_y**2.+b_z**2.)
            xnew(k+1)=xnew(k+1)+b_mag*dsin(theta)*d_x/sqr_abc+
     &b_x*(dcos(theta)-1)
            xnew(k+2)=xnew(k+2)+b_mag*dsin(theta)*d_y/sqr_abc+
     &b_y*(dcos(theta)-1)
            xnew(k+3)=xnew(k+3)+b_mag*dsin(theta)*d_z/sqr_abc+
     &b_z*(dcos(theta)-1)
C           xnew(k+1)=dcos(theta)*b_x-b_mag*dsin(theta)*d_x/sqr_abc
C           xnew(k+2)=dcos(theta)*b_y-b_mag*dsin(theta)*d_y/sqr_abc
C           xnew(k+3)=dcos(theta)*b_z-b_mag*dsin(theta)*d_z/sqr_abc
          endif
        enddo
C        xnew(1)=xnew(1)+step*(rnd(seed)-0.5D0)*2.0D0
C        xnew(2)=xnew(2)+step*(rnd(seed)-0.5D0)*2.0D0
C        xnew(3)=xnew(3)+step*(rnd(seed)-0.5D0)*2.0D0
         xnew(ndim-2)=xnew(ndim-2)+step*(rnd(seed))
         xnew(ndim-1)=xnew(ndim-1)+step*(rnd(seed))
         xnew(ndim)=xnew(ndim)+step*(rnd(seed))
       return
       end

       subroutine alpha_beta(dih,ndim,seed)
       implicit none
       integer pdb_rec,dim3,FIX_GFAC,DG_SW
       integer ndim,seed,i,n,jmax
       real*8 dih(ndim),pei,p,theta,rnd,vt(2000),vtr(2000),DG,GFAC
       real*8 a,b,c,d,e,f,g,h,dih_step,factor,SCALEFAC,TEMP,RADIUS,NONE
       parameter (pei=3.14159D0)
       common /CHIMA4/ dih_step,a,b,c,d,e,f,g,h
       common /dih1/ factor
       common /vt1/vt,vtr
       common /vt2/ jmax
       common /accuracy3/ RADIUS,NONE
       common /TEMP1/ SCALEFAC,TEMP
       common /pdbrec/ pdb_rec,dim3
       common /guide1/ DG,DG_SW
       common /guide2/ GFAC,FIX_GFAC
       p=rnd(seed)
       if (p.lt.a/(a+b+c+d+e+f+g+h))then
         do i=1,ndim
C           if (vt(i).gt.0.D0) then
             dih(i)=0.0D0
C           endif
         enddo
       else if(p.lt.(a+b)/(a+b+c+d+e+f+g+h))then
C         n=int(ndim/2)
         n=int(dble(ndim-2)*rnd(seed))+2
         dih(n-2)=-12.5D0*pei/180.D0
         dih(n-1)=97.5D0*pei/180.D0
         dih(n)=-0.999D0*pei
         dih(n+1)=-97.5D0*pei/180.D0
         dih(n+2)=12.5D0*pei/180.D0
       else if(p.lt.(a+b+c)/(a+b+c+d+e+f+g+h))then
C         n=int(ndim/2)
         n=int(dble(ndim-2)*rnd(seed))+2
         dih(n-2)=12.5D0*pei/180.D0
         dih(n-1)=-97.5D0*pei/180.D0
         dih(n)=-0.999D0*pei
         dih(n+1)=97.5D0*pei/180.D0
         dih(n+2)=-12.5D0*pei/180.D0
       else if(p.lt.(a+b+c+d)/(a+b+c+d+e+f+g+h))then
         do i=1,ndim
           theta=1.91986214D0+
     &(2.8274334D0-1.91986214D0)*rnd(seed)/1.D0
           dih(i)=theta
         enddo
       else if(p.lt.(a+b+c+d+e)/(a+b+c+d+e+f+g+h))then
         do i=1,ndim
           theta=-1.91986214D0+
     &(1.91986214D0-2.8274334D0)*rnd(seed)/1.D0
           dih(i)=theta
         enddo                    
       else if(p.lt.(a+b+c+d+e+f)/(a+b+c+d+e+f+g+h))then
C         n=int(ndim/2)
         n=int(dble(ndim-2)*rnd(seed))+2
         dih(n-2)=-56.25D0*pei/180.D0
         dih(n-1)=110.D0*pei/180.D0
         dih(n)=-110.D0*pei/180.D0
         dih(n+1)=56.25D0*pei/180.D0
       else if(p.lt.(a+b+c+d+e+f+g)/(a+b+c+d+e+f+g+h))then
C         n=int(ndim/2
         n=int(dble(ndim-2)*rnd(seed))+2
         dih(n-2)=56.25D0*pei/180.D0
         dih(n-1)=-110.D0*pei/180.D0
         dih(n)=110.D0*pei/180.D0
         dih(n+1)=-56.25D0*pei/180.D0
       else        
         n=int(dble(ndim)*rnd(seed))
         do i=1,ndim
           if (rnd(seed).gt.dexp(vt(i)))then
C         dih(n)=dih(n)+(pei*rnd(seed)*(-1.D0**int(rnd(seed*10))))
             theta=dih_step*(-1.D0**((rnd(seed)*100.D0)))*pei/180.D0
             dih(i)=dih(i)+theta
             if(i.lt.ndim-2)then
               dih(i+2)=dih(i)-theta
             endif
          endif
        enddo
       endif
       return
       end
