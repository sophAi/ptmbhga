        subroutine tran(xinib,xnew,ndim,seed)
        implicit none
        integer ndim,seed,i,j,k,jmax,l,FIX_GFAC,DG_SW
        real*8 xinib(ndim),rnd,pei,dist(ndim),theta,astep,SCALEFAC,TEMP
        real*8 vmax,vmin,dmax,phi,vt(2000),vtr(2000),DG,GFAC
        real*8 xnew(ndim),step,dummy,radius,non,dguess,gtol_input
        parameter(pei=3.141592654D0)
        common /accuracy3/ radius,non
        common /bmin1/ dguess,astep,step,gtol_input
        common /vt1/ vt,vtr
        common /vt2/ jmax
        common /vt3/ vmin,vmax,dmax
        common /TEMP1/ SCALEFAC,TEMP
        common /guide1/ DG,DG_SW
        common /guide2/ GFAC,FIX_GFAC
C        vmax=-1.0D6
C        vmin=1.0D6
C        dmax=-1.0D0
C        do j=1,ndim/3
C          write(*,*) "vtr",j,"=",vtr(j),",vt=",vt(j)
C        enddo
C        pause
C        do j=1,ndim/3
C          k=3*j    
C          dist(j)=dsqrt(xinib(k-2)**2+xinib(k-1)**2+xinib(k)**2)
C          if (dist(j).gt.dmax) dmax=dist(j)
C          if (vt(j).gt.vmax) then
C            vmax=vt(j)
C            jmax=j
C          endif
C          if (vt(j).lt.vmin) vmin=vt(j)
C        enddo
        do i=1,ndim/3
          k=3*i
c          write(*,*) "vtr(max)=",vtr(i),i
          if ((vt(i).gt.astep*vmin).and.(i.eq.jmax)) then
C          if(vtr(i).eq.1.D0.and.vt(i).gt.astep*vmin)then
C            write(*,*) "jmax=",i,jmax,",vt(i)=",vt(i),vtr(i)
C            pause
            theta=pei*rnd(seed)
C            phi=pei*2.0D0*rnd(seed)*((-1.D0)**(int(rnd(seed)*10)))
            phi=pei*rnd(seed)*2.D0
            xnew(k-2)=dmax*dsin(theta)*dcos(phi)
            xnew(k-1)=dmax*dsin(theta)*dsin(phi)
            xnew(k)=dmax*dcos(theta)
          else 
            xnew(k-2)=xinib(k-2)+step*(rnd(seed)-0.5D0)*2.D0      
            xnew(K-1)=xinib(k-1)+step*(rnd(seed)-0.5D0)*2.D0
            xnew(k)=xinib(k)+step*(rnd(seed)-0.5D0)*2.D0
            dummy=xnew(k-2)**2+xnew(k-1)**2+xnew(k)**2
C            dummy=dsqrt(dummy)
            if (dummy.gt.radius) then
              xnew(k-2)=xnew(k-2)*dsqrt(radius/dummy)
              xnew(k-1)=xnew(k-1)*dsqrt(radius/dummy)
              xnew(k)=xnew(k)*dsqrt(radius/dummy)
            endif
          endif
        enddo
        return
        end
 
        subroutine tran_n(xinib,xnew,ndim,seed,
     &RMS_SW)
        implicit none
        integer ndim,seed,j1,j2,j3,RMS_TRAN_SW
        real*8 xinib(ndim),xnew(ndim),rnd,astep,step,RMS_FAC
        real*8 gtol_input,dguess
        logical RMS_SW
        common/bmin1/ dguess,astep,step,gtol_input 
C        common /bmin3/ GMAX,GMAX_SW
C        common /RMS1/ RMS
        common /RMS3/ RMS_FAC
        common /RMS4/ RMS_TRAN_SW
14      j2=int(dble(ndim)*rnd(seed))+1
        if(j2.gt.ndim.or.j2.le.0)then
          goto 14
        endif
        do j1=1,ndim
          xnew(j1)=xinib(j1)
        enddo
        do j1=1,ndim
          if (j2.ne.0)then
15          j3=int(dble(ndim)*rnd(seed))+1
            if(j3.gt.ndim.or.j3.le.0)goto 15
            xnew(j3)=xinib(j3)+step*(2.D0*rnd(seed)-1.D0)
            j2=j2-1
          else
            goto 12
          endif
        enddo
12      if(RMS_SW) RMS_TRAN_SW=1
        return
        end
