      subroutine mutation(xinib,xnew,ndim,seed)
      implicit none
      integer A,ndim,seed,j1,j2,j3,j4
     &,FIX_GFAC,dim3,pdb_rec,
     &DG_SW,r,s,rndi,MIN_S,MIN_T
      real*8 xinib(ndim),xnew(ndim),RADIUS,NONE,SCALEFAC,DG,
     &GFAC,TEMP
      real*8 p,exp1,rnd,esys,RMS,grad(ndim),xtemp(ndim),
     &min_pot,esys2
      logical sel,pass,EVAP,CFLAG
      common /alloy1/ A
      common /accuracy3/ RADIUS,NONE
      common /TEMP1/ SCALEFAC,TEMP
      COmmon /pdbrec/ pdb_rec,dim3
      common /guide1/ DG,DG_SW
      common /guide2/ GFAC,FIX_GFAC
C      if(DG_SW.eq.1)then
        min_pot=0.D0
        do J1=1,ndim
          xnew(J1)=xinib(J1)
          xtemp(J1)=xinib(J1)
        enddo
        call POTENTIAL(ndim,xtemp,grad,esys,RMS,.false.,EVAP)
        min_pot=esys
        do J1=1,A
          do J2=A+1,ndim/3
            xtemp(J2*3-2)=xinib(J1*3-2)
            xtemp(J2*3-1)=xinib(J1*3-1)
            xtemp(J2*3)=xinib(J1*3)
            xtemp(J1*3-2)=xinib(J2*3-2)
            xtemp(J1*3-1)=xinib(J2*3-1)
            xtemp(J1*3)=xinib(J2*3)
C            call POTENTIAL(ndim,xtemp,grad,esys,RMS,.false.,EVAP)
            call LOCAL_MIN(XTEMP,esys,NDIM,SEED,CFLAG,MIN_S,
     &MIN_T,EVAP)
C            write(*,*) J1,J2,",mutation=",esys,",min_pot=",min_pot
            if (esys.lt.min_pot)then
C              write(*,*) J1,J2,",mutation=",esys,",min_pot=",min_pot
              min_pot=esys
              do J3=1,ndim
                xnew(J3)=xtemp(J3)
              enddo
            endif
            xtemp(J2*3-2)=xinib(J2*3-2)
            xtemp(J2*3-1)=xinib(J2*3-1)
            xtemp(J2*3)=xinib(J2*3)
            xtemp(J1*3-2)=xinib(J1*3-2)
            xtemp(J1*3-1)=xinib(J1*3-1)
            xtemp(J1*3)=xinib(J1*3) 
          enddo
C          pause
        enddo
C        call POTENTIAL(ndim,xnew,grad,esys,RMS,.false.,EVAP)
C       write(*,*) "after mutation,min=",esys
C        pause
         return
C      else
        if(FIX_GFAC.eq.1) GFAC=TEMP
        do j1=1,ndim
          xnew(j1)=xinib(j1)
        enddo
        j3=0
        if(dim3.eq.0) then       !for 3N   
          do j1=1,ndim/3
            p=rnd(seed)
            j2=j1*3
            pass=.false.
C          if(sel) pause"after finish2"        
C          write(*,*) "step 1,j1=",j1,xnew(j2-2),xnew(j2-1),xnew(j2)
              j4=j3
              sel=.false.
          enddo
        endif
c      endif
c      pause
      return
      end
        
