C     The GA operators are good random move method!!
      subroutine step_move(xinib_org1,xinib_org2,xinib_new1,xinib_new2,
     &ndim,seed)
      implicit none
      integer kcount,ndim,seed,sw_ga,sw_eular2
      real*8 xinib_org1(ndim),xinib_org2(ndim),xinib_new1(ndim),p,
     &xinib_new2(ndim),eular_fac2,rnd
      real*8 mv1,mv2,mv3,mv4,mv5,mv6,mv7,mv8,mv9,mv10
      common/move4/ mv1,mv2,mv3,mv4,mv5,mv6,mv7,mv8,mv9,mv10
      common/SW_GA2/ sw_ga
      common/eularangle2/sw_eular2
      common/eularangle4/eular_fac2
      p=rnd(seed)
C      write(*,*) mv1,mv2,mv3,mv4,mv5,mv6,mv7,mv8,mv9,mv10,p,ndim
C      pause
c      do kcount=1,ndim
C        write(*,*) kcount,"xinib_org1=",xinib_org1(kcount)
c        write(*,*) kcount,"xinib_org2=",xinib_org2(kcount)
C      enddo
C      pause
C      call eularangles_tran(xinib_org,ndim,seed)
      if (p.lt.mv1)then
C     3N Particles moving
        sw_ga=0
        if (sw_eular2.eq.1)then
          call eularangles_tran(xinib_org1,ndim,seed,eular_fac2)
        endif
        call tran(xinib_org1,xinib_new1,ndim,seed)
      else if (p.lt.(mv1+mv2)) then
C     3N Particles Mutation
        sw_ga=0
        call mutation(xinib_org1,xinib_new1,ndim,seed)   
      else if (p.lt.(mv1+mv2+mv3)) then
C     N dimension movining
        sw_ga=0
        if (sw_eular2.eq.1)then
          call eularangles_tran(xinib_org1,ndim,seed,eular_fac2)
        endif  
        call tran_n(xinib_org1,xinib_new1,ndim,seed,.false.)
      else if (p.lt.(mv1+mv2+mv3+mv4)) then
C     N dimension mutation
        if (sw_eular2.eq.1)then
          call eularangles_tran(xinib_org1,ndim,seed,eular_fac2)
        endif
        sw_ga=0
        call inversion(xinib_org1,xinib_new1,ndim,seed)
      else if (p.lt.(mv1+mv2+mv3+mv4+mv5)) then
        sw_ga=0
        call dih_tran(xinib_org1,xinib_new1,ndim,seed)
      else if (p.lt.(mv1+mv2+mv3+mv4+mv5+mv6)) then
        if (sw_eular2.eq.1)then
          call eularangles_tran(xinib_org1,ndim,seed,eular_fac2)
        endif
        sw_ga=0
       call tran_n(xinib_org1,xinib_new1,ndim,seed,.true.)
      else if (p.lt.(mv1+mv2+mv3+mv4+mv5+mv6+mv7)) then
        if (sw_eular2.eq.1)then
          call eularangles_tran(xinib_org1,ndim,seed,eular_fac2)
          call eularangles_tran(xinib_org2,ndim,seed,eular_fac2)
        endif
        sw_ga=0
        call arithmetic(xinib_org1,xinib_org2,xinib_new1,ndim) 
      else if (p.lt.(mv1+mv2+mv3+mv4+mv5+mv6+mv7+mv8)) then
        if (sw_eular2.eq.1)then
          call eularangles_tran(xinib_org1,ndim,seed,eular_fac2)
          call eularangles_tran(xinib_org2,ndim,seed,eular_fac2)
        endif
        sw_ga=0
        call geometic(xinib_org1,xinib_org2,xinib_new1,ndim)
      else if (p.lt.(mv1+mv2+mv3+mv4+mv5+mv6+mv7+mv8+mv9)) then
        if (sw_eular2.eq.1)then
          call eularangles_tran(xinib_org1,ndim,seed,eular_fac2)
          call eularangles_tran(xinib_org2,ndim,seed,eular_fac2)
        endif
        sw_ga=1
        call crossing(xinib_org1,xinib_org2,xinib_new1,xinib_new2,ndim,
     &seed)
      else
        if (sw_eular2.eq.1)then
          call eularangles_tran(xinib_org1,ndim,seed,eular_fac2)
          call eularangles_tran(xinib_org2,ndim,seed,eular_fac2)
        endif
        sw_ga=1
        call twopoint(xinib_org1,xinib_org2,xinib_new1,xinib_new2,ndim,
     &seed)
      endif
      return
      end

