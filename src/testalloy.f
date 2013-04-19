       Program testalloy
       implicit none
       integer i,A,SW,JMAX
       real*8 x(1000),enr,eshlon_a,eta_a,p_a,q_a,rzero_a
       real*8 eshlon_b,eta_b,p_b,q_b,rzero_b
       real*8 eshlon_ab,eta_ab,p_ab,q_ab,rzero_ab,VTR(1000)
       real*8 radius,none,GRAD(1000),VT(1000),VMIN,VMAX,DMAX
       character kinda*2,kindb*2,kindab*2
       LOGICAL EVAP,GRADT
       COMMON /parametera/ eshlon_a,eta_a,p_a,q_a,rzero_a
       COMMON /parameterb/ eshlon_b,eta_b,p_b,q_b,rzero_b
       COMMON /parameterab/ eshlon_ab,eta_ab,p_ab,q_ab,rzero_ab
       COMMON /accuracy3/ RADIUS,NONE
       COMMON /alloy1/ A
       COMMON /cgsw/ SW
       COMMON /GRAD1/ GRAD
       COMMON /GRAD2/ GRADT
       COMMON /vt1/ VT,VTR
       COMMON /vt2/ JMAX
       COMMON /vt3/ VMIN,VMAX,DMAX
       RADIUS=100.D0
       open(10,file="parameter1.dat",status="old")
       read(10,*) kinda,eshlon_a,eta_a,p_a,q_a,rzero_a
       close(10)
       open(11,file="alloy.dat",status="old")
       read(11,*) A
       close(11)
       open(12,file="parameter2.dat",status="old")
       read(12,*) kindb,eshlon_b,eta_b,p_b,q_b,rzero_b
       close(12)
       open(13,file="parameter3.dat",status="old")
       read(13,*) kindab,eshlon_ab,eta_ab,p_ab,q_ab,rzero_ab
       close(13)
       open(20,file="0038",status="old")
       do i=1,38*3,3
         read(20,*) x(i),x(i+1),x(i+2)
         write(*,*) x(i),x(i+1),x(i+2)
       enddo
       call alloy(38*3,x,enr,evap)
       write(*,*) "ENERGY=",enr,evap,A
       close(20)
       stop
       end
