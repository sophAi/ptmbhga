       Subroutine simplex_min(xinib,min_pot,eff_co,seed,cflag,EVAP)
       implicit none
       integer eff_co,seed,M
       real*8 xinib(eff_co),min_pot,RMS
       real*8 x(eff_co),y(eff_co+1),p(eff_co+1,eff_co)
       real*8 DELTAX,FTOL,EREAL,GRAD(eff_co)
       integer i,j,mini
       common /accuracy5/ FTOL
       logical EVAP,cflag
       mini=1
       EVAP=.false.
       cflag=.false.
*** The parameter of simplex method ***
C       DELTAX=4.D0
       DELTAX=0.4D0
       M=eff_co+1
***************************************
       do i=1,M
         do j=1,eff_co
           if(i.eq.1) then
             p(i,j)=xinib(j)
             x(j)=p(i,j) 
           else if(j.eq.(i-1)) then
             p(i,j)=p(1,j)+DELTAX
             x(j)=p(i,j)
           else
             p(i,j)=p(1,j)
             x(j)=p(i,j)
           endif
         enddo 
         call POTENTIAL(eff_co,x,GRAD,EREAL,RMS,.false.,EVAP)
         y(i)=EREAL
       enddo 
       call amoeba(p,y,FTOL,eff_co,cflag)
       do i=2,M
         if(y(i).lt.y(mini)) then
           mini=i
         endif
       enddo
       do i=1,eff_co
         xinib(i)=p(mini,i)
       enddo
       min_pot=y(mini)
       end
       
       

