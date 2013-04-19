	real function rnd(seed) 
        implicit none
        integer seed
        real*8 ran,w,y1,y2,ran1,ran2
        logical gauss
        common/gauss1/ gauss
        if (.not.gauss) then
          rnd=ran(seed)
        else
          w=2.D0
          do while(w.ge.1.D0)
            ran1=ran(seed+1)
            ran1=ran1*2.D0-1.D0
            ran2=ran(seed-1)
            ran2=ran2*2.D0-1.D0
            w=ran1**2+ran2**2
          enddo   
          w=dsqrt((-2.D0*dlog(w))/w)
          y1=ran1*w
          y2=ran2*w
          pick=ran(seed+3)
          if(pick.ge.0.5D0) then
             rnd=y1
          else
             rnd=y2
          endif   
        endif
        return
        end
