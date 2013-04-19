        program angle_cal
        integer load
        real*8 degree,rad,pei
        parameter(pei=3.1415926D0)
        write(*,*) "P.J.Hsu's Degree<<-->>Rad calculator" 
14      write(*,*) "1.degree-->>rad"
        write(*,*) "2.rad-->>degree"
        read(*,*) load
        if (load.eq.1) then
          write(*,*) "Input angle in <degree>="
          read(*,*) degree
          rad=degree*pei/180.D0
          write(*,*) "=",rad,"<rad>"
        else if (load.eq.2)then
          write(*,*) "Input angle in <rad>="
          read(*,*) rad
          degree=rad*180.D0/pei
          write(*,*) "=",degree,"<degree>"
        else
          goto 14
        endif
        write(*,*) degree,"<degree>  =",rad,"<rad>"
        stop
        end
               
          
