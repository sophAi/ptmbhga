       program gauss_test
       integer seed
       real*8 rnd,value
       logical gauss
       common/gauss1/ gauss
       gauss=.true.
       seed=1000
       value=rnd(seed)
       write(*,*) value

       stop
       end
       
