C                       <<<   readDATA.f   >>>                 2002/3/26
      
        subroutine ReaRan(L_n)
        implicit double precision (A-H,O-Z)
        parameter (icutR2=100000)
        common /Rea1/arange,frange,drange,icutR
        write(*,*) "Doing moment"
        open(UNIT=8,FILE='En.dat')
        read(8,*) L_n
        close(8)
        write(*,*) 'Moment numbers =',L_n
        icutR=icutR2
c        pause ' You have "range" ?'
        open(11,file= 'range')
        read(11,*) arange
        frange = 0.d0
        read(11,*,end=10) frange
10      write(*,*) ' frange =',frange
        write(*,*) ' arange =',arange
        close(11)
        acut=dble(icutR)
        drange=(arange-frange)/acut
        write(*,*)
        write(*,*) ' Your range is from',frange,'to',arange
        write(*,*) '      and "drange" is:',drange,' !'
        write(*,*)
        return
        end

        subroutine ReaEn(ndim)
        implicit double precision(A-H,O-Z)
        double precision moment(20)
        character name*6
c        common /allf/L_n
        common /Rea2/moment
        common /name1/ name
        name="moment"
c        pause ' You have "En.dat" ?'
        write(*,*)
        open(UNIT=8,FILE='En.dat')
        read(8,*) L_n
c        write(*,*) 'Moment numbers =',L_n
        do i = 1, L_n
        read(8,*) moment(i)
        write(*,*) 'moment(',i,')=',moment(i)
        enddo
        close(8)
        write(*,*)
        return
        end
