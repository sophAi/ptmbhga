        program testfunc
        implicit none
        character pos*4,kind*2
        real*8 outx(2000),pot,func,dist(2000),tol
        integer atom_num,t1,t2,t3,t4,j1,j2,j3,eff_co,t,same(200,200)
        real*8 ecoh,zb,p,q,rzero,eta,eshlon,r,d,rms
        real*8 LL,ene,enef,enef2,simp,cg,dist2(200,200)
        integer j4,j5,same2(200,200,200)
        logical EVAP
        common/accuracy1/ r,d
        common/accuracy3/ enef,enef2
        common/parameter1/ eshlon,eta,p,q,rzero
        write(*,*) "test energy for Many-body SMA potential"
        write(*,*) " ---------------------------------------"
        open(69,file="parameter1.dat",status="old")
        read(69,*) kind,eshlon,eta,p,q,rzero
        close(69) 
        write(*,*) " kind eshlon  eta     p     q   rzero"
        write(*,*) kind,eshlon,eta,p,q,rzero      
        open(44,file="accuracy.dat",status="old")
        read(44,*) LL,ene,enef,enef2,simp,cg
        close(44)                              
        write(*,*) "Please input the source file:"
        read(*,*) atom_num
        write(*,*) "Please input the tolerence:"
        read(*,*) tol
        eff_co=atom_num*3
        t=0
        t1=atom_num/1000
        t2=atom_num/100
        t2=mod(t2,10)
        t3=atom_num/10
        t3=mod(t3,10)
        t4=mod(atom_num,10)
        pos=char(t1+48)//char(t2+48)//char(t3+48)//char(t4+48)
        open(1,file=pos,status="old") 
        j2=0
        write(*,*) "-------------------------------------------"
        do j1=1,atom_num*3,3
          j2=j2+1
          read(1,*) outx(j1),outx(j1+1),outx(j1+2)
          write(*,*) kind,j2,"'atom",(outx(j1+j3),j3=0,2)
        enddo
        call centre(outx,atom_num)
        do j1=1,atom_num
          j2=j1*3
          dist(j1)=outx(j2)**2+outx(j2-1)**2+outx(j2-2)**2
          write(*,*) j1,dist(j1)
        enddo
        do j1=1,atom_num
          j4=j1*3
          j3=0
          do j2=j1+1,atom_num
            j5=j2*3
            dist2(j1,j2)=(outx(j4)-outx(j5))**2+(outx(j4-1)-
     &outx(j5-1))**2+(outx(j4-2)-outx(j5-2))**2
            dist2(j1,j2)=dsqrt(dist2(j1,j2))
            dist2(j2,j1)=dist2(j1,j2)
            if (abs(dist(j1)-dist(j2)).le.tol.and.j1.ne.j2)then
C              write(*,*) j1,j2,(dist(j1)+dist(j2))/2
              j3=j3+1
              same(j1,j3)=j2
            endif
          enddo
        write(*,*) j1," the same",(same(j1,j2),j2=1,j3)         
        enddo
        pause
        do j1=1,atom_num
          do j2=1,atom_num
            if (j1.ne.j2)then
              j4=0
              do j3=j2+1,atom_num
                if (abs(dist2(j1,j2)-dist2(j1,j3)).lt.tol)then
                  j4=j4+1
                  same2(j1,j2,j4)=j3
                endif
              enddo
              write(*,*) j1,j2,",the same2",(same2(j1,j2,j5),j5=1,j4)
            endif
          enddo
          pause
        enddo
        pot=func(outx,atom_num*3,t) 
        write(*,*) "-------------------------------------------"
        write(*,*) " Test total potential=",pot," (ev)"
        write(*,*) " The Eb(per atom)=",(pot/atom_num)," (ev)"
        stop
        end
            
