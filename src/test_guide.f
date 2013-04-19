      program test_guide
      integer j1,j2,j3,j4,j5,gbin_total,guide(10000000),dim3,atom_num
      integer seed,gbin,gbin1,gbin2,gbin3,thredu,thredd,guide_total
      real xinib1,xinib2,xinib3,xinib(1000000),dg,radius,rnd,pot
      real dist(1000000),o_x,o_y,o_z
      character name1*6,ch*2
      logical gauss
      common /gauss1/gauss
      gauss=.false.
      o_x=0.D0
      o_y=0.D0
      o_z=0.D0
      write(*,*) "Please select your configuration type"
      write(*,*) "1.3N"
      write(*,*) "2.N"
      read(*,*) dim3
      write(*,*) "Please input number of atoms"
      read(*,*) atom_num
      radius=1.D0+(3.0D0*atom_num/17.77153175D0)**(1.0D0/3.0D0)
      radius=radius*2.0D0**(1.0D0/6.0D0)
      radius=radius*radius
      write(*,*) "Please input file name(plus guide.all):"
      read(*,*) name1
      open(10,file=name1//"guide.all",status="old")
      gbin_total=0
11    read(10,*,end=12)
      gbin_total=gbin_total+1
      goto 11
12    rewind(10)
      if(dim3.eq.1)then
        gbin=int(real(gbin_total)**(1./3.))
      endif
      if(dim3.eq.2)then
        gbin=gbin_total
      endif
      guide_total=0
      j2=0
      do j1=1,gbin_total
        read(10,*) j5,guide(j1)
        if(guide(j1).ne.0) j2=j2+1
        guide_total=guide_total+guide(j1)
      enddo
      write(*,*) "Please input the thredshold(default=",guide_total/j2,
     &")"
      write(*,*) "Lower thredshold:"
      read(*,*) thredd
      write(*,*) "Higher thredshold:(can be 0 for no upper bond)"
      read(*,*) thredu
      dg=(2.*sqrt(radius))/real(gbin)
      j2=0
      open(13,file="test_guide.xyz",status="replace")
      open(14,file="test_guide.num",status="replace")
      open(15,file="test_guide.dat",status="replace")
      open(18,file="test_gauss.dat",status="replace")
      do j1=1,gbin_total
        if(guide(j1).gt.thredd.and.guide(j1).le.thredu)then
          if (dim3.eq.1)then
            j2=j2+1
            j3=j2*3
            do j4=1,guide(j1)
              if(j4.ge.100) goto 15
              gbin3=int(j1/(gbin*gbin))
              gbin2=mod(j1,(gbin*gbin))
              gbin1=mod(gbin2,gbin)-1
              gbin2=int(gbin2/gbin)
              xinib(j3-2)=real(gbin1)*dg-sqrt(radius)
              xinib(j3-1)=real(gbin2)*dg-sqrt(radius)
              xinib(j3)=real(gbin3)*dg-sqrt(radius)
            xinib1=xinib(j3-2)+(dg*rnd(seed)*(-1)**int(rnd(seed)))
            xinib2=xinib(j3-1)+(dg*rnd(seed)*(-1)**int(rnd(seed)))
            xinib3=xinib(j3)+(dg*rnd(seed)*(-1)**int(rnd(seed)))
            write(14,*) xinib1,xinib2,xinib3
            enddo
15          write(*,*) j2,j1,guide(j1),xinib(j3-2),xinib(j3-1),xinib(j3)
     &,gbin1,gbin2,gbin3
            dist(j2)=sqrt((xinib(j3-2)-o_x)**2+(xinib(j3-1)-o_y)**2+
     &(xinib(j3)-o_z)**2)
            write(15,*) j1,guide(j1)
            write(18,*) dist(j2),guide(j1)
          endif
          if (dim3.eq.2)then
            j2=j2+1
            do j4=1,guide(j1)
              xinib(j2)=real((j1-1))*dg-sqrt(radius)
              xinib1=xinib(j2)+(dg*rnd(seed)*(-1)**int(rnd(seed)))
            enddo       
          endif  
        endif
      enddo
      write(13,*) j2
      write(13,*) 
      do j1=1,j2
        j3=j1*3
        write(13,*) "H",xinib(j3-2),xinib(j3-1),xinib(j3)
      enddo
      open(17,file="test_guide_g.xyz",status="replace")
      open(16,file=name1//".xyz",status="old")
      read(16,*) j1,pot
      read(16,*)
      do j1=1,atom_num
        j2=j1*3
        read(16,*) ch,xinib1,xinib2,xinib3 
        gbin1=int((xinib1+sqrt(radius))/dg)
        gbin2=int((xinib2+sqrt(radius))/dg)
        gbin3=int((xinib3+sqrt(radius))/dg)
        gbin_total=(1+gbin1+gbin2*gbin+gbin3*gbin*gbin)
        write(17,*) gbin_total,10000
      enddo
      write(*,*) "Radius=",radius,",DG=",dg
      close(10)
      close(13)    
      close(14)
      close(15)
      close(16)
      close(17)
      close(18)
      stop
      end
