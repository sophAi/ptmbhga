       program testdist
       integer i,j,atom_num1,a1,flag1,cycle_num,np(1000),type1
       real*8 k1a,k1b,k1c,x1(3000),y1(3000),z1(3000),dist1a,dist1b
       real*8 energy1,R01a,ratio1a,ratio1b,avdist1a,avdist1b,avdist1c
       real*8 m1a(3000),m1b(3000),m1c(3000),R01c,ratio1c,cent_dist(1000)
       real*8 distm1a(3000),centx,ctnty,centz,dist1c,R01b,add_fac
       real*8 ratiocent,avdistcent(1000),maxtemp,kx(1000),total_bond
       character name_xyz1*12,atom1a*2,atom1b*2,quick_name*2
       character yn*1,name_dist1*8
       character atom_name1a*2,atom_name1b*2
       k1=0.D0
       avdist1=0.D0
       centx=0.D0
       centy=0.D0
       centz=0.D0
       write(*,*) "Quick input?(y/n)"
       read(*,*) yn
       if(yn.eq."y")then
         add_fac=0.01D0
         cycle_num=50
         R01a=2.556D0
         ratio1a=1.D0
         R01b=2.884D0
         ratio1b=1.D0
         R01c=2.556D0
         ratio1c=1.D0
         ratiocent=1.2D0
         write(11,*) cent_dist(i),avdistcent(i),atom_num1
         yn="n"
       write(*,*) "Please input the number of cluster #1(ex:01,02,03..)"
         read(*,*) quick_name
         name_xyz1="0038CA"//quick_name//".xyz"
         goto 30
       endif
       write(*,*) "Please input the add factor(default=0.01D0)"
       read(*,*) add_fac
       write(*,*) "please input the cycle number(default=50)"
       read(*,*) cycle_num
       write(*,*) "Please input the file name of cluster #1:(*.xyz)"
       read(*,*) name_xyz1
       write(*,*) "please input R0 of atom a-a(default=2.556D0)"
       read(*,*) R01a
       write(*,*) "Please input the ratio of atom a-a(default=1.D0)"
       read(*,*) ratio1a
       write(*,*) "please input R0 of atom b-b(default=2.884D0)"
       read(*,*) R01b
       write(*,*) "Please input the ratio of atom b-b(default=1.D0)"
       read(*,*) ratio1b
       write(*,*) "please input R0 of atom a-b(default=2.556D0)"
       read(*,*) R01c
       write(*,*) "Please input the ratio of atom a-b(default=1.D0)"
       read(*,*) ratio1c
       write(*,*) 
     &"Please input the ratio of internal strain(default=1.2D0)"
       read(*,*) ratiocent
30     name_dist1=name_xyz1
       open(10,file=name_dist1//".dist",status="replace")
       open(1,file=name_xyz1,status="old")
       read(1,*) atom_num1,energy1
       a1=1
       flag1=0
       do i=1,atom_num1 
          read(1,*) atom1a,x1(i),y1(i),z1(i)
          centx=centx+x1(i)
          centy=centy+y1(i)
          centz=centz+z1(i)
          if(atom1a.eq.atom1b.and.flag1.eq.0)then 
            a1=a1+1
            atom_name1a=atom1a
          else
            atom_name1b=atom1a
            if(i.ne.1) flag1=1
          endif
          atom1b=atom1a
       enddo
       centx=centx/atom_num1
       centy=centy/atom_num1
       centz=centz/atom_num1
       close(1)
       write(*,*)"Read cluster #1 A=",a1,",atom1=",atom_name1a,",R0=",
     &R01a,",atom2=",atom_name1b,",R0=",R01b
      write(10,*)"Read cluster #1 A=",a1,",atom1=",atom_name1a,",R0=",
     &R01a,",atom2=",atom_name1b,",R0=",R01b
       do n=1,cycle_num
         avdist1a=0.D0
         avdist1b=0.D0
         avdist1c=0.D0
         k1a=0.D0
         k1b=0.D0
         k1c=0.D0
         do i=1,atom_num1
           m1a(i)=0.D0
           m1b(i)=0.D0
           m1c(i)=0.D0
           distm1a(i)=0.D0
           do j=1,atom_num1
             if(i.ne.j)then
               dist1a=(x1(i)-x1(j))**2+(y1(i)-y1(j))**2+(z1(i)-z1(j))**2
               dist1a=dsqrt(dist1a)
               if((i.le.a1).and.(j.le.a1))then 
                 if(dist1a.le.(ratio1a*R01a))then
                   avdist1a=avdist1a+dist1a
                   distm1a(i)=distm1a(i)+dist1a
                   k1a=k1a+1.D0
                   m1a(i)=m1a(i)+1.D0
                 endif
               else if((i.gt.a1).and.(j.gt.a1))then
                 if(dist1a.le.(ratio1b*R01b))then
                   avdist1a=avdist1a+dist1a
                   distm1a(i)=distm1a(i)+dist1a
                   k1a=k1a+1.D0
                   m1a(i)=m1a(i)+1.D0
                 endif
               else
                 if(dist1a.le.(ratio1c*R01c))then
                   avdist1a=avdist1a+dist1a
                   distm1a(i)=distm1a(i)+dist1a
                   k1a=k1a+1.D0
                   m1a(i)=m1a(i)+1.D0
                 endif
               endif
             endif
           enddo
           distm1a(i)=distm1a(i)/m1a(i)
         enddo
         do i=1,a1
           do j=i+1,a1
             dist1b=(x1(i)-x1(j))**2+(y1(i)-y1(j))**2+(z1(i)-z1(j))**2
             dist1b=dsqrt(dist1b)
             if(dist1b.le.(ratio1a*R01a))then
               avdist1b=avdist1b+dist1b
               k1b=k1b+1.D0
               m1b(i)=m1b(i)+1.D0
             endif
           enddo
         enddo
         do i=a1+1,atom_num1
           do j=i+1,atom_num1
             dist1c=(x1(i)-x1(j))**2+(y1(i)-y1(j))**2+(z1(i)-z1(j))**2
             dist1c=dsqrt(dist1c)
             if(dist1c.le.(ratio1b*R01b))then
               avdist1c=avdist1c+dist1c
               k1c=k1c+1.D0
               m1c(i)=m1c(i)+1.D0
             endif
           enddo
         enddo
         avdist1a=avdist1a/k1a
         avdist1b=avdist1b/k1b
         avdist1c=avdist1c/k1c
         write(10,*) "<STEP>=",n
         write(10,*) "For ",name_xyz1
         write(10,*) "Energy is ",energy1
         write(10,*) "Average distance is ",avdist1a
         write(10,*) "Center of mass=",centx,centy,centz
         open(11,file=name_xyz1,status="old")
         read(11,*) atom_num1,energy1
         write(10,*) 
     &"# name  cent_dist  #nearest_bond  #aa  #bb  avg_bondist" 
         do i=1,atom_num1
            read(11,*) atom1a,x1(i),y1(i),z1(i)
         cent_dist(i)=(x1(i)-centx)**2+(y1(i)-centy)**2+(z1(i)-centz)**2
            cent_dist(i)=dsqrt(cent_dist(i))
         write(10,*) i," ",atom1a,cent_dist(i),m1a(i),m1b(i),m1c(i)
     &,distm1a(i)
         enddo
         close(11)
        write(10,*) atom_name1a,"=",avdist1b,",HIT=",k1b,
     &",",atom_name1b,"=",avdist1c,",HIT=",k1c
         write(10,*) atom_name1a,"/",atom_name1b," or ",atom_name1b,"/",
     &atom_name1a
         write(10,*) avdist1b/avdist1c," or ",avdist1c/avdist1b
         write(10,*) "============================================="
         write(10,*)
C         pause
         ratio1a=ratio1a+add_fac
         ratio1b=ratio1b+add_fac
         write(10,*) "============================================="
         write(10,*) "Ratio of cluster 1 ",atom_name1a,"=",ratio1a,
     &",dist=",ratio1a*R01a,",R0=",R01a
         write(10,*) "Ratio of cluster 1 ",atom_name1b,"=",ratio1b,
     &",dist=",ratio1b*R01b,",R0=",R01b
       enddo
       close(10)
       open(11,file=name_dist1//"a.txt",status="replace")
       open(12,file=name_dist1//"b.txt",status="replace")
       open(13,file=name_dist1//".txt",status="replace")
       do i=1,atom_num1
         kx(i)=0.D0
         do j=1,atom_num1
           if(i.ne.j)then
             dist1a=(x1(i)-x1(j))**2+(y1(i)-y1(j))**2+(z1(i)-z1(j))**2
             dist1a=dsqrt(dist1a)
             if((i.le.a1).and.(j.le.a1))then
               if(dist1a.le.(ratiocent*R01a))then
                 avdistcent(i)=avdistcent(i)+dist1a
                 kx(i)=kx(i)+1.D0
               endif
             else if((i.gt.a1).and.(j.gt.a1))then
               if(dist1a.le.(ratiocent*R01b))then
                 avdistcent(i)=avdistcent(i)+dist1a
                 kx(i)=kx(i)+1.D0
               endif
             else
               if(dist1a.le.(ratiocent*R01c))then
                 avdistcent(i)=avdistcent(i)+dist1a
                 kx(i)=kx(i)+1.D0
               endif
             endif
           endif
         enddo
         avdistcent(i)=avdistcent(i)/kx(i)
       enddo
       do i=1,atom_num1
         np(i)=i
       enddo
       do i=1,atom_num1
         do j=i+1,atom_num1
           if(cent_dist(np(i)).gt.cent_dist(np(j)))then
             maxtemp=np(i)
             np(i)=np(j)
             np(j)=maxtemp
           endif
         enddo
       enddo
       do i=1,atom_num1
         total_bond=avdistcent(np(i))*kx(np(i))
         if(np(i).le.a1)then
           type1=1
           write(11,*)
     &cent_dist(np(i)),avdistcent(np(i)),kx(np(i)),total_bond,
     &atom_num1-a1,type1 
         else
           type1=2
           write(12,*)
     &cent_dist(np(i)),avdistcent(np(i)),kx(np(i)),total_bond,
     &atom_num1-a1,type1
         endif
         write(13,*)
     &cent_dist(np(i)),avdistcent(np(i)),kx(np(i)),total_bond,
     &atom_num1-a1,type1
       enddo
       close(11)
       close(12)
       close(13)
       write(*,*) 
     &"The Output file is ",name_dist1,".dist , ",name_dist1,".txt"
       write(*,*) name_dist1,"a.txt and ",name_dist1,"b.txt"
       write(*,*) "The *.ist file can be use to study the internal" 
       write(*,*) "strain by xmgr,x is the distance from center of mass"
       write(*,*) "y is the average bond distance per atom"
       write(*,*) "COMPLETED!!"
       stop
       end
