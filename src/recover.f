       program recover
       implicit none
       integer j0,j1,j2,j3,t1,t2,t3,t4,ndim,ens_num,temp1
       real*8 pot(1000),xini(1000,3000),temp2,temp(3000)
       character name*6,kinda*2
20     format(F25.13,F15.10)
       write(*,*) "Please input the particle number"
       read(*,*) ndim
       write(*,*) "Please input the meterial"
       read(*,*) kinda
       ndim=ndim*3
       t1=(ndim)/3000
       t2=(ndim)/300
       t2=mod(t2,10)
       t3=(ndim)/30
       t3=mod(t3,10)
       t4=mod(((ndim)/3),10)
       name=char(t1+48)//char(t2+48)//char(t3+48)//char(t4+48)//kinda
       open(13,file=name//"_recover.rec",status="old")
       read(13,*) temp2
       ens_num=0
21     read(13,*,end=22)
         ens_num=ens_num+1
       goto 21
22     rewind(13)
       close(13)
       ens_num=ens_num/ndim
       write(*,*) "Reading "
     &,name//"_recover.rec file,the ensamble number is ",ens_num
       open(12,file=name//".xyz",status="old")
       read(12,*) temp1,temp2
       read(12,*) 
       do j1=1,ndim/3
         j2=j1*3
         read(12,*) kinda,temp(j2-2),temp(j2-1),temp(j2)
       enddo
       write(*,*)
       write(*,*) "The ",name//".xyz file has been loaded successfully"
       write(*,*) 
       write(*,*)"What ensamble do you want to replace?(1 ~",ens_num,")"
       read(*,*) j3
       open(10,file=name//"_recover.rec",status="old")
       open(11,file=name//"_recover.new",status="replace")
       read(10,*) j0
       write(11,*) j0
       do j1=1,ens_num
         do j2=1,ndim
           if(j1.eq.j3)then
             write(11,20) temp2,temp(j2)
           else
             read(10,20) pot(j1),xini(j1,j2)
             write(11,20) pot(j1),xini(j1,j2)
           endif
         enddo
       enddo
       write(*,*)
       write(*,*) "You data have been wrote in ",
     &name//"_recover.new file"
       write(*,*)
       write(*,*) "Please rename this file to ",
     &name//"_recover.rec "
       close(10)
       close(11)
       close(12)
       STOP
       END

