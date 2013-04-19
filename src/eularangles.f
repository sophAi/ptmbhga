       subroutine eularangles_gp(xnew,ndim,seed)
       implicit none
       integer ndim,seed,i,j,k,num,p1,p2,p3,p4,DG_RECOVER,dim3,pdb_rec,
     &ndim_temp,dim3_temp
       real*8 xold(ndim),xnew(ndim),rnd,pei,r
       real*8 angle1,angle2,angle3,fac
       character name2*4,name6*6
       parameter(pei=3.14159D0)
       common/accuracy1/ r
       common/seeding/ num
       common/seeding2/ p1,p2,p3,p4
       common/name1/ name6
       common/pdbrec/ pdb_rec,dim3
       angle1=pei*rnd(seed)
       angle2=pei*rnd(seed)
       angle3=pei*rnd(seed)
       if(num.eq.0) then
         do i=1,ndim
           xnew(i)=r*rnd(seed)
         enddo
         goto 20
       endif
       name2=char(p1+48)//char(p2+48)//char(p3+48)//char(p4+48)
       open(61,file=name2,status="old")
       do i=1,(num*3),3     
         read(61,*) (xold(i+j),j=0,2)
         xnew(i)=((dcos(angle3)*dcos(angle1)-dcos(angle2)*dsin(angle1)
     &*dsin(angle3))*xold(i)
     &+(dcos(angle3)*dsin(angle1)+dcos(angle2)*dcos(angle1)*dsin
     &(angle3))*xold(i+1)
     &+(dsin(angle3)*dsin(angle2))*xold(i+2))
  
         xnew(i+1)=(((-1.)*dsin(angle3)*dcos(angle1)-dcos(angle2)*
     &dsin(angle1)*dcos(angle3))*xold(i)
     &+((-1.)*dsin(angle3)*dsin(angle1)+dcos(angle2)*dcos(angle1)
     &*dcos(angle3))*xold(i+1)
     &+(dcos(angle3)*dsin(angle2))*xold(i+2))
 
         xnew(i+2)=((dsin(angle2)*dsin(angle1))*xold(i)
     &+((-1.)*dsin(angle2)*dcos(angle1))*xold(i+1)
     &+dcos(angle2)*xold(i+2))
       enddo
       close(61)
       do k=1,(ndim-(num*3)),3
         xnew(ndim-k+1)=r*rnd(seed)
         xnew(ndim-k)=r*rnd(seed)
         xnew(ndim-k-1)=r*rnd(seed)
       enddo
20     return
       end
 
       subroutine eularangles_alloy(xnew,ndim,seed)
       implicit none
       integer ndim,seed,i,j,k,num,p1,p2,p3,p4,DG_RECOVER,dim3,pdb_rec
     &,ndim_temp,dim3_temp
       real*8 xold(ndim),xnew(ndim),rnd,pei,r
       real*8 angle1,angle2,angle3,fac
       character name2*4,name6*6
       parameter(pei=3.14159D0)
       common/accuracy1/ r
       common/seeding/ num
       common/seeding2/ p1,p2,p3,p4
       common/name1/ name6
       common/pdbrec/ pdb_rec,dim3
       angle1=pei*rnd(seed)
       angle2=pei*rnd(seed)
       angle3=pei*rnd(seed)
       if(num.eq.0) then
         do i=1,ndim
           xnew(i)=r*rnd(seed)
         enddo
         goto 20
       endif
       name2=char(p1+48)//char(p2+48)//char(p3+48)//char(p4+48)
       open(61,file=name2,status="old")
       do i=1,(num*3),3
         read(61,*) (xold(i+j),j=0,2)
         xnew(i)=((dcos(angle3)*dcos(angle1)-dcos(angle2)*dsin(angle1)
     &*dsin(angle3))*xold(i)
     &+(dcos(angle3)*dsin(angle1)+dcos(angle2)*dcos(angle1)*dsin
     &(angle3))*xold(i+1)
     &+(dsin(angle3)*dsin(angle2))*xold(i+2))

         xnew(i+1)=(((-1.)*dsin(angle3)*dcos(angle1)-dcos(angle2)*
     &dsin(angle1)*dcos(angle3))*xold(i)
     &+((-1.)*dsin(angle3)*dsin(angle1)+dcos(angle2)*dcos(angle1)
     &*dcos(angle3))*xold(i+1)
     &+(dcos(angle3)*dsin(angle2))*xold(i+2))

         xnew(i+2)=((dsin(angle2)*dsin(angle1))*xold(i)
     &+((-1.)*dsin(angle2)*dcos(angle1))*xold(i+1)
     &+dcos(angle2)*xold(i+2))
       enddo
       close(61)
       do k=1,(ndim-(num*3)),3
         xnew(ndim-k+1)=r*rnd(seed)
         xnew(ndim-k)=r*rnd(seed)
         xnew(ndim-k-1)=r*rnd(seed)
       enddo
20     return
       end 

       subroutine dih_chima(xnew,ndim,seed)
       integer ndim,seed,num,p1,p2,p3,p4,DG_RECOVER,dim3,pdb_rec
     &,ndim_temp,dim3_temp,J1,J2,J3,J4,gmin(ndim),
     &gmin_temp,J1_min,J2_min
       real*8 xnew(ndim),r,RADIUS,NONE,DG,pei,rnd
       parameter(pei=3.1415962D0)
       character name6*6
       common /accuracy1/ r
       common /accuracy3/ RADIUS,NONE
       common /seeding/ num
       common /seeding2/ p1,p2,p3,p4
       common /pdbrec/ pdb_rec,dim3
       if(DG_RECOVER.eq.1)then
         open(49,file=name6//"_guide.all",status="old")
         read(49,*) ndim_temp,dim3_temp
         if(ndim_temp.ne.ndim.or.dim3_temp.ne.dim3)
     &then
           write(*,*) "Waring!!READING DATA ERROR"
           goto 22
         endif
         do J1=1,ndim
           gmin_temp=0
         enddo
         close(49)
         do J1=1,ndim
           xnew(J1)=dble(gmin(J1)-1)*DG-dsqrt(RADIUS)
         enddo  
       endif
22     do J1=1,ndim
         xnew(J1)=0.0001D0*rnd(seed)
       enddo
       return
       end 
    
       subroutine dim_chima(xnew,ndim,seed)
       integer ndim,seed,num,p1,p2,p3,p4,DG_RECOVER,dim3,pdb_rec
     &,ndim_temp,dim3_temp,J1,J2,J3,J4,gmin(ndim),
     &gmin_temp,J1_min,J2_min,eff_co
       real*8 xnew(ndim),r,RADIUS,NONE,DG,dih(ndim/3-2),x(6000),rnd
       common /accuracy1/ r
       common /accuracy3/ RADIUS,NONE
       common /seeding/ num
       common /seeding2/ p1,p2,p3,p4
       common /pdbrec/ pdb_rec,dim3
       eff_co=ndim/3-2
       do J1=1,ndim/3-2
         dih(J1)=r*rnd(seed)
       enddo
       call dih_gen(eff_co,dih,x)
       do J1=1,ndim
         xnew(J1)=x(J1)
       enddo
       return
       end
        
