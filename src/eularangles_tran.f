       subroutine eularangles_tran(xnew,ndim,seed,eular_fac)
       implicit none
       integer seed,ndim,i
       real*8 xnew(ndim),xold(ndim),pei,angle1,angle2,angle3,rnd,
     &eular_fac
       parameter(pei=3.14159D0)
C       angle1=angle1+rnd(seed)
C       angle2=angle2+rnd(seed)
C       angle3=angle3+rnd(seed)
       angle1=pei*eular_fac*rnd(seed)
       angle2=pei*eular_fac*rnd(seed)
       angle3=pei*eular_fac*rnd(seed)
C       angle3=pei*rnd(seed)
       do i=1,ndim
         xold(i)=xnew(i)
       enddo
       do i=1,ndim,3
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
       return
       end
