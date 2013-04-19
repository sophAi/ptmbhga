       function func(x,ndim,t)
       implicit none
       integer ndim
       real*8 x(ndim),func,attpot,attpot2,repot,repot2
       real*8 temp(ndim),at_dis,pote,att_pot,rep_pot,L,d
       integer i,j,k,t
       real*8 ecoh,zb,p,q,rzero,eta,eshlon,r,lo,lo2
       common/accuracy1/ r,d
       common/accuracy3/lo,lo2
       common/parameter1/ eshlon,eta,p,q,rzero
***if d eq 1 ,use paper r=...***
       if(d.eq.(1.))then
         L=((ndim*2.)**(1./3.))*r
         lo2=r
       else
         L=lo
       endif
*----- These are the parameters of the Metal Cluster Energy ---       
       pote=0.D0
       if(ndim.gt.3) then
         do i=1,ndim
           temp(i)=x(i)
         enddo
       endif
       attpot2=0.D0
       repot2=0.D0
       do i=1,ndim,3
         attpot=0.D0
         repot=0.D0
         do j=1,ndim,3
           t=0
           at_dis=0.D0
           if (j.ne.i) then
             do k=0,2,1
               at_dis=at_dis+(temp(i+k)-temp(j+k))**2
             enddo
             at_dis=dsqrt(at_dis)
             if(at_dis.ge.L)then
               t=1
               return
             endif            
             attpot=attpot+att_pot(at_dis,rzero,q)
             repot=repot+rep_pot(at_dis,rzero,p)    
           endif
         enddo
         attpot2=attpot2+dsqrt(attpot)
         repot2=repot2+repot
       enddo
       pote=(eshlon*repot2)-(eta*attpot2)
       func=pote
       if(abs(func).ge.9999999) then
         t=1
         return
       endif
       end
       function att_pot(x,rzero,q)
       implicit none
       real*8 x,att_pot,attrct
       real*8 rzero,q,coff
*************??************
       if(x.lt.1.d-15) then
         att_pot=99999.D+15
       else
         coff=(1.-x)
         attrct=dexp(2.*q*coff)
         att_pot=attrct
       endif
       end
       function rep_pot(x,rzero,p)
       implicit none
       real*8 x,rep_pot,repuls
       real*8 rzero,p,coff
***********??*************************
       if (x.lt.1.d-15) then
         rep_pot=99999.D+15
       else
         coff=(1.-x)
         repuls=dexp(p*coff)
         rep_pot=repuls
       endif
       end
