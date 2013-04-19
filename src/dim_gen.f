       subroutine dim_gen(xnew,ndim,seed)
       implicit none
       integer ndim,seed,j1,type
       real*8 xnew(ndim),vt(2000),vtr(2000)
       common/engtype/ type
       common/vt1/ vt,vtr
       do j1=1,ndim
         vtr(j1)=0
         vt(j1)=0
       enddo
       if(type.eq.1)then
         call eularangles_gp(xnew,ndim,seed)
       endif
       if(type.eq.2)then
         call eularangles_alloy(xnew,ndim,seed)
       endif
       if(type.eq.3)then
       
       endif
       if(type.eq.4)then
   
       endif
       if(type.eq.5)then
         call moment_gen(xnew,ndim,seed)
       endif
       if(type.eq.6)then
         
       endif
       if(type.eq.7)then
         call dih_chima(xnew,ndim,seed)
       endif
       if(type.eq.8)then
         call dim_chima(xnew,ndim,seed)
       endif
       if(type.eq.9)then
         call eularangles_alloy(xnew,ndim,seed)
       endif  
       return
       end
