      subroutine startend(myid,nproc,is1,is2,istart,iend)
      implicit none
      integer myid,nproc,is1,is2,istart,iend,ilength
      integer iblock,ir
      ilength=is2-is1+1
      iblock=int(ilength/nproc)
      ir=ilength-iblock*nproc
      if(myid.lt.ir) then
         istart=is1+myid*(iblock+1)
         iend=istart+iblock
      else
         istart=is1+myid*iblock+ir
         iend=istart+iblock-1
      endif
      if(ilength.lt.1) then
         istart=1
         iend=0
      endif
      return
      end

