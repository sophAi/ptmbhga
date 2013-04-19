*========================================================================
* File Name : rnd.f 
* Copyright (C) 2008-2011 Po-Jen Hsu <xanadu8850@pchome.com.tw>
* Creation Date : 19-04-2010
* Last Modified : Thu 16 Sep 2010 02:13:50 PM CST
* License : GPL (see bottom)
* Encoding : utf-8
* Project : sophAi
* Description :
* ========================================================================
      FUNCTION rnd(no_value)   
      implicit none
      include "mpif.h"
!This is the most uniform random number generator in numerical recipes F77
!The codes are modified from ran2
!It also include seed generator for self-generating random number
!The function out put value from 0.D0 to 1.D0 for double precision
!No need to input argument,just using ran()
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      integer no_value
      REAL*8 rnd,AM,EPS,RNMX
      REAL*8 seed_from_time,time_scale,seed_scale
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.D0/dble(IM1),
     *IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     *IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2D-7,
     *RNMX=1.D0-EPS,time_scale=1D3,seed_scale=13D4)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      seed_from_time=mpi_wtime()*time_scale
      idum=
     &-dint((seed_from_time-dble(dint(seed_from_time)))*seed_scale) !Self generated Seed number, must < 0 to re-initialized 
C      write(*,*)dble(seed_from_time),dble(dint(seed_from_time))  !for test
C     &,seed_from_time-dble(dint(seed_from_time)),idum
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1     
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      rnd=dmin1(AM*dble(iy),RNMX)
      return
      END

C      function rndi(rndi_range)  
C      implicit none
C      integer rndi,rndi_range !assign the range from 0 to rani_range
C      real*8 rnd
C      rndi=1+dint(dble(rndi_range)*rnd())
C      end

* ======================GNU General Public License=======================
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*  
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
* =======================================================================
