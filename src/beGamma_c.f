       double precision function  BeGamma(x,ndim)     ! 2002/3/25
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       double precision x(ndim),landa(20)
       double precision moment(20)
       common /Rea1/arange,frange,drange,icutR
       common /Rea2/moment
       common /den/landa
       common /allf/L_n
       L_n = ndim
          do i=1,L_n
          landa(i)=x(i)
          enddo
       asum = 0.D0
          DO 100 i = 0, icutR-1
          xj = dble(i) *  drange + frange
          call makegama(x,xj,gama,ndim)
          xk = xj +  drange
          call makegama(x,xk,bgama,ndim)
          asum=asum+ drange *(gama+bgama)/2.d0
100       continue
       BeGamma=DLOG(asum)
c          do 102 j = 1, L_n
c          BeGamma=BeGamma+x(j)*moment(j)
c102       continue
104    return
       end

       subroutine  makegama(x,xj,gama,ndim)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION x(ndim)
       DOUBLE PRECISION moment(20)
c       common /allf/L_n
       common /Rea2/moment
       L_n=ndim
       Axj=1.D0
       gama = 0.d0
          do 11 k = 1, L_n
          Axj=Axj*xj
c          gama = gama - x(k) * Axj
          gama = gama - x(k) * (Axj-moment(k))
11        continue
       gama = DEXP(gama)
       return
       end
 
       subroutine answer(ndim,upot,outx)
       implicit double precision (a-h,o-z)
       dimension outx(6000)
       open(3,file="landa")
       write(*,*) "The coordinate of atoms:"
       write(3,*) ndim
c       write(3,*) 
       k=0 
       do i=1, ndim
         k=k+1
         write(3,*) outx(i)
         write(*,*) "moment",k,(outx(i))
       enddo 
       write(3,*)
       write(3,*) 'The min_pot=',upot
       close(3)
       return
       end

       subroutine ADFUNCLA(NDIM,X,GRAD)
       implicit double precision (a-h,o-z)
       double precision moment(20)
       dimension x(NDIM),GRAD(NDIM)
       common /Rea2/moment
       do i=1,NDIM
         GRAD(i)=moment(i)-aMeanEn(x,i)
       enddo
       return
       end

       double precision function aMeanEn(landa2,n2)
       implicit double precision (A-H,O-Z)
       double precision landa(20),landa2(20)
       common /den/landa
       common /den2/n
       external densityFxn,densityF
       n=n2
       do i=1,L_n
         landa(i)=landa2(i)
       enddo
       call myInt(densityF,aMeanEn0)
       call myInt(densityFxn,aMeanEn)
       aMeanEn=aMeanEn/aMeanEn0
       return
       end

       double precision function densityF(x)
       implicit double precision (A-H,O-Z)
       double precision landa(20),moment(20)
       common /den/landa
       common /allf/L_n
       common /Rea2/moment
c       write(*,*) ' Now is in densityF !'
       Ax=1.D0
       densityF=0.d0
       do k=1,L_n
         Ax=Ax*x
         densityF=densityF-landa(k)*Ax
       enddo
       densityF=DEXP(densityF)
       return
       end
       
       double precision function densityFxn(x)
       implicit double precision (A-H,O-Z)
       external densityF
       common /den2/n
       densityFxn=FwithXn(x,densityF,n)
       return
       end

       double precision function FwithXn(x,func,n)
       implicit double precision(A-H,O-Z)
       external func
       FwithXn=func(x)
       do i=1,n
         FwithXn=FwithXn*x
       enddo
       return
       end
       
       subroutine myInt(func,s)
       implicit double precision (a-h,o-z)
       common /den2/n
       common /Rea1/arange,frange,drange,icutR
       external func
       s=0.d0
       do i=1,icutR-1
         xi=dble(i)*drange+frange
         s1=func(xi)
         xj=xi+drange
         s2=func(xj)
         s=s+(s1+s2)/2.d0
       enddo
       end

