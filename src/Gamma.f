       SUBROUTINE MOMENT_GEN(X,NDIM,SEED)
       INTEGER NDIM,SEED,J1
       REAL*8 X(NDIM),RND,FAC
       PARAMETER(FAC=1.D0)
       DO J1=1,NDIM
         X(J1)=FAC*RND(SEED)
       ENDDO
       RETURN
       END

       FUNCTION MOMENTFUNC(NDIM,X)
       implicit double precision (A-H,O-Z)
       INTEGER NDIM,J1
       REAL*8 MOMENTFUNC,X(NDIM),BeGamma
c       COMMON /cgsw/ SW
*************************************
*  Put your energy subroutine here  *
*  This function is for derivative  *
*************************************
c       CALL MOMENT(NDIM,X,EREAL)
       MOMENTFUNC=BeGamma(X,NDIM)
       DO J1=1,NDIM
         write(*,*) J1,X(J1)
       ENDDO
      WRITE(*,*) "EREAL=",MOMENTFUNC
C        WRITE(*,*) "EVAP EREAL IN GPFUNC=",EREAL
c       MOMENTFUNC=EREAL
       END

*  Before use BeGamma must call ReaEn

        double precision function BeGamma(landa2,n)
        implicit double precision (A-H,O-Z)
        double precision landa(200),landa2(n)
        parameter (arange=0.d0,frange=15.d0)
        common /den/landa
        external del_denF
c        write(*,*) ' Now is in BeGamma !'
        do i=1,n
         landa(i)=landa2(i)
        enddo
c        write(*,*) landa2,n
        a=arange
        b=frange
        s=0.d0
        call qtrap(del_denF,a,b,s)
        BeGamma=DLOG(s)
c        write(*,*) BeGamma,s
C        if (BeGamma.lt.0.d0) then
C         BeGamma=1.d6
C         write(*,*) ' BeGamma is less than 0...'
C        endif
        return
        end

        subroutine ReaEn(L_n,upot,outx,io)
        implicit double precision(A-H,O-Z)
        double precision moment(200),outx(6000)
        integer pdb_rec,dim3,io,nproc,myid,t1,t2,t3,t4
        character name*6
        common nproc,myid
        common /Rea2/moment
        common /pdbrec/ pdb_rec,dim3
        common /name1/name
        if(io.eq.1)then
          dim3=1
c        pause ' You have "En.dat" ?'
c        write(*,*)
          open(UNIT=8,FILE='En.dat',status="old")
          read(8,*) L_n
c        write(*,*) 'Moment numbers =',L_n
          do i = 1, L_n
            read(8,*) moment(i)
c        write(*,*) 'moment(',i,')=',moment(i)
          enddo
          close(8)
          t1=(L_n)/1000
          t2=(L_n)/100
          t2=mod(t2,10)
          t3=(L_n)/10
          t3=mod(t3,10)
          t4=mod(((L_n)),10)
          name=char(t1+48)//char(t2+48)//char(t3+48)//char(t4+48)//"Mo"
c        write(*,*)
       endif
       if(io.eq.0)then
         open(9,file="landa",status="replace")
         write(*,*) "The components of the moment function:"
         write(9,*) L_n
c       write(3,*) 
         k=0
         write(*,*) 'The minimun value=',upot
         do i=1, L_n
           k=k+1
           write(9,*) outx(i)
           write(*,*) "moment",k,(outx(i))
         enddo
         write(9,*)
         write(9,*) 'The minimun value=',upot
       close(9)

       endif
        
        return
        end

        double precision function del_denF(x,L_n)
        implicit double precision (A-H,O-Z)
        double precision landa(200),moment(200)
        common /den/landa
        common /Rea2/moment
c        write(*,*) ' Now is in del_denF !'
        Ax=1.D0
        s=0.d0
        do k=1,L_n
         Ax=Ax*x
         s=s-landa(k)*(Ax-moment(k))
        enddo
        del_denF=DEXP(s)
c        write(*,*) 'del_denF,s=',del_denF,s
        return
        end

      SUBROUTINE qtrap(func,a,b,s)
      implicit double precision (a-h,o-z)
      EXTERNAL func
c      PARAMETER (EPS=1.e-6,JMAX=20)
      PARAMETER (EPS=1.e-6)
      common /quse/JMAX
c      write(*,*) 'JMAX=',JMAX
      olds=-1.e30
      do 11 j=1,JMAX
        call trapzd(func,a,b,s,j)
c        call trainf(func,b,1.d6,s2,j)
c        s=s1+s2
        if (abs(s-olds).lt.EPS*abs(olds)) return
        if (s.eq.0..and.olds.eq.0..and.j.gt.6) return
        olds=s
11    continue
c      pause 'too many steps in qtrap'
      END

      SUBROUTINE trapzd(func,a,b,s,n)
      implicit double precision (a-h,o-z)
      EXTERNAL func
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
c     write(*,*) it
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END


        subroutine ADFUNCLA(L_n,landa,dmom)
        implicit double precision (a-h,o-z)
        dimension dmom(L_n)
        double precision moment(200),landa(L_n)
        common /Rea2/moment
        im=1
        do i=1,L_n
         dmom(i)=moment(i)-aMeanEn(landa,i,L_n)
         if (dabs(dmom(i)).gt.1.d6) then
          im=0
          write(*,*) 'dmom is too big !',dmom(i)
         endif
c         dmom(i)=-dmom(i)
        enddo
        return
        end


        double precision function aMeanEn(landa,n,L_n)
        implicit double precision (A-H,O-Z)
        double precision moment(L_n),landa(L_n),landa2(200)
        parameter (arange=0.d0,frange=15.d0)
        external densityFxn,densityF
        common /den/landa2
        common /den2/n2
c        write(*,*) ' Now is in aMeanEn !' 
        do i=1,L_n
         landa2(i)=landa(i)
        enddo
        a=arange
        b=frange
        call qtrap(densityF,a,b,s1)
        n2=n
        call qtrap(densityFxn,a,b,s2)
        aMeanEn=s2/s1
        return
        end

        double precision function densityF(L_n,x)
        implicit double precision (A-H,O-Z)
        double precision landa(200)
        common /den/landa
c        write(*,*) ' Now is in densityF !'
        Ax=1.D0
        densityF=0.d0
        do k=1,L_n
         Ax=Ax*x
         densityF=densityF-landa(k)*Ax
        enddo
        densityF=DEXP(densityF)
        return
        end

        double precision function densityFxn(L_n,x)
        implicit double precision (A-H,O-Z)
        double precision landa(200)
        common /den/landa
        common /den2/n
c        write(*,*) ' Now is in densityFxn !'
        Ax=1.D0
        densityF=0.d0
        do k=1,L_n
         Ax=Ax*x
         densityF=densityF-landa(k)*Ax
        enddo
        xn=1.d0
        do i=1,n
         xn=xn*x
        enddo
        densityFxn=DEXP(densityF)*xn
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

