      subroutine file_io(atom_num,num.ens_num,min,name,io)
      include "mpif.h"
      integer ens_num,par_num,MAXGEN,samebc,num,cent2,asw,io
      integer MUPDATE,MAXIT,FALSEIT,type,rec_sw,recover_rec
      integer TAKESTEP,RESETSEED,ndim,naccept,recover_sw
      real*8 simp,cg,bc,accrat,qrate,tmp
      real*8 aa,bb,cc,dd,ee,ff,LL,ene,enef,enef2,dx,recover_ratio
      real*8 DGUESS,ASTEP,STEP,GTOL,rec_ratio
      char name*6
      common nproc,myid
      common /seeding/ num
      common /same1/ samebc,MAXGEN
      common /accuracy1/ LL
      common /accuracy2/ ene
      common /accuracy3/ enef,enef2
      common /accuracy4/ dx,bc
      common /accuracy5/ simp
      common /accuracy6/ cg
      common /anatsw/ asw
      common /bmin1/ DGUESS,ASTEP,STEP,GTOL
      common /bmin2/ MAXIT,MUPDATE,TAKESTEP,RESETSEED,FALSEIT
      common /engtype/ type
      common /cent/ cent2
      common /rec_sw1/ rec_sw
      common /rec_sw2/ rec_ratio
      common /recover1/ recover_sw,recover_rec
      common /recover2/ recover_ratio
      common /mpigen1/ aa,bb,cc,dd,ee,ff
      common /mpigen2/ par_num
      if(io.eq.1)then
        open(1,file="accuracy.dat",status="old")
        read(1,*) LL,ene,enef,enef2,simp,cg
        close(1)
        open(2,file="atom_num.dat")
        read(2,*) atom_num,num
        close(2)
        open(3,file="config.dat")
        read(3,*) ens_num,min,aa,bb,cc,dd,ee,ff,par_num,samebc,MAXGEN
        close(3)
        open(4,file="dfunc.dat",status="old")
        read(4,*) dx,bc,asw
        close(4)
        open(5,file="rec.dat",status="old")
        read(5,*) rec_sw,rec_ratio
        close(5)
        open(6,file="recover.dat",status="old")
        read(6,*) recover_sw,recover_rec,recover_ratio
        close(6)
        open(7,file="type.dat",status="old")
        read(7,*) type,cent2
        close(7)
        if (myid.eq.0)then
          write(*,*) "MPI-GA,VERSION,8.1"
          write(*,*) "@ You use ",nproc," node(s)"
          write(*,*) "@ ",atom_num," particles"
          write(*,*) "@ seeding number= ",num,",(0 for no seeding)"
          write(*,*) "@ ",ens_num," ensembles"
          if (min.eq.1)then
            write(*,*) "@ using Simplex minimization"
          else if(min.eq.2)then
            write(*,*) "@ using Conjugate Gradient minimization"
          else if(min.eq.3)then  
            write(*,*) "@ using Conjugate Gradient then Simplex"
          else if(min.eq.4)then
            write(*,*) "@ using Simplex then Conjugate Gradient"
          else
            write(*,*) "@ using LBFGS minimization,MAXIT=",MAXGEN
            if (TAKESTEP.eq.1)then
              write(*,*) "@ Take step!!"
            else
              write(*,*) "@ Don't take step!!"
            endif
            if (RESETSEED.eq.1)then
              write(*,*) "@ Reset to seeding number!!"
            else
              write(*,*) "@ Don't reset to seeding number!!"
            endif
          endif
          if (ga.eq.1)then
            write(*,*) "@ using NON-MPI GA "
          else
            write(*,*) "@ using MPI GA "
          endif
          write(*,*) "@ The ratio of GA:"
          write(*,*) "@ Inversion: ",aa
          write(*,*) "@ Arithmetic:",bb
          write(*,*) "@ Geometic:  ",cc
          write(*,*) "@ Crossing:  ",dd
          write(*,*) "@ 2-Point:   ",ee
          write(*,*) "@ Angular move:",ff
          write(*,*) "@",par_num," parents"
          write(*,*) "@ When ",samebc," lowest sets are found.."
          write(*,*) "@",MAXGEN," maximun of generations"
          if (rec_sw.eq.0)then
            write(*,*) "@ Normal recording"
          endif 
          if (rec_sw.eq.1)then
            write(*,*) "@ Record nothing!"
          endif
        endif
      else
        open(9,file=name//".inf",status="replace")
        write(9,*) "<<<Config information v8.0>>>"
        write(9,*) "@ ",atom_num," particles"
        write(9,*) "@ seeding number= ",num,",(0 for no seeding)"
        if (min.eq.1)then
          write(9,*) "@ using Simplex minimization"
        else if(min.eq.2)then  
          write(9,*) "@ using Conjugate Gradient minimization"
        else if(min.eq.3)then  
          write(9,*) "@ using Conjugate Gradient then Simplex"
        else if(min.eq.4)then
          write(9,*) "@ using Simplex then Conjugate Gradient"
        else
          write(9,*) "@ using LBFGS minimization,MAXIT",MAXIT
          if(TAKESTEP.eq.1)then
            write(9,*) "@ Take step!!"
          else
            write(9,*) "@ Don't take step!!"
          endif
          if(RESETSEED.eq.1)then
            write(9,*) "@ Reset to seeding number!!"
          else
            write(9,*) "@ Don't reset to seeding number!!"
          endif
          write(9,*) "@ Fail LBFGS=",FALSEIT
        endif
        write(9,*) "@ ",ens_num," ensemble"
        write(9,*) "@ The ratio of GA:"
        write(9,*) "@ Inversion: ",aa
        write(9,*) "@ Arithmetic:",bb
        write(9,*) "@ Geometic:  ",cc
        write(9,*) "@ Crossing:  ",dd
        write(9,*) "@ 2-Point:   ",ee
        write(9,*) "@ Angular Move:",ff
        write(9,*) "@",par_num," parents"
        write(9,*) "@ Same number=",samebc
        write(9,*) "@ Max gens=",MAXGEN
        write(9,*) "@ Width of initial box is",LL
        write(9,*) "@ BC in func=",enef
        write(9,*) "@ Accuracy parameter in func=",enef2
        write(9,*) "@ The energy sort accuracy is ",ene
        write(9,*) "@ The dX in dfunc is ",dx
        write(9,*) "@ BC in dfunc is ",bc
        write(9,*) "@ Tolerance of Simplex is ",simp
        write(9,*) "@ Tolerance of Conjugate Gradient is ",cg
        write(9,*) "@ Tolerance of LBFGS(GTOL)=",GTOL
        write(9,*) "@ MCSCRH UPDATE =",MUPDATE
        write(9,*) "@ Takestep in bmin:Astep=",ASTEP,",Step=",STEP
        write(9,*) "@ Input DGUESS=",DGUESS
        if(rec_sw.eq.0)then
          write(9,*) "@ Normal recording"
        endif
        if(rec_sw.eq.1)then
          write(9,*) "@ Step recording,recording ratio=",rec_ratio
        endif
        if(rec_sw.eq.2)then
          write(9,*) "@ Record nothing!"
        endif
        if(recover_rec.eq.1) then
          write(9,*) "@ Start backup,backup ratio=",recover_ratio
        endif
        write(9,*) "For ",name 
******INPUT THE POTENTIAL DETIALS********************
      endif  
      stop
