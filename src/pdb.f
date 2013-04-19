        subroutine pdb(ndim,X)
C        include "mpif.h"
        implicit none
        integer ndim,nproc,myid,J1,J2,J3
        real*8 X(30000),o_x(1000),o_y(1000),o_z(1000)
        real*8 h_x(1000),h_y(1000),h_z(1000)
        character pdb_name*6
        common nproc,myid
        common /pdbname1/ pdb_name
        common /pdb2/ o_x,o_y,o_z,h_x,h_y,h_z
1       format(A6,I5,A11,I4,A4,F8.3,F8.3,F8.3)
2       format(A6,I5,I5)
3       format(A6,I5,I5,I5)
4       format(A17)
5       format(A16)
6       format(A17)
7       format(A3)
        open(14,file=pdb_name//".pdb",status="unknown")
        write(14,4) "HEADER    PROTEIN"
        write(14,5) "COMPND    "//pdb_name
        write(14,6) "AUTHOR    P.J.Hsu"
        do J1=1,ndim/3
          J3=J1*3
          write(14,1) 
     &"ATOM  ",J1,"  CA  UNK  ",J1,"    ",X(J3-2),X(J3-1),X(J3)
        enddo
        J2=ndim/3
        do J1=1,ndim/3-2
          J2=J2+1
          write(14,1)
     &"ATOM  ",J2,"  O   UNK  ",J2,"    ",o_x(J1),o_y(J1),o_z(J1)
          J2=J2+1
          write(14,1)
     &"ATOM  ",J2,"  H   UNK  ",J2,"    ",h_x(J1),h_y(J1),h_z(J1)
        enddo       
C       do J1=1,ndim/3-1
C         if (J1.eq.1)then
C           write(14,2) "CONECT",J1,J1+1
C         else if (J1.eq.(ndim/3-1))then
C           write(14,3)"CONECT",J1,J1+1,J1-1
C           write(14,2) "CONECT",J1+1,J1
C         else
C           write(14,3)"CONECT",J1,J1+1,J1-1
C         endif
C       enddo
        write(14,7) "END"
        close(14)
        return
        end
