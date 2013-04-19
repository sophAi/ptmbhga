        PROGRAM XYZTOPDB
        IMPLICIT NONE
        INTEGER J1,NT,J2,J3,J4,J5,J6
        REAL*8 ENERGY,X(3000),Y(3000),Z(3000)
        character *2 TYPE1
        character name1*6
        write(*,*) "Please input the *.xyz file name:(input *)"
        read(*,*) name1
        write(*,*) "Please input the initial number(default=0)"
        read(*,*) J2
        write(*,*) "To:"
        read(*,*) J3
        DO J4=J2,J3
          J5=J4/10
          J6=mod(J4,10)
          write(*,*)
          write(*,*) "Opening ",name1//char(J5+48)//char(J6+48),
     &".xyz file....."
          open(50,file=name1//char(J5+48)//char(J6+48)//
     &".xyz",status="old")
          open(51,file=name1//char(J5+48)//char(J6+48)//
     &".pdb",status="replace")
          open(52,file=name1//char(J5+48)//char(J6+48)//
     &"Cu.pdb",status="replace")
          read(50,*) NT,ENERGY
          write(51,*)'HEADER      ',name1,ENERGY
          write(52,*)'HEADER      ',name1,ENERGY
          do 5 J1=1,NT
            read(50,*) TYPE1,X(J1),Y(J1),Z(J1)
C            type1="CA"
C       write(51,*)'HEADER',name1,ENERGY
            if(TYPE1.eq."Cu")then
          write(52,33)'ATOM', J1,     TYPE1, 'UNK', J1,X(J1),Y(J1),Z(J1)
            endif
          write(51,33)'ATOM', J1,     TYPE1, 'UNK', J1,X(J1),Y(J1),Z(J1)
5         continue
33        FORMAT      (a4,I7, 2X,a2,2X,  a3, I6,4X,3(F8.3))
          write(51,34)'END   '
          write(52,34)'END   '
34        FORMAT(a6,3(I5))
          write(*,*) name1//char(J5+48)//char(J6+48),
     &".xyz has been successfully opened"
          write(*,*) "Transform to :"
          write(*,*) "==>> ",
     &name1//char(J5+48)//char(J6+48),".pdb"
          write(*,*) "==>> ",name1//char(J5+48)//char(J6+48),
     &"Cu.pdb"
          write(*,*) "Writing succeed!"
          write(*,*) 
          CLOSE(50)
          CLOSE(51)
          CLOSE(52)
        ENDDO
        STOP
        END

