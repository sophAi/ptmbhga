        subroutine seed_io(kind,kinda,name2,io)
        include "mpif.h"
        integer io,type,nproc,myid
        character kind*2,kinda*2,kindab*2,name2*2
        common/kinda/ kinda
        common/kindb/ kindb
        common/kind/ kindab
        common/engtype/ type
        if(myid.eq.0)then
          if(type.eq.1)then
            kind=kinda
          else if(type.eq.2)then
            kind=kindab
          else if(type.eq.3)then
            kind="Be"
            kinda="Be"
          else if(type.eq.4)then
            kind="Al"
            kinda="Al"
          else if(type.eq.5)then
            kind="Mo"
            kinda="Mo"
          else if(type.eq.6)then
            kind="Cu"
            kinda="Cu"
          else if(type.eq.7)then
            kind="Cu"
            kinda="Cu"
          endif
        endif
        
