	subroutine random(x)
        implicit none 
	integer j,r,atoms,seed
	real*8 x
        common /atoms /atoms
        common /seed /seed
	external rnd
        r=1 
c        do j=1,atoms*3
        x=r*rnd(seed)
c	end do
	return
	end
