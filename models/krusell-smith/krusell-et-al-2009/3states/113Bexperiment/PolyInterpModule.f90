!	-----------------------------------------------------------------
!	File name: PolyInterpModule.f90
!
!	collection of polynomial interpolation routines from 
!	Numerical Recipe in C
!	-----------------------------------------------------------------


module PolyInterpModule

contains

!	1 dimensional polynomial interpolation

real(8) function polyinterp1(x, xa, ya)

	implicit none

	integer :: i, m, ns, n
	real(8), intent(in) :: x, xa(:), ya(:)
	real(8) :: den, dif, dift, ho, hp, w, y, dy, c(size(xa)), d(size(xa))

	n = size(xa)
	ns = 1
	dif = dabs(x - xa(1))

	do i = 1, n
		dift = dabs(x - xa(i))
		if (dift < dif) then
			ns = i
			dif = dift
		end if
		c(i) = ya(i)
		d(i) = ya(i)
	end do

	y = ya(ns)
	ns = ns - 1
	do m = 1, n-1
		do i = 1, n-m
			ho = xa(i) - x
			hp = xa(i+m) - x
			w = c(i+1) - d(i)
			den = ho - hp
			if (den == 0 )	pause 'Error in polyinterp1'
			den = w/den
			d(i) = hp*den
			c(i) = ho*den
		end do
		if (2*ns < n-m) then
			dy = c(ns+1)
		else
			dy = d(ns)
			ns = ns - 1
		end if
		y = y + dy
	end do

	polyinterp1= y

end function



!	2 dimensional polynomial interpolation

real(8) function polyinterp2(x1, x2, x1a, x2a, ya)

	implicit none

	integer :: n1, j
	real(8), intent(in) :: x1, x2, x1a(:), x2a(:), ya(:,:)
	real(8) :: ymtmp(size(x1a)), y

	n1 = size(x1a)
	
	do j = 1, n1
		ymtmp(j) = polyinterp1(x2, x2a, ya(j,:))
	end do
	y = polyinterp1(x1, x1a, ymtmp)

	polyinterp2 = y

end function


!	3 dimensional polynomial interpolation: recursively call polyinterp2

real(8) function polyinterp3(x1, x2, x3, x1a, x2a, x3a, ya)

	implicit none

	integer :: n1, n2, i, j
	real(8), intent(in) :: x1, x2, x3, x1a(:), x2a(:), x3a(:), ya(:,:,:)
	real(8) :: ymtmp(size(x1a),size(x2a)), y

	n1 = size(x1a)
	n2 = size(x2a)
	
	do i = 1, n1
		do j = 1, n2
			ymtmp(i,j) = polyinterp1(x3, x3a, ya(i,j,:))
		end do
	end do

	y = polyinterp2(x1, x2, x1a, x2a, ymtmp);

	polyinterp3 = y

end function


!	4 dimensional polynomial interpolation: recursively call polyinterp3

real(8) function polyinterp4(x1, x2, x3, x4, x1a, x2a, x3a, x4a, ya)

	implicit none

	integer :: n1, n2, n3, i, j, kk
	real(8), intent(in) :: x1, x2, x3, x4, x1a(:), x2a(:), x3a(:),		&
						&	x4a(:),	ya(:,:,:,:)
	real(8) :: ymtmp(size(x1a),size(x2a),size(x3a)), y

	n1 = size(x1a)
	n2 = size(x2a)
	n3 = size(x3a)
		
	do i = 1, n1
		do j = 1, n2
			do kk = 1, n3
				ymtmp(i,j,kk) = polyinterp1(x4, x4a, ya(i,j,kk,:))
			end do
		end do
	end do

	y = polyinterp3(x1, x2, x3, x1a, x2a, x3a, ymtmp)

	polyinterp4 = y

end function



!	5 dimensional polynomial interpolation: recursively call polyinterp4

real(8) function polyinterp5(x1, x2, x3, x4, x5, x1a, x2a, x3a, x4a, x5a, ya)

	implicit none

	integer :: n1, n2, n3, n4, i, j, kk, l
	real(8), intent(in) :: x1, x2, x3, x4, x5, x1a(:), x2a(:), x3a(:),	&
						&	x4a(:), x5a(:), ya(:,:,:,:,:)
	real(8) :: ymtmp(size(x1a),size(x2a),size(x3a),size(x4a)), y

	n1 = size(x1a)
	n2 = size(x2a)
	n3 = size(x3a)
	n4 = size(x4a)
	
	do i = 1, n1
		do j = 1, n2
			do kk = 1, n3
				do l = 1, n4
					ymtmp(i,j,kk,l) = polyinterp1(x5, x5a, ya(i,j,kk,l,:))
				end do
			end do
		end do
	end do

	y = polyinterp4(x1, x2, x3, x4, x1a, x2a, x3a, x4a, ymtmp)

	polyinterp5 = y

end function


!	6 dimensional polynomial interpolation: recursively call polyinterp5

real(8) function polyinterp6(x1, x2, x3, x4, x5, x6, x1a, x2a, x3a,	x4a,	&
							&	x5a, x6a, ya)

	implicit none

	integer :: n1, n2, n3, n4, n5, i, j, kk, l, m
	real(8), intent(in) :: x1, x2, x3, x4, x5, x6, x1a(:), x2a(:), x3a(:),	&
						&	x4a(:), x5a(:), x6a(:), ya(:,:,:,:,:,:)
	real(8) :: ymtmp(size(x1a),size(x2a),size(x3a),size(x4a),size(x5a)), y

	n1 = size(x1a)
	n2 = size(x2a)
	n3 = size(x3a)
	n4 = size(x4a)
	n5 = size(x5a)

	do i = 1, n1
		do j = 1, n2
			do kk = 1, n3
				do l = 1, n4
					do m = 1, n5
						ymtmp(i,j,kk,l,m) = polyinterp1(x6, x6a, ya(i,j,kk,l,m,:))
					end do
				end do
			end do
		end do
	end do
	
	y = polyinterp5(x1, x2, x3, x4, x5, x1a, x2a, x3a, x4a, x5a, ymtmp)

	polyinterp6 = y

end function

end module