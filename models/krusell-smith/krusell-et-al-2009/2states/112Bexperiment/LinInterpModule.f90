!	-----------------------------------------------------------------
!	File name: LinInterpModule.f90
!
!	collection of linear interpolation routines
!	-----------------------------------------------------------------



module LinInterpModule

contains


! locate x in xgrid 

integer function locate(x, xgrid)

	implicit none

	integer :: n, ju, jm, jl, jndx
	real(8), intent(in) :: x, xgrid(:)
	
	n = size(xgrid)
	jl = 0
	ju = n+1
	do while (ju-jl > 1)
		jm = (ju + jl)/2
		if (x >= xgrid(jm)) then
			jl = jm
		else
			ju = jm
		end if
	end do

	if (x == xgrid(1)) then 
		jndx = 1
	else if (x == xgrid(n)) then
		jndx = n-1
	else 
		jndx = jl
	end if

	locate = jndx

end function


!	one dimensional linear interpolation 

real(8) function lininterp1(x, xgrid, ZZ)

	implicit none

	integer :: nx, xndx
	real(8) :: t, xz
	real(8), intent(in) :: x, xgrid(:), ZZ(:)
	

	nx = size(xgrid)
	xndx = locate(x, xgrid)
	
	if (xndx == 0 .or. xndx == nx) then
		print*, "out of bound in lininterp1"
	end if

	t = (x - xgrid(xndx))/(xgrid(xndx+1) - xgrid(xndx))
	xz = (1-t)*ZZ(xndx) + t*ZZ(xndx+1)

	lininterp1 = xz

end function


!	two dimensional interpolation

real(8) function lininterp2(x, y, xgrid, ygrid, ZZ)

	implicit none

	integer :: nx, ny, xndx, yndx
	real(8) :: t, u, xyz
	real(8), intent(in) :: x, y, xgrid(:), ygrid(:), ZZ(:,:)
	
	nx = size(ZZ,1)
	ny = size(ZZ,2)
	
	xndx = locate(x, xgrid)
	yndx = locate(y, ygrid)

	if (xndx == 0 .or. xndx == nx) then
		print*, "the first argment is out of range in lininterp2"
	else if (yndx == 0 .or. yndx == ny) then
		print*, "the second argment is out of range in lininterp2"
	end if

	t = (x - xgrid(xndx))/(xgrid(xndx+1) - xgrid(xndx))
	u = (y - ygrid(yndx))/(ygrid(yndx+1) - ygrid(yndx))

	xyz = (1-t)*(1-u)*ZZ(xndx,yndx) + t*(1-u)*ZZ(xndx+1,yndx)	&
		+ t*u*ZZ(xndx+1,yndx+1) + (1-t)*u*ZZ(xndx,yndx+1)

	lininterp2 = xyz

end function


!	3 dimensional linear interpolation
!	using continuous convex combination of neighboring grid points 

real(8) function lininterp3(x, y, z, xgrid, ygrid, zgrid, ZZ)

	implicit none

	integer :: nx, ny, nz, xndx, yndx, zndx
	real(8) :: t, u, v, xyz
	real(8), intent(in) :: x, y, z, xgrid(:), ygrid(:), zgrid(:), ZZ(:,:,:)

	nx = size(ZZ,1)
	ny = size(ZZ,2)
	nz = size(ZZ,3)

	xndx = locate(x, xgrid)
	yndx = locate(y, ygrid)
	zndx = locate(z, zgrid)

	if (xndx == 0 .or. xndx == nx) then
		print*, "the first argment is out of range in lininterp3"
	else if (yndx == 0 .or. yndx == ny) then
		print*, "the second argment is out of range in lininterp3"
	else if (zndx == 0 .or. zndx == nz) then 
		print*, "the third argment is out of range in lininterp3"
	end if
    

	t = (x - xgrid(xndx))/(xgrid(xndx+1) - xgrid(xndx))
	u = (y - ygrid(yndx))/(ygrid(yndx+1) - ygrid(yndx))
	v = (z - zgrid(zndx))/(zgrid(zndx+1) - zgrid(zndx))

	xyz = t*u*v*ZZ(xndx+1, yndx+1, zndx+1)				&
		+ (1-t)*u*v*ZZ(xndx, yndx+1, zndx+1)			&
		+ t*(1-u)*v*ZZ(xndx+1, yndx, zndx+1)			&
		+ (1-t)*(1-u)*v*ZZ(xndx, yndx, zndx+1)			&
		+ t*u*(1-v)*ZZ(xndx+1, yndx+1, zndx)			&
		+ (1-t)*u*(1-v)*ZZ(xndx, yndx+1, zndx)			&
		+ t*(1-u)*(1-v)*ZZ(xndx+1, yndx, zndx)			&
		+ (1-t)*(1-u)*(1-v)*ZZ(xndx, yndx, zndx)

	lininterp3 = xyz

end function


!	4 dimensional linear interpolation

real(8) function lininterp4(a, b, c, d, agrid, bgrid, cgrid, dgrid,	ZZ)

	implicit none

	integer :: na, nb, nc, nd, andx, bndx, cndx, dndx
	real(8) :: s, t, u, v, val
	real(8), intent(in) :: a, b, c, d, agrid(:), bgrid(:), cgrid(:),	&
						&	dgrid(:), ZZ(:,:,:,:)

	na = size(ZZ,1)
	nb = size(ZZ,2)
	nc = size(ZZ,3)
	nd = size(ZZ,4)
		
	andx = locate(a, agrid)
	bndx = locate(b, bgrid)
	cndx = locate(c, cgrid)
	dndx = locate(d, dgrid)

	if (andx == 0 .or. andx == na) then
		print*, "the first argment is out of range in lininterp4"
	else if (bndx == 0 .or. bndx == nb) then
		print*, "the second argment is out of range in lininterp4"
	else if (cndx == 0 .or. cndx == nc) then
		print*, "the third argment is out of range in lininterp4"
	else if (dndx == 0 .or. dndx == nd) then
		print*, "the fourth argment is out of range in lininterp4"
	end if

	s = (a - agrid(andx))/(agrid(andx+1) - agrid(andx))
	t = (b - bgrid(bndx))/(bgrid(bndx+1) - bgrid(bndx))
	u = (c - cgrid(cndx))/(cgrid(cndx+1) - cgrid(cndx))
	v = (d - dgrid(dndx))/(dgrid(dndx+1) - dgrid(dndx))

	val = v*(u*(t*(s*ZZ(andx+1,bndx+1,cndx+1,dndx+1)		&
		+ (1-s)*ZZ(andx,bndx+1,cndx+1,dndx+1))				&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx+1,dndx+1)			&
		+ (1-s)*ZZ(andx,bndx,cndx+1,dndx+1)))				&
		+ (1-u)*(t*(s*ZZ(andx+1,bndx+1,cndx,dndx+1)			&
		+ (1-s)*ZZ(andx,bndx+1,cndx,dndx+1))				&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx,dndx+1)				&
		+ (1-s)*ZZ(andx,bndx,cndx,dndx+1))))				&
		+ (1-v)*(u*(t*(s*ZZ(andx+1,bndx+1,cndx+1,dndx)		&
		+ (1-s)*ZZ(andx,bndx+1,cndx+1,dndx))				&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx+1,dndx)				&
		+ (1-s)*ZZ(andx,bndx,cndx+1,dndx)))					&
		+ (1-u)*(t*(s*ZZ(andx+1,bndx+1,cndx,dndx)			&
		+ (1-s)*ZZ(andx,bndx+1,cndx,dndx))					&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx,dndx)				&
		+ (1-s)*ZZ(andx,bndx,cndx,dndx))))

	lininterp4 = val

end function


!	5 dimensional linear interpolation

real(8) function lininterp5(a, b, c, d, e, agrid, bgrid, cgrid, dgrid,	&
							&	egrid, ZZ)

	implicit none

	integer :: na, nb, nc, nd, ne, andx, bndx, cndx, dndx, endx
	real(8) :: s, t, u,	v, w, val
	real(8), intent(in) :: a, b, c, d, e, agrid(:), bgrid(:), cgrid(:),		&
						&	dgrid(:), egrid(:), ZZ(:,:,:,:,:)

	na = size(ZZ,1)
	nb = size(ZZ,2)
	nc = size(ZZ,3)
	nd = size(ZZ,4)
	ne = size(ZZ,5)
	
	andx = locate(a, agrid)
	bndx = locate(b, bgrid)
	cndx = locate(c, cgrid)
	dndx = locate(d, dgrid)
	endx = locate(e, egrid)

	if (andx == 0 .or. andx == na) then
		print*, "the first argment is out of range in lininterp5"
	else if (bndx == 0 .or. bndx == nb) then
		print*, "the second argment is out of range in lininterp5"
	else if (cndx == 0 .or. cndx == nc) then
		print*, "the third argment is out of range in lininterp5"
	else if (dndx == 0 .or. dndx == nd) then
		print*, "the fourth argment is out of range in lininterp5"
	else if (endx == 0 .or. endx == ne) then
		print*, "the fifth argment is out of range in lininterp5"
	end if

	s = (a - agrid(andx))/(agrid(andx+1) - agrid(andx))
	t = (b - bgrid(bndx))/(bgrid(bndx+1) - bgrid(bndx))
	u = (c - cgrid(cndx))/(cgrid(cndx+1) - cgrid(cndx))
	v = (d - dgrid(dndx))/(dgrid(dndx+1) - dgrid(dndx))
	w = (e - egrid(endx))/(egrid(endx+1) - egrid(endx))

	val = w*(v*(u*(t*(s*ZZ(andx+1,bndx+1,cndx+1,dndx+1,endx+1)		&
		+ (1-s)*ZZ(andx,bndx+1,cndx+1,dndx+1,endx+1))				&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx+1,dndx+1,endx+1)				&
		+ (1-s)*ZZ(andx,bndx,cndx+1,dndx+1,endx+1)))				&
		+ (1-u)*(t*(s*ZZ(andx+1,bndx+1,cndx,dndx+1,endx+1)			&
		+ (1-s)*ZZ(andx,bndx+1,cndx,dndx+1,endx+1))					&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx,dndx+1,endx+1)				&
		+ (1-s)*ZZ(andx,bndx,cndx,dndx+1,endx+1))))					&
		+ (1-v)*(u*(t*(s*ZZ(andx+1,bndx+1,cndx+1,dndx,endx+1)		&
		+ (1-s)*ZZ(andx,bndx+1,cndx+1,dndx,endx+1))					&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx+1,dndx,endx+1)				&
		+ (1-s)*ZZ(andx,bndx,cndx+1,dndx,endx+1)))					&
		+ (1-u)*(t*(s*ZZ(andx+1,bndx+1,cndx,dndx,endx+1)			&
		+ (1-s)*ZZ(andx,bndx+1,cndx,dndx,endx+1))					&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx,dndx,endx+1)					&
		+ (1-s)*ZZ(andx,bndx,cndx,dndx,endx+1)))))					&
		+ (1-w)*(v*(u*(t*(s*ZZ(andx+1,bndx+1,cndx+1,dndx+1,endx)	&
		+ (1-s)*ZZ(andx,bndx+1,cndx+1,dndx+1,endx))					&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx+1,dndx+1,endx)				&
		+ (1-s)*ZZ(andx,bndx,cndx+1,dndx+1,endx)))					&
		+ (1-u)*(t*(s*ZZ(andx+1,bndx+1,cndx,dndx+1,endx)			&
		+ (1-s)*ZZ(andx,bndx+1,cndx,dndx+1,endx))					&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx,dndx+1,endx)					&
		+ (1-s)*ZZ(andx,bndx,cndx,dndx+1,endx))))					&
		+ (1-v)*(u*(t*(s*ZZ(andx+1,bndx+1,cndx+1,dndx,endx+1)		&
		+ (1-s)*ZZ(andx,bndx+1,cndx+1,dndx,endx))					&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx+1,dndx,endx)					&
		+ (1-s)*ZZ(andx,bndx,cndx+1,dndx,endx)))					&
		+ (1-u)*(t*(s*ZZ(andx+1,bndx+1,cndx,dndx,endx)				&
		+ (1-s)*ZZ(andx,bndx+1,cndx,dndx,endx))						&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx,dndx,endx)					&
		+ (1-s)*ZZ(andx,bndx,cndx,dndx,endx)))))

	lininterp5 = val

end function


!	for int6lin val = x*(val in int5lin with sixth index fndx+1)
!				+ (1-x)*(val in int5lin with sixth index fndx)

end module