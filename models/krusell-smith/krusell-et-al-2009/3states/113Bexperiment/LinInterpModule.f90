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

real(8) function lininterp3(x, y, z, xxgrid, yygrid, zzgrid, ZZ)

	implicit none

	integer :: nxx, nyy, nzz, xndx, yndx, zndx
	real(8) :: t, u, v, xyz
	real(8), intent(in) :: x, y, z, xxgrid(:), yygrid(:), zzgrid(:), ZZ(:,:,:)

	nxx = size(ZZ,1)
	nyy = size(ZZ,2)
	nzz = size(ZZ,3)

	xndx = locate(x, xxgrid)
	yndx = locate(y, yygrid)
	zndx = locate(z, zzgrid)

	if (xndx == 0 .or. xndx == nxx) then
		print*, "the first argment is out of range in lininterp3"
	else if (yndx == 0 .or. yndx == nyy) then
		print*, "the second argment is out of range in lininterp3"
	else if (zndx == 0 .or. zndx == nzz) then 
		print*, "the third argment is out of range in lininterp3"
	end if
    

	t = (x - xxgrid(xndx))/(xxgrid(xndx+1) - xxgrid(xndx))
	u = (y - yygrid(yndx))/(yygrid(yndx+1) - yygrid(yndx))
	v = (z - zzgrid(zndx))/(zzgrid(zndx+1) - zzgrid(zndx))

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

real(8) function lininterp4(aa, bb, cc, dd, aagrid, bbgrid, ccgrid, ddgrid,	ZZ)

	implicit none

	integer :: naa, nbb, ncc, ndd, andx, bndx, cndx, dndx
	real(8) :: s, t, u, v, val
	real(8), intent(in) :: aa, bb, cc, dd, aagrid(:), bbgrid(:), ccgrid(:),	&
						&	ddgrid(:), ZZ(:,:,:,:)

	naa = size(ZZ,1)
	nbb = size(ZZ,2)
	ncc = size(ZZ,3)
	ndd = size(ZZ,4)
		
	andx = locate(aa, aagrid)
	bndx = locate(bb, bbgrid)
	cndx = locate(cc, ccgrid)
	dndx = locate(dd, ddgrid)

	if (andx == 0 .or. andx == naa) then
		print*, "the first argment is out of range in lininterp4"
	else if (bndx == 0 .or. bndx == nbb) then
		print*, "the second argment is out of range in lininterp4"
	else if (cndx == 0 .or. cndx == ncc) then
		print*, "the third argment is out of range in lininterp4"
	else if (dndx == 0 .or. dndx == ndd) then
		print*, "the fourth argment is out of range in lininterp4"
	end if

	s = (aa - aagrid(andx))/(aagrid(andx+1) - aagrid(andx))
	t = (bb - bbgrid(bndx))/(bbgrid(bndx+1) - bbgrid(bndx))
	u = (cc - ccgrid(cndx))/(ccgrid(cndx+1) - ccgrid(cndx))
	v = (dd - ddgrid(dndx))/(ddgrid(dndx+1) - ddgrid(dndx))

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

real(8) function lininterp5(aa, bb, cc, dd, ee, aagrid, bbgrid, ccgrid, ddgrid,	&
							&	eegrid, ZZ)

	implicit none

	integer :: naa, nbb, ncc, ndd, nee, andx, bndx, cndx, dndx, endx
	real(8) :: s, t, u,	v, w, val
	real(8), intent(in) :: aa, bb, cc, dd, ee, aagrid(:), bbgrid(:), ccgrid(:),		&
						&	ddgrid(:), eegrid(:), ZZ(:,:,:,:,:)

	naa = size(ZZ,1)
	nbb = size(ZZ,2)
	ncc = size(ZZ,3)
	ndd = size(ZZ,4)
	nee = size(ZZ,5)
	
	andx = locate(aa, aagrid)
	bndx = locate(bb, bbgrid)
	cndx = locate(cc, ccgrid)
	dndx = locate(dd, ddgrid)
	endx = locate(ee, eegrid)

	if (andx == 0 .or. andx == naa) then
		print*, "the first argment is out of range in lininterp5"
	else if (bndx == 0 .or. bndx == nbb) then
		print*, "the second argment is out of range in lininterp5"
	else if (cndx == 0 .or. cndx == ncc) then
		print*, "the third argment is out of range in lininterp5"
	else if (dndx == 0 .or. dndx == ndd) then
		print*, "the fourth argment is out of range in lininterp5"
	else if (endx == 0 .or. endx == nee) then
		print*, "the fifth argment is out of range in lininterp5"
	end if

	s = (aa - aagrid(andx))/(aagrid(andx+1) - aagrid(andx))
	t = (bb - bbgrid(bndx))/(bbgrid(bndx+1) - bbgrid(bndx))
	u = (cc - ccgrid(cndx))/(ccgrid(cndx+1) - ccgrid(cndx))
	v = (dd - ddgrid(dndx))/(ddgrid(dndx+1) - ddgrid(dndx))
	w = (ee - eegrid(endx))/(eegrid(endx+1) - eegrid(endx))

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
		+ (1-v)*(u*(t*(s*ZZ(andx+1,bndx+1,cndx+1,dndx,endx)			&
		+ (1-s)*ZZ(andx,bndx+1,cndx+1,dndx,endx))					&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx+1,dndx,endx)					&
		+ (1-s)*ZZ(andx,bndx,cndx+1,dndx,endx)))					&
		+ (1-u)*(t*(s*ZZ(andx+1,bndx+1,cndx,dndx,endx)				&
		+ (1-s)*ZZ(andx,bndx+1,cndx,dndx,endx))						&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx,dndx,endx)					&
		+ (1-s)*ZZ(andx,bndx,cndx,dndx,endx)))))

	lininterp5 = val

end function


!	6 dimensional linear interpolation

real(8) function lininterp6(aa, bb, cc, dd, ee, ff, aagrid, bbgrid, ccgrid, ddgrid,	&
							&	eegrid, ffgrid, ZZ)

	implicit none

	integer :: naa, nbb, ncc, ndd, nee, nff, andx, bndx, cndx, dndx, endx, fndx
	real(8) :: s, t, u,	v, w, x, val
	real(8), intent(in) :: aa, bb, cc, dd, ee, ff, aagrid(:), bbgrid(:), ccgrid(:),		&
						&	ddgrid(:), eegrid(:), ffgrid(:), ZZ(:,:,:,:,:,:)

	naa = size(ZZ,1)
	nbb = size(ZZ,2)
	ncc = size(ZZ,3)
	ndd = size(ZZ,4)
	nee = size(ZZ,5)
	nff = size(ZZ,6)

	andx = locate(aa, aagrid)
	bndx = locate(bb, bbgrid)
	cndx = locate(cc, ccgrid)
	dndx = locate(dd, ddgrid)
	endx = locate(ee, eegrid)
	fndx = locate(ff, ffgrid)

	if (andx == 0 .or. andx == naa) then
		print*, "the first argment is out of range in lininterp6"
	else if (bndx == 0 .or. bndx == nbb) then
		print*, "the second argment is out of range in lininterp6"
	else if (cndx == 0 .or. cndx == ncc) then
		print*, "the third argment is out of range in lininterp6"
	else if (dndx == 0 .or. dndx == ndd) then
		print*, "the fourth argment is out of range in lininterp6"
	else if (endx == 0 .or. endx == nee) then
		print*, "the fifth argment is out of range in lininterp6"
	else if (fndx == 0 .or. endx == nff) then
		print*, "the sixth argment is out of range in lininterp6"
	end if

	s = (aa - aagrid(andx))/(aagrid(andx+1) - aagrid(andx))
	t = (bb - bbgrid(bndx))/(bbgrid(bndx+1) - bbgrid(bndx))
	u = (cc - ccgrid(cndx))/(ccgrid(cndx+1) - ccgrid(cndx))
	v = (dd - ddgrid(dndx))/(ddgrid(dndx+1) - ddgrid(dndx))
	w = (ee - eegrid(endx))/(eegrid(endx+1) - eegrid(endx))
	x = (ff - ffgrid(fndx))/(ffgrid(fndx+1) - ffgrid(fndx))

	val = x*(w*(v*(u*(t*(s*ZZ(andx+1,bndx+1,cndx+1,dndx+1,endx+1,fndx+1)	&
		+ (1-s)*ZZ(andx,bndx+1,cndx+1,dndx+1,endx+1,fndx+1))				&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx+1,dndx+1,endx+1,fndx+1)				&
		+ (1-s)*ZZ(andx,bndx,cndx+1,dndx+1,endx+1,fndx+1)))					&
		+ (1-u)*(t*(s*ZZ(andx+1,bndx+1,cndx,dndx+1,endx+1,fndx+1)			&
		+ (1-s)*ZZ(andx,bndx+1,cndx,dndx+1,endx+1,fndx+1))					&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx,dndx+1,endx+1,fndx+1)				&
		+ (1-s)*ZZ(andx,bndx,cndx,dndx+1,endx+1,fndx+1))))					&
		+ (1-v)*(u*(t*(s*ZZ(andx+1,bndx+1,cndx+1,dndx,endx+1,fndx+1)		&
		+ (1-s)*ZZ(andx,bndx+1,cndx+1,dndx,endx+1,fndx+1))					&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx+1,dndx,endx+1,fndx+1)				&
		+ (1-s)*ZZ(andx,bndx,cndx+1,dndx,endx+1,fndx+1)))					&
		+ (1-u)*(t*(s*ZZ(andx+1,bndx+1,cndx,dndx,endx+1,fndx+1)				&
		+ (1-s)*ZZ(andx,bndx+1,cndx,dndx,endx+1,fndx+1))					&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx,dndx,endx+1,fndx+1)					&
		+ (1-s)*ZZ(andx,bndx,cndx,dndx,endx+1,fndx+1)))))					&
		+ (1-w)*(v*(u*(t*(s*ZZ(andx+1,bndx+1,cndx+1,dndx+1,endx,fndx+1)		&
		+ (1-s)*ZZ(andx,bndx+1,cndx+1,dndx+1,endx,fndx+1))					&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx+1,dndx+1,endx,fndx+1)				&
		+ (1-s)*ZZ(andx,bndx,cndx+1,dndx+1,endx,fndx+1)))					&
		+ (1-u)*(t*(s*ZZ(andx+1,bndx+1,cndx,dndx+1,endx,fndx+1)				&
		+ (1-s)*ZZ(andx,bndx+1,cndx,dndx+1,endx,fndx+1))					&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx,dndx+1,endx,fndx+1)					&
		+ (1-s)*ZZ(andx,bndx,cndx,dndx+1,endx,fndx+1))))					&
		+ (1-v)*(u*(t*(s*ZZ(andx+1,bndx+1,cndx+1,dndx,endx+1,fndx+1)		&
		+ (1-s)*ZZ(andx,bndx+1,cndx+1,dndx,endx,fndx+1))					&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx+1,dndx,endx,fndx+1)					&
		+ (1-s)*ZZ(andx,bndx,cndx+1,dndx,endx,fndx+1)))						&
		+ (1-u)*(t*(s*ZZ(andx+1,bndx+1,cndx,dndx,endx,fndx+1)				&
		+ (1-s)*ZZ(andx,bndx+1,cndx,dndx,endx,fndx+1))						&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx,dndx,endx,fndx+1)					&
		+ (1-s)*ZZ(andx,bndx,cndx,dndx,endx,fndx+1))))))					&
	    + (1-x)*(w*(v*(u*(t*(s*ZZ(andx+1,bndx+1,cndx+1,dndx+1,endx+1,fndx)  &
		+ (1-s)*ZZ(andx,bndx+1,cndx+1,dndx+1,endx+1,fndx))					&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx+1,dndx+1,endx+1,fndx)				&
		+ (1-s)*ZZ(andx,bndx,cndx+1,dndx+1,endx+1,fndx)))					&
		+ (1-u)*(t*(s*ZZ(andx+1,bndx+1,cndx,dndx+1,endx+1,fndx)				&
		+ (1-s)*ZZ(andx,bndx+1,cndx,dndx+1,endx+1,fndx))					&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx,dndx+1,endx+1,fndx)					&
		+ (1-s)*ZZ(andx,bndx,cndx,dndx+1,endx+1,fndx))))					&
		+ (1-v)*(u*(t*(s*ZZ(andx+1,bndx+1,cndx+1,dndx,endx+1,fndx)			&
		+ (1-s)*ZZ(andx,bndx+1,cndx+1,dndx,endx+1,fndx))					&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx+1,dndx,endx+1,fndx)					&
		+ (1-s)*ZZ(andx,bndx,cndx+1,dndx,endx+1,fndx)))						&
		+ (1-u)*(t*(s*ZZ(andx+1,bndx+1,cndx,dndx,endx+1,fndx)				&
		+ (1-s)*ZZ(andx,bndx+1,cndx,dndx,endx+1,fndx))						&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx,dndx,endx+1,fndx)					&
		+ (1-s)*ZZ(andx,bndx,cndx,dndx,endx+1,fndx)))))						&
		+ (1-w)*(v*(u*(t*(s*ZZ(andx+1,bndx+1,cndx+1,dndx+1,endx,fndx)		&
		+ (1-s)*ZZ(andx,bndx+1,cndx+1,dndx+1,endx,fndx))					&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx+1,dndx+1,endx,fndx)					&
		+ (1-s)*ZZ(andx,bndx,cndx+1,dndx+1,endx,fndx)))						&
		+ (1-u)*(t*(s*ZZ(andx+1,bndx+1,cndx,dndx+1,endx,fndx)				&
		+ (1-s)*ZZ(andx,bndx+1,cndx,dndx+1,endx,fndx))						&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx,dndx+1,endx,fndx)					&
		+ (1-s)*ZZ(andx,bndx,cndx,dndx+1,endx,fndx))))						&
		+ (1-v)*(u*(t*(s*ZZ(andx+1,bndx+1,cndx+1,dndx,endx,fndx)			&
		+ (1-s)*ZZ(andx,bndx+1,cndx+1,dndx,endx,fndx))						&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx+1,dndx,endx,fndx)					&
		+ (1-s)*ZZ(andx,bndx,cndx+1,dndx,endx,fndx)))						&
		+ (1-u)*(t*(s*ZZ(andx+1,bndx+1,cndx,dndx,endx,fndx)					&
		+ (1-s)*ZZ(andx,bndx+1,cndx,dndx,endx,fndx))						&
		+ (1-t)*(s*ZZ(andx+1,bndx,cndx,dndx,endx,fndx)						&
		+ (1-s)*ZZ(andx,bndx,cndx,dndx,endx,fndx))))))





	lininterp6 = val

end function





end module