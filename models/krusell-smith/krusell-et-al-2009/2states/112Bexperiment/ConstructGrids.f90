!	----------------------------------------------------------------------
!	File name: ConstructGrids.f90
!	----------------------------------------------------------------------



subroutine IndividualCapitalGrid()

	use Globals

	integer i
	real(8) step, lamax, lamin, lagrid(na)

	! step for grid construction
  
  	lamax=dlog(amax+5.0D0)
	lamin=dlog(amin+5.0D0)

	step=(lamax-lamin)/(na-1)
	lagrid(1)=lamin
	do i=2,na
	lagrid(i)=lagrid(i-1)+step
	end do

	do i=1,na
	agrid(i)=dexp(lagrid(i))-5.0D0
    end do
  



	! save grid

	open(1, file='Output\agrid.txt', status='unknown')
	write(1, '(f12.6)') agrid

end subroutine



subroutine AggregateCapitalGrid()

	use Globals

!	integer i
!	real(8) step

	! step for grid construction
    

!	step = (kmax - kmin)/(nk - 1)

	
	! aggregate capital grid

!	kgrid(1) = kmin
!	do i = 2, nk
!		kgrid(i) = kgrid(i-1) + step
!	end do

	kgrid(1)=11.10
	kgrid(2)=11.52
	kgrid(3)=11.94
	kgrid(4)=12.36

	kmax=kgrid(4)
	kmin=kgrid(1)
	kss=(kmax+kmin)/2.0



	! save grid

	open(1, file='Output\kgrid.txt', status='unknown')
	write(1, '(f12.6)') kgrid


end subroutine



subroutine IdiosyncraticShockGrid()

	use Globals
	
	implicit none
	
	epsgrid(1) = 0.0
	epsgrid(2) = 1.0


	! save grids

	open(1, file='Output\epsgrid.txt', status='unknown')

	write(1,'(f12.6)') epsgrid
	

end subroutine
	



subroutine AggregateShockGrid()

	use Globals
	
	implicit none
	

	! zgrid: 2 states

	zgrid(1) = 0.99
	zgrid(2) = 1.01



	! save grids

	open(1, file='Output\zgrid.txt', status='unknown')

	write(1,'(f12.6)') zgrid



end subroutine

subroutine BetaGrid()

	use Globals
	
	implicit none
	

	bta(1) = 0.9858
	bta(2) = 0.9894
	bta(3) = 0.9930



	! save grids

	open(1, file='Output\bta.txt', status='unknown')

	write(1,'(f12.6)') bta



end subroutine


subroutine PGrid()

	use Globals

	integer i


	PGgrid (1) = 0.00
    PGgrid (2) = 0.20
	PGgrid (3) = 0.40
    PGgrid (4) = 0.60
    PGgrid (5) = 0.80
	PGgrid (6) = 1.00
  
	
    PBgrid (1) = 0.00
    PBgrid (2) = 0.20
    PBgrid (3) = 0.40
    PBgrid (4) = 0.60
    PBgrid (5) = 0.80
	PBgrid (6) = 1.00

  

	! save grid

	open(1, file='Output\PGgrid.txt', status='unknown')
	write(1, '(f12.6)') PGgrid

	! save grid

	open(1, file='Output\PBgrid.txt', status='unknown')
	write(1, '(f12.6)') PBgrid

end subroutine