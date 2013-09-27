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



	kgrid(1)=11.1
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
	
	epsgrid(1) = 0.0D0
	epsgrid(2) = 0.0D0
	epsgrid(3) = 0.0D0
	epsgrid(4) = 1.0D0


	! save grids

	open(1, file='Output\epsgrid.txt', status='unknown')

	write(1,'(f12.6)') epsgrid
	

end subroutine
	



subroutine AggregateShockGrid()

	use Globals
	
	implicit none
	

	! zgrid: 3 states 1:B 2:D 3:G

	zgrid(1) = 0.99D0
	zgrid(2) = 0.99D0
	zgrid(3) = 1.01D0



	! save grids

	open(1, file='Output\zgrid.txt', status='unknown')

	write(1,'(f12.6)') zgrid



end subroutine

subroutine BetaGrid()

	use Globals
	
	implicit none
	

	bta(1) = 0.9823D0
	bta(2) = 0.9879D0
	bta(3) = 0.9935D0



	! save grids

	open(1, file='Output\bta.txt', status='unknown')

	write(1,'(f12.6)') bta

end subroutine


subroutine UiGrid()

	use Globals
	
	implicit none
	

	ui(1) = 0.038D0
	ui(2) = 0.391D0  ! the first period after being fired high benefits
	ui(3) = 0.391D0  
	ui(4) = 0.0D0
	


	! save grids

	open(1, file='Output\ui.txt', status='unknown')

	write(1,'(f12.6)') ui



end subroutine