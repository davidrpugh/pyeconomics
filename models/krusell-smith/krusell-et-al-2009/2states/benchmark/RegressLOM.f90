!	----------------------------------------------------------------------
!	File name: RegressLOM.f90
!
!	regresses simulated data to get new law of motion for aggregate
!	capital and price equations.
!	----------------------------------------------------------------------


subroutine RegressLOM()

	use Globals
	use Numerical_Libraries

	implicit none

	integer, parameter:: nobs = Nperiod-Nskip
	integer tg, tb
	real(8) RHSG(goodtime,1), RHSB(badtime,1), LHSG(goodtime), LHSB(badtime)
	real(8) sstG, sseG, sstB, sseB
	
	
	!	construct regressors
	tg=0
	tb=0
	do time = Nskip+1, Nperiod

	if (iZdata(time)==2) then
	tg=tg+1
	RHSG(tg,1) = dlog(Kdata(time))
	LHSG(tg) = dlog(Kdata(time+1))
	else
	tb=tb+1
	RHSB(tb,1) = dlog(Kdata(time))
	LHSB(tb) = dlog(Kdata(time+1))
	end if   
	end do

	!	law of motion for aggregate capital

	call drlse(goodtime, LHSG, 1, RHSG, goodtime, 1, NewKcoefG, sstG, sseG)
	call drlse(badtime, LHSB, 1, RHSB, badtime, 1, NewKcoefB, sstB, sseB)


	open(1, file='Output\SSTG.txt', status='unknown')
    open(2, file='Output\SSEG.txt', status='unknown')
	open(3, file='Output\SSTB.txt', status='unknown')
    open(4, file='Output\SSEB.txt', status='unknown')
	write(1, '(F12.6)') sstG
	write(2, '(F12.6)') sseG
	write(3, '(F12.6)') sstB
	write(4, '(F12.6)') sseB
end subroutine

