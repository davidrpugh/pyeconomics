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
	integer tg, tb, td
	real(8) RHSG(goodtime,1), RHSB(badtime,1), LHSG(goodtime), LHSB(badtime), RHSD(downtime,1),  LHSD(downtime)
	real(8) sstG, sseG, sstB, sseB, sstD, sseD, sigG, sigB, sigD
	
	
	!	construct regressors
	tg=0
	tb=0
	td=0
	do time = Nskip+1, Nperiod

	if (iZdata(time)==3) then
	tg=tg+1
	RHSG(tg,1) = dlog(Kdata(time))
	LHSG(tg) = dlog(Kdata(time+1))
	else if (iZdata(time)==1) then
	tb=tb+1
	RHSB(tb,1) = dlog(Kdata(time))
	LHSB(tb) = dlog(Kdata(time+1))
	else if (iZdata(time)==2) then
	td=td+1
	RHSD(td,1) = dlog(Kdata(time))
	LHSD(td) = dlog(Kdata(time+1))

	end if   
	end do

	!	law of motion for aggregate capital

	call drlse(goodtime, LHSG, 1, RHSG, goodtime, 1, NewKcoefG, sstG, sseG)
	call drlse(badtime, LHSB, 1, RHSB, badtime, 1, NewKcoefB, sstB, sseB)
	call drlse(downtime, LHSD, 1, RHSD, downtime, 1, NewKcoefD, sstD, sseD)

	sigG=dsqrt(sseG/real((goodtime-2)))*100.0
	sigB=dsqrt(sseB/real((badtime-2)))*100.0
	sigD=dsqrt(sseD/real((downtime-2)))*100.0



	open(1, file='Output\SSTG.txt', status='unknown')
    open(2, file='Output\SSEG.txt', status='unknown')
	open(3, file='Output\SSTB.txt', status='unknown')
    open(4, file='Output\SSEB.txt', status='unknown')
	open(5, file='Output\SSTD.txt', status='unknown')
    open(6, file='Output\SSED.txt', status='unknown')
    open(7, file='Output\sigG.txt', status='unknown')
	open(8, file='Output\sigB.txt', status='unknown')
    open(9, file='Output\sigD.txt', status='unknown')

	write(1, '(F12.6)') sstG
	write(2, '(F12.6)') sseG
	write(3, '(F12.6)') sstB
	write(4, '(F12.6)') sseB
	write(5, '(F12.6)') sstD
	write(6, '(F12.6)') sseD
	write(7, '(F12.6)') sigG
	write(8, '(F12.6)') sigB
	write(9, '(F12.6)') sigD


end subroutine

