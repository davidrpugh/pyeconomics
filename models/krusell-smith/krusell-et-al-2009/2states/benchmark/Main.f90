!	----------------------------------------------------------------------
!	File name : Main.f90
!
!	Main program that solves the model with aggregate fluctuations and stochastic beta.
!	----------------------------------------------------------------------


program Fluctuations

	use Globals
	use dfport

	implicit none

	
	real(8) errK, errR, errW, err
	character(8) starttime, finishtime, solvetime, simultime
	character(10) startdate, finishdate, today


	!	open a log file

	open(0, file='Output\Logfile.txt', status = 'unknown')


	!	set the starting time

	call date(startdate)
	starttime = clock()


	!	initialize parameters and global variables

	call InitializeParameters()
	call InitializeCoefficients()


	!	construct grids

	call IndividualCapitalGrid()
	call IdiosyncraticShockGrid()
	call AggregateCapitalGrid()
	call AggregateShockGrid()
	call BetaGrid()


	!	initialize value function and decision rule for asset

	call IdiosyncraticTransition()
    call AggregateTransition()
	call BetaTransition()

    !   generate aggregate shock for simulation

	call AggregateShockSeries()
	call InitializeValue()


	do iterLOM = 1, maxiterLOM
     
	sfac=0.5


        print*, "*********iterLOM is*********",iterLOM
        print*, "KcoefG(1)=", KcoefG(1)
		print*, "KcoefG(2)=", KcoefG(2)
		print*, "KcoefB(1)=", KcoefB(1)
		print*, "KcoefB(2)=", KcoefB(2)

		
		!	solve value function

		call date(today)
		solvetime = clock()

		write(0, '(//A,I2)') "Coefficients of Equations in Iteration ", iterLOM
		write(0, '(/A, 2F10.5)') "Law of Motion for KapitalG:", KcoefG
		write(0, '(/A, 2F10.5)') "Law of Motion for KapitalB:", KcoefB
		write(0, '(//A,I2,A,A,A)') "SolveValueFunction in Iteration ",		&
								  iterLOM, " started at ", today, solvetime
	
		call SolveValueFunction()


		!	generate artificial data through simulation

		call date(today)
		simultime = clock()
		write(0, '(/A,I2,A,A,A)') "SimulateData in Iteration ",			&
								 iterLOM, " started at ", today, simultime

		call SimulateData()


		!	updata laws of motion through regressions of simulated data

		call RegressLOM()

		errK = maxval(dabs(NewKcoefG - KcoefG))+maxval(dabs(NewKcoefB - KcoefB))
		

		KcoefG = sfac*NewKcoefG + (1.0D0 - sfac)*KcoefG
		KcoefB = sfac*NewKcoefB + (1.0D0 - sfac)*KcoefB		

		!	convergence test

		if (errK < tol_LOM) then
			write(0, '(//A,I2)') "Error in LOM Iteration ", iterLOM
			write(0, '(/A,F10.6)') "errK = ", errK
			write(0, '(//A,I2)') "Law of Motion converged at iteration ", iterLOM
			exit 
		else
			write(0, '(//A,I2)') "Errors in LOM Iteration ", iterLOM
			write(0, '(/A,F10.6)')"errK = ", errK 
		end if

	end do


	!	final simulation

	call date(today)
	simultime = clock()
	write(0, '(//A,A,A)') "FinalSimulation started at ", today, simultime

	call FinalSimulation()


	!	set the finishing time

	call date(finishdate)
	finishtime = clock()

	
	!	print the starting and finishing times

	write(0, '(/)') 
	write(0, '(A,A,A)') "The program started at ", startdate, starttime
	write(0, '(A,A,A)') "The program finished at ", finishdate, finishtime

end program Fluctuations
