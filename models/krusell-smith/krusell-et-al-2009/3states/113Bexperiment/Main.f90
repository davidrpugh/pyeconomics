!	----------------------------------------------------------------------
!	File name : Main.f90
!
!	----------------------------------------------------------------------


program Transition

	use Globals
	use dfport

	implicit none

	character(8) starttime, finishtime, solvetime, simultime
	character(10) startdate, finishdate, today
	integer iterK, tim
    real(8) ErrK 
   

    call date(startdate)
	starttime = clock()

!   Initialize parameters

	call InitializeParameters()
	call InitializeK()

!   Construct grids	
  
	call IndividualCapitalGrid()
    call AggregateCapitalGrid()
	call BetaGrid()
    call PrGrid()
	call IdiosyncraticShockGrid()
	call AggregateShockGrid()


!   Transitions	
	call IdiosyncraticTransition()
	call AggregateTransition()
    call BetaTransition()
	
	
!   Initialize
    call GenerateSeries()
!	call InitializeValue()
    call InitializeValueRead()
	call InitializeValue2()
	call Computegbar()



	open(0, file='Output\Logfile.txt', status = 'unknown')

!   main loop
	do iterLOM = 1, maxiterLOM
	 if (iterLOM<10) then
	  sfac=0.20
	  else if (iterLOM<18) then
	  sfac=0.40
	  else
	  sfac=0.20
	  end if
	

	write(0, '(//A,I2)') "Coefficients of Equations in Iteration ", iterLOM
	write(0, '(/A, 2F10.5)') "Law of Motion for Kapital:", Kcoef
	write(0, '(//A,I2,A,A,A)') "SolveValueFunction in Iteration ",		&
								iterLOM, " started at ", today, solvetime


!   Solve value function
	call date(today)
	solvetime = clock()

    call SolveValueFunction()
	call IterateBack()
    call SimulateData()
    
	Kdata(1) = K(1)

	ErrK = maxval(dabs(K-Kdata))
	errKdata = dabs(K-Kdata)
	
	open(1, file='Output\errKdata.txt', status='unknown')
	close(1, status='delete')

	open(1, file='Output\errKdata.txt', status='unknown')
    write(1, '(f12.6)') errKdata

     
	  if (ErrK>TolErrK) then
			K = Kdata*sfac+K*(1.0D0-sfac)
			write(0, '(//A,I2)') "Errors in K Iteration ", iterLOM
			write(0, '(/A,F10.6)')"ErrK = ", ErrK 
            open(1, file='Output\K.txt', status='unknown')
        	write(1, '(F15.6)') K

		else

			write(0, '(//A,I2)') "Error in K Iteration ", iterLOM
			write(0, '(/A,F10.6)') "ErrK = ", ErrK
			write(0, '(//A,I2)') "K series converged at iteration ", iterLOM
			exit 
		end if

	 call UpDateLOM()
	 Kcoef = sfac*NewKcoef + (1.0D0 - sfac)*Kcoef
	 
	 do tim=1, N1	   
	  irate2(tim) = zseries(tim)*alpha*((K(tim)/((1.0-useries(tim))*hour))**(alpha-1.0))-delta
	  wage2(tim) = zseries(tim)*(1-alpha)*((K(tim)/((1.0-useries(tim))*hour))**(alpha))
	 end do

	end do


   	!	final simulation
	call date(today)
	simultime = clock()
	write(0, '(//A,A,A)') "FinalSimulation started at ", today, simultime

	call FinalSimulation()
	call ComputeWelfare()
	


	!	set the finishing time

	call date(finishdate)
	finishtime = clock()

	
	!	print the starting and finishing times

	write(0, '(/)') 
	write(0, '(A,A,A)') "The program started at ", startdate, starttime
	write(0, '(A,A,A)') "The program finished at ", finishdate, finishtime

end Program Transition
