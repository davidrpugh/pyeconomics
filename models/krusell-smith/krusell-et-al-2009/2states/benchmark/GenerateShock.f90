!	----------------------------------------------------------------------
!	File name: GenerateShock.f90
!   
!   Subroutines to generate aggregate and idiosyncratic (two kinds of) shocks
!	----------------------------------------------------------------------


subroutine AggregateShockSeries()

	use Numerical_Libraries
	use Globals

	integer, parameter:: zseed = 1
	integer jz, kz
	real(8) ushock


	! set the seed number for random number generator

	call rnset(zseed)


	! generate Zdata
       
	! The first period is always a good period 
	   iZdata(1) = nz
	   Zdata(1) = zgrid(nz)
	   jz = nz

	goodtime=0
	badtime=0

    ! Generate until Nperiod+1 for the use in the NextIdiosyncraticShock
	do time = 2, Nperiod+1

		ushock = drnunf()

		do kz = 1, nz

			if (ushock <= CtrZ(jz,kz)) then
				jz = kz
				exit
			end if

		end do

		iZdata(time) = jz
		Zdata(time) = zgrid(jz)
		if (jz==1 .and. time>Nskip .and. time<Nperiod+1) then
		badtime=badtime+1
		else if (jz==2 .and. time>Nskip .and. time<Nperiod+1) then
		goodtime=goodtime+1
		end if

	end do

	open(1, file='Output\Zdata.txt', status='replace')
	write(1, '(f12.6)') Zdata


end subroutine




!Given z, z', generate epsilon'. 
!
!**Be very careful to give an appropriate initial epsilon in the simulation(CrEps(indiv,1)). 
!  It has to be conditioned on the fact that the first period is always good

subroutine NextIdiosyncraticShock(CrEps)

	use Numerical_Libraries
	use Globals

	integer j, CrEps(Nindiv,2)

	real(8) ushock


	do indiv = 1, Nindiv

		ushock = drnunf()

		do j = 1, neps

			if (ushock <= CtrEps(iZdata(time),iZdata(time+1),CrEps(indiv,1),j)) then
				CrEps(indiv,2) = j
				exit
			end if
		
		end do

	end do


end subroutine

subroutine NextBeta(CrBta)

	use Numerical_Libraries
	use Globals

	integer j, CrBta(Nindiv,2)
	real(8) ushock


	do indiv = 1, Nindiv

		ushock = drnunf()
        
		do j = 1, nbta

			if (ushock <= CtrBta(CrBta(indiv,1),j)) then
				CrBta(indiv,2) = j
				exit
			end if
		
		end do

	end do


end subroutine
