
subroutine SaveCoefficients()

	use Globals

	open(1, file='Output\Coefficients.txt', status='unknown')
	write(1, '(2F10.6)') Kcoef(1), Kcoef(2)

end subroutine	
	
	
subroutine SaveValueFunction()

	use Globals


	open(1, file='Output\V.txt', status='unknown')
	open(2, file='Output\AS.txt', status='unknown')

    do ia = 1, na
	do iPg = 1, nP
	do iPb = 1, nP		
    do ibta = 1, nbta
	do ik = 1, nk

	

		 write(1, '(F15.6)') V(ia,iPg,iPb,ibta,ik)
		 write(2, '(F15.6)') AS(ia,iPg,iPb,ibta,ik)
		
    end do
	end do
	end do
	end do
	end do 


    open(3, file='Output\V2.txt', status='unknown')
	open(4, file='Output\AS2.txt', status='unknown')

    do time = 1, N1
    do ia = 1, na
	do iPg = 1, nP
	do iPb = 1, nP		
    do ibta = 1, nbta


	

		 write(3, '(F15.6)') V2(time, ia,iPg,iPb,ibta)
		 write(4, '(F15.6)') AS2(time, ia,iPg,iPb,ibta)
		
    end do
	end do
	end do
	end do
	end do 

end subroutine

subroutine SaveSimulatedTimeSeriesData()

	use Globals

	open(1, file='Output\Kdata.txt', status='unknown')
	open(2, file='Output\Gini.txt', status='unknown')
	open(3, file='Output\wealth1.txt', status='unknown')
	open(4, file='Output\wealth5.txt', status='unknown')
	open(5, file='Output\wealth10.txt', status='unknown')
	open(6, file='Output\wealth20.txt', status='unknown')  
	open(7, file='Output\wealth30.txt', status='unknown')
	open(8, file='Output\Negratio.txt', status='unknown')  

	open(9, file='Output\PanelA .txt', status='unknown')
	open(10, file='Output\PanelPg.txt', status='unknown')
    open(11, file='Output\PanelPb.txt', status='unknown')
	open(12, file='Output\PanelBta.txt', status='unknown')



	write(1, '(F12.6)') Kdata
	write(2, '(F12.6)') Gini
	write(3, '(F12.6)') wealth1
	write(4, '(F12.6)') wealth5
	write(5, '(F12.6)') wealth10
	write(6, '(F12.6)') wealth20
	write(7, '(F12.6)') wealth30
	write(8, '(F12.6)') Negratio

	write(9, '(F12.6)') PanelA
	write(10, '(F12.6)') PanelPg
	write(11, '(F12.6)') PanelPb
	write(12, '(I1)') PanelBta
	

end subroutine