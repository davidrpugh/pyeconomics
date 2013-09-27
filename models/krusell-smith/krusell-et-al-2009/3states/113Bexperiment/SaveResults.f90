
subroutine SaveCoefficients()

	use Globals

	open(1, file='Output\Coefficients.txt', status='unknown')
	write(1, '(2F10.6)') Kcoef(1), Kcoef(2)

end subroutine	
	
	
subroutine SaveValueFunction()

	use Globals


	open(1, file='Output\Vl.txt', status='unknown')
	open(2, file='Output\ASlo.txt', status='unknown')
	open(3, file='Output\Vm.txt', status='unknown')
	open(4, file='Output\ASmi.txt', status='unknown')
	open(5, file='Output\Vh.txt', status='unknown')
	open(6, file='Output\AShi.txt', status='unknown')
	open(19, file='Output\ASlo31.txt', status='unknown')
	open(20, file='Output\ASlo2.txt', status='unknown')
	open(21, file='Output\ASmi31.txt', status='unknown')
	open(22, file='Output\ASmi2.txt', status='unknown')
	open(23, file='Output\AShi31.txt', status='unknown')
	open(24, file='Output\AShi2.txt', status='unknown')



    do ia = 1, na

		 write(19, '(F15.6)') (ASlo(ia,1,1,1,1,1,2)+ASlo(ia,3,1,1,1,1,2))/2.0	
		 write(20, '(F15.6)') ASlo(ia,2,1,1,1,1,2)
		 write(21, '(F15.6)') (ASmi(ia,1,1,1,1,1,2)+ASmi(ia,3,1,1,1,1,2))/2.0	
		 write(22, '(F15.6)') ASmi(ia,2,1,1,1,1,2)
		 write(23, '(F15.6)') (AShi(ia,1,1,1,1,1,2)+AShi(ia,3,1,1,1,1,2))/2.0	
		 write(24, '(F15.6)') AShi(ia,2,1,1,1,1,2)

	do iPge = 1, nP
	do iPbe = 1, nP		
 	do iPgs = 1, nP
	do iPbs = 1, nP	
	do iPbf = 1, nP	
	do ik = 1, nk

	

		 write(1, '(F15.6)') Vl(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik)
		 write(2, '(F15.6)') ASlo(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik)
		 write(3, '(F15.6)') Vm(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik)
		 write(4, '(F15.6)') ASmi(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik)
		 write(5, '(F15.6)') Vh(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik)
		 write(6, '(F15.6)') AShi(ia,iPge,iPbe,iPgs,iPbs,iPbf,ik)	

	end do
	end do	
    end do
	end do
	end do
	end do
	end do 


    open(7, file='Output\V2l.txt', status='unknown')
	open(8, file='Output\AS2l.txt', status='unknown')
    open(9, file='Output\V2m.txt', status='unknown')
	open(10, file='Output\AS2m.txt', status='unknown')   
    open(11, file='Output\V2h.txt', status='unknown')
	open(12, file='Output\AS2h.txt', status='unknown')
	open(13, file='Output\AS2l31.txt', status='unknown')
	open(14, file='Output\AS2l2.txt', status='unknown')
	open(15, file='Output\AS2m31.txt', status='unknown')
	open(16, file='Output\AS2m2.txt', status='unknown')
	open(17, file='Output\AS2h31.txt', status='unknown')
	open(18, file='Output\AS2h2.txt', status='unknown')

    do time = 1, N1

		 write(13, '(F15.6)') (AS2l(time,25,1,1,1,1,1)+AS2l(time,25,3,1,1,1,1))/2.0	
		 write(14, '(F15.6)') AS2l(time,25,2,1,1,1,1)
		 write(15, '(F15.6)') (AS2m(time,25,1,1,1,1,1)+AS2m(time,25,3,1,1,1,1))/2.0	
		 write(16, '(F15.6)') AS2m(time,25,2,1,1,1,1)
		 write(17, '(F15.6)') (AS2h(time,25,1,1,1,1,1)+AS2h(time,25,3,1,1,1,1))/2.0	
		 write(18, '(F15.6)') AS2h(time,25,2,1,1,1,1)


    do ia = 1, na
	do iPge = 1, nP
	do iPbe = 1, nP		
 	do iPgs = 1, nP
	do iPbs = 1, nP	
	do iPbf = 1, nP	


		 write(7, '(F15.6)') V2l(time, ia,iPge,iPbe,iPgs,iPbs,iPbf)
		 write(8, '(F15.6)') AS2l(time, ia,iPge,iPbe,iPgs,iPbs,iPbf)
		 write(9, '(F15.6)') V2m(time, ia,iPge,iPbe,iPgs,iPbs,iPbf)
		 write(10, '(F15.6)') AS2m(time, ia,iPge,iPbe,iPgs,iPbs,iPbf)
		 write(11, '(F15.6)') V2h(time, ia,iPge,iPbe,iPgs,iPbs,iPbf)
		 write(12, '(F15.6)') AS2h(time, ia,iPge,iPbe,iPgs,iPbs,iPbf)


    end do
	end do		
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
	open(10, file='Output\PanelBta.txt', status='unknown')
	open(11, file='Output\PanelPge.txt', status='unknown')
    open(12, file='Output\PanelPbe.txt', status='unknown')
	open(13, file='Output\PanelPgs.txt', status='unknown')
    open(14, file='Output\PanelPbs.txt', status='unknown')
	open(15, file='Output\PanelPbf.txt', status='unknown')
 

	write(1, '(F12.6)') Kdata
	write(2, '(F12.6)') Gini
	write(3, '(F12.6)') wealth1
	write(4, '(F12.6)') wealth5
	write(5, '(F12.6)') wealth10
	write(6, '(F12.6)') wealth20
	write(7, '(F12.6)') wealth30
	write(8, '(F12.6)') Negratio

	write(9, '(F12.6)') PanelA
	write(10, '(I1)') PanelBta
	write(11, '(F12.6)') PanelPge
	write(12, '(F12.6)') PanelPbe	
	write(13, '(F12.6)') PanelPgs
	write(14, '(F12.6)') PanelPbs
	write(15, '(F12.6)') PanelPbf
end subroutine