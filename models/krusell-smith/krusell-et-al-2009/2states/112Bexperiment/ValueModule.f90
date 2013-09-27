!	----------------------------------------------------------------------
!	File name : Value.f90
!	----------------------------------------------------------------------	*/


module ValueModule

	use Globals


contains


real(8) function Utility(ja)

	implicit none

	real(8), parameter:: tiny = 1.0D-20
	real(8) cm 
	real(8) ja
	
	cm = income-ja
	if (cm>0) then
	Utility = dlog(cm)
	else
	Utility = -1000
	end if

end function


real(8) function ValueF(ja)

	implicit none

	real(8) ja, EV1, EV1pr, EV1dpr

	call splint(agrid,EV(:,iPg,iPb,ibta,ik),EV2p(:,iPg,iPb,ibta,ik),ja,EV1,EV1pr,EV1dpr,na)

	ValueF = Utility(ja) + bta(ibta)*EV1

end function

real(8) function Utility2(ja)

	implicit none

	real(8), parameter:: tiny = 1.0D-20
	real(8) cm 
    real(8) ja

	
	cm = income-ja
	if (cm>0) then
	Utility2 = dlog(cm)
	else
	Utility2 = -1000
	end if
	


end function


real(8) function ValueF2(ja)

	implicit none


	real(8) ja, EV1, EV1pr, EV1dpr

	call splint(agrid,EV2(:,iPg,iPb,ibta),EV2pp(:,iPg,iPb,ibta),ja,EV1,EV1pr,EV1dpr,na)

	ValueF2 = Utility2(ja) + bta(ibta)*EV1


end function
end module
