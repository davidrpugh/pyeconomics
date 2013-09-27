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

	call splint(agrid,EV(:,ieps,ibta,ik,iz),EV2p(:,ieps,ibta,ik,iz),ja,EV1,EV1pr,EV1dpr,na)

	ValueF = Utility(ja) + bta(ibta)*EV1

end function


end module
