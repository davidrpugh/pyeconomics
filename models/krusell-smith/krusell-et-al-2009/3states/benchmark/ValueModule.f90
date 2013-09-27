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
	
	cm = max(income - ja, tiny)

	Utility = dlog(cm)


end function


real(8) function ValueF(ja)

	implicit none

	real(8) ja, EV1, EV1pr, EV1dpr

	call splint(agrid,EV(:,ieps,ibta,ik,iz),EV2p(:,ieps,ibta,ik,iz),ja,EV1,EV1pr,EV1dpr,na)

	ValueF = Utility(ja) + bta(ibta)*EV1

end function


real(8) function ValueFp(ja)

	implicit none

	real(8) ja, EV1, EV1pr, EV1dpr, cm
	real(8), parameter:: tiny = 1.0D-20

	call splint(agrid,EV(:,ieps,ibta,ik,iz),EV2p(:,ieps,ibta,ik,iz),ja,EV1,EV1pr,EV1dpr,na)

	cm = max(income - ja, tiny)
	ValueFp = -1.0/cm + bta(ibta)*EV1pr

end function

real(8) function ValueFp2(ja)

	implicit none

	real(8) ja, EV1, EV1pr, EV1dpr, cm
	real(8), parameter:: tiny = 1.0D-20

	call splint(agrid,EV(:,ieps,ibta,ik,iz),EV2p(:,ieps,ibta,ik,iz),ja,EV1,EV1pr,EV1dpr,na)

	cm = max(income - ja, tiny)
	ValueFp2 = -1.0/(cm*cm) + bta(ibta)*EV1dpr

end function
end module
