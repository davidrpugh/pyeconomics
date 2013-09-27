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

	if (ibta==1) then
	call splint(agrid,EVl(:,iPge,iPbe,iPgs,iPbs,iPbf,ik),EV2pl(:,iPge,iPbe,iPgs,iPbs,iPbf,ik),ja,EV1,EV1pr,EV1dpr,na)
	else if (ibta==2) then
	call splint(agrid,EVm(:,iPge,iPbe,iPgs,iPbs,iPbf,ik),EV2pm(:,iPge,iPbe,iPgs,iPbs,iPbf,ik),ja,EV1,EV1pr,EV1dpr,na)
	else
	call splint(agrid,EVh(:,iPge,iPbe,iPgs,iPbs,iPbf,ik),EV2ph(:,iPge,iPbe,iPgs,iPbs,iPbf,ik),ja,EV1,EV1pr,EV1dpr,na)
	end if


	ValueF = Utility(ja) + bta(ibta)*EV1

end function


real(8) function ValueFp(ja)

	implicit none

	real(8) ja, EV1, EV1pr, EV1dpr, cm
	real(8), parameter:: tiny = 1.0D-20

	if (ibta==1) then
	call splint(agrid,EVl(:,iPge,iPbe,iPgs,iPbs,iPbf,ik),EV2pl(:,iPge,iPbe,iPgs,iPbs,iPbf,ik),ja,EV1,EV1pr,EV1dpr,na)
	else if (ibta==2) then
	call splint(agrid,EVm(:,iPge,iPbe,iPgs,iPbs,iPbf,ik),EV2pm(:,iPge,iPbe,iPgs,iPbs,iPbf,ik),ja,EV1,EV1pr,EV1dpr,na)
	else
	call splint(agrid,EVh(:,iPge,iPbe,iPgs,iPbs,iPbf,ik),EV2ph(:,iPge,iPbe,iPgs,iPbs,iPbf,ik),ja,EV1,EV1pr,EV1dpr,na)
	end if

	cm = max(income - ja, tiny)
	ValueFp = -1.0/cm + bta(ibta)*EV1pr

end function

real(8) function ValueFp2(ja)

	implicit none

	real(8) ja, EV1, EV1pr, EV1dpr, cm
	real(8), parameter:: tiny = 1.0D-20

	if (ibta==1) then
	call splint(agrid,EVl(:,iPge,iPbe,iPgs,iPbs,iPbf,ik),EV2pl(:,iPge,iPbe,iPgs,iPbs,iPbf,ik),ja,EV1,EV1pr,EV1dpr,na)
	else if (ibta==2) then
	call splint(agrid,EVm(:,iPge,iPbe,iPgs,iPbs,iPbf,ik),EV2pm(:,iPge,iPbe,iPgs,iPbs,iPbf,ik),ja,EV1,EV1pr,EV1dpr,na)
	else
	call splint(agrid,EVh(:,iPge,iPbe,iPgs,iPbs,iPbf,ik),EV2ph(:,iPge,iPbe,iPgs,iPbs,iPbf,ik),ja,EV1,EV1pr,EV1dpr,na)
	end if

	cm = max(income - ja, tiny)
	ValueFp2 = -1.0/(cm*cm) + bta(ibta)*EV1dpr

end function




real(8) function Utility2(ja)

	implicit none

	real(8), parameter:: tiny = 1.0D-20
	real(8) cm 
    real(8) ja

	
	cm = max(income - ja, tiny)

	Utility2 = dlog(cm)



end function


real(8) function ValueF2(ja)

	implicit none


	real(8) ja, EV1, EV1pr, EV1dpr
	
	call splint(agrid,EV2(:,iPge,iPbe,iPgs,iPbs,iPbf,ibta),EV2pp(:,iPge,iPbe,iPgs,iPbs,iPbf,ibta),ja,EV1,EV1pr,EV1dpr,na)


	ValueF2 = Utility2(ja) + bta(ibta)*EV1


end function


real(8) function ValueF2p(ja)

	implicit none

	real(8) ja, EV1, EV1pr, EV1dpr, cm
	real(8), parameter:: tiny = 1.0D-20

	call splint(agrid,EV2(:,iPge,iPbe,iPgs,iPbs,iPbf,ibta),EV2pp(:,iPge,iPbe,iPgs,iPbs,iPbf,ibta),ja,EV1,EV1pr,EV1dpr,na)

	cm = max(income - ja, tiny)
	ValueF2p = -1.0/cm + bta(ibta)*EV1pr

end function

real(8) function ValueF2p2(ja)

	implicit none

	real(8) ja, EV1, EV1pr, EV1dpr, cm
	real(8), parameter:: tiny = 1.0D-20

	call splint(agrid,EV2(:,iPge,iPbe,iPgs,iPbs,iPbf,ibta),EV2pp(:,iPge,iPbe,iPgs,iPbs,iPbf,ibta),ja,EV1,EV1pr,EV1dpr,na)

	cm = max(income - ja, tiny)
	
	ValueF2p2 = -1.0/(cm*cm) + bta(ibta)*EV1dpr

end function

end module
