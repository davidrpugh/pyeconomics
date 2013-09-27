!	----------------------------------------------------------------------
!	File name: Update.f90
!
!	
!	----------------------------------------------------------------------
!CALL RLSE (NOBS, Y, NIND, X, LDX, INTCEP, B, SST, SSE)

 ! Arguments:
 !    NOBS   - Number of observations.  (Input)
 !     Y      - Vector of length NOBS containing the dependent
 !            (response) variable.  (Input)
 !    NIND   - Number of independent (explanatory) variables.  (Input)
 !    X      - NOBS by NIND matrix containing the independent
 !             (explanatory) variables.  (Input)
 !    LDX    - Leading dimension of X exactly as specified in the
 !             dimension statement in the calling program.  (Input)
 !    INTCEP - Intercept option.  (Input)
 !              INTCEP  Action
 !                 0     An intercept is not in the model.
 !                 1     An intercept is in the model.
 !    B      - Vector of length INTCEP+NIND containing a
 !             least-squares solution for the regression coefficients.
 !             (Output)
 !             For INTCEP = 0, the fitted value for observation I is
 !                 B(1)*X(I,1)+B(2)*X(I,2)+...+B(NIND)*X(I,NIND)
 !             For INTCEP = 1, the fitted value for observation I is
 !                 B(1)+B(2)*X(I,1)+...+B(NIND+1)*X(I,NIND)
 !    SST    - Total sum of squares.  (Output)
 !             If INTCEP = 1, the total sum of squares is corrected for
 !             the mean.
 !    SSE    - Sum of squares for error.  (Output)

subroutine UpdateLOM()

	use Globals
	use Numerical_Libraries

	implicit none

	integer, parameter:: nobs = Nperiod-N1
	real(8) RHS(nobs,2)
	real(8) sst, sse
	
	
	!	construct common regressors

	RHS(:,1) = dlog(Kdata(N1+1:Nperiod))
	  

	!	law of motion for aggregate capital

	call drlse(nobs, dlog(Kdata(N1+2:Nperiod+1)), 1, RHS, nobs, 1, NewKcoef, sst, sse)

	print*, sst
	print*, sse

	open(1, file='Output\SST.txt', status='unknown')
    open(2, file='Output\SSE.txt', status='unknown')

	write(1, '(F12.6)') sst
	write(2, '(F12.6)') sse

end subroutine


