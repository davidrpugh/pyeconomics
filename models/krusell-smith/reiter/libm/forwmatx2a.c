/* LEAVE THIS LINE ALWAYS UNCHANGED:  */
#include "cfunc.h"


#define ADD(fr_,to_,v_) VECEL(iTo,lauf) = to_; VECEL(iFr,lauf) = fr_;  VECEL(Val,lauf++) = v_


/* LEAVE THIS LINE ALWAYS UNCHANGED:  */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  int i,j,n,nnz,lauf;

  DECL_MAT(pos);
  DECL_MAT(frac);
  DECL_MAT(D);

  DECL_OUTPUT(iFr);
  DECL_OUTPUT(iTo);
  DECL_OUTPUT(Val);

  CHECK_ARGN(2,3);

  GET_ARG_MAT(pos,1);
  GET_ARG_MAT(frac,2);


  /*  number of intervals destination distribution */
  n = pos_nr-1;

  nnz = 2*n-1;  /* since pos(1) = 1 and pos(end)=n!!  */
  CREATE_OUTPUT(1,iFr,nnz,1,REAL);
  CREATE_OUTPUT(2,iTo,nnz,1,REAL);
  CREATE_OUTPUT(3,Val,nnz,1,REAL);

  lauf=1;
  for(i=1;i<=n;i++){
    int p1 = (int) VECEL(pos,i);
    int p2 = (int) VECEL(pos,i+1);
    double f1 = VECEL(frac,i);
    double f2 = VECEL(frac,i+1);
    if(p1==p2){
      ADD(p1,i,f2-f1);
    } else {
      ADD(p1,i,1-f1);
      for(j=p1+1;j<p2;j++){
	ADD(j,i,1);
      }
      ADD(p2,i,f2);
    }
  }

}
