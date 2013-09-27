/* LEAVE THIS LINE ALWAYS UNCHANGED:  */
#include "cfunc.h"


/* LEAVE THIS LINE ALWAYS UNCHANGED:  */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  int i,j,n;

  DECL_MAT(pos);
  DECL_MAT(frac);
  DECL_MAT(D);

  DECL_OUTPUT(D2);

  CHECK_ARGN(3,1);

  GET_ARG_MAT(pos,1);
  GET_ARG_MAT(frac,2);
  GET_ARG_MAT(D,3);


  /*  number of intervals destination distribution */
  n = pos_nr-1;

  CREATE_OUTPUT(1,D2,n,1,REAL);

  for(i=1;i<=n;i++){
    int p1 = (int) VECEL(pos,i);
    int p2 = (int) VECEL(pos,i+1);
    double f1 = VECEL(frac,i);
    double f2 = VECEL(frac,i+1);
    if(p1==p2){
      VECEL(D2,i) = VECEL(D,p1)*(f2-f1);
    } else {
      VECEL(D2,i) = VECEL(D,p1)*(1-f1);
      for(j=p1+1;j<p2;j++)
	VECEL(D2,i) += VECEL(D,j);
      VECEL(D2,i) += VECEL(D,p2)*f2;
    }
  }

}
