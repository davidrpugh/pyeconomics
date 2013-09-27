#if 0
#include <mcc.h>
#include <assert.h>
#endif
#include <mex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>


#define GET_ARG_MAT(x,i) x = mxGetPr(prhs[(i)-1])-1; x##_nr = mxGetM(prhs[(i)-1]); x##_nc = mxGetN(prhs[(i)-1])
#define GET_GLOBAL_MAT(x) x##_ = mexGetArray(#x,"global"); x##_nr = mxGetM(x##_); x##_nc = mxGetN(x##_); x = mxGetPr(x##_)-1

#if 1
#define DECL_MAT(x) mxArray *x##_; double *x; int x##_nr,x##_nc
#endif

#define CREATE_MAT(x,nr,nc) x##_ = mxCreateDoubleMatrix(nr,nc,mxREAL);x##_nr = mxGetM(x##_); x##_nc = mxGetN(x##_)

#define DECL_SCALAR(x) mxArray *x##_; double x
#define DECL_OUTPUT(x) double *x; int x##_nr,x##_nc

#define GET_GLOBAL_SCALAR(x) x##_ = mexGetArray(#x,"global"); x = mxGetScalar(x##_)


#define GET_ARG_SCALAR(i) mxGetScalar(prhs[(i)-1])
#define CREATE_OUTPUT(i,name,nr,nc,type) {plhs[(i)-1] = mxCreateDoubleMatrix(nr,nc, mx##type); name = mxGetPr(plhs[(i)-1])-1;name##_nr=nr;name##_nc=nc;}

#define MATEL(m,i,j) (m[(i) + ((j)-1)*m##_nr])
#define VECEL(m,i) (m[i])

#define NARGIN nrhs
#define NARGOUT nlhs

#define CHECK_ARGN(i,j) {\
  char str[100];\
  if(nrhs!=(i)) {sprintf(str,"%d inputs required.",i);mexErrMsgTxt(str);}\
  if(!(nlhs==(j)|| (j==1 && nlhs==0))) {sprintf(str,"%d outputs required.",j);mexErrMsgTxt(str);}\
}


/*  STANDARD DEFINITIONS:  */


#define assert(cond,str) if(!(cond)) {mexErrMsgTxt(#str);}
#define INCR(i,a,b) for((i)=(a);(i)<=(b);(i)++)
#define CHECKALLOC(m) {if (!(m)) {printf("allocation failure in matrix()\n");exit(1);}}
#define CREATE_VECTOR2D(name,type,nrow,ncols) \
{\
	int _i_;\
        name=(type **) malloc((nrow)*sizeof(type*)) - 1;\
	for(_i_=1;_i_<=nrow;_i_++) name[_i_]=(type *)malloc((ncols[_i_]+1)*sizeof(type));\
}
#define FREE_VECTOR2D(name,nrow)  {int _i_;\
	for(_i_=1;_i_<=nrow;_i_++) free(name[_i_]); free(name+1);}
#define NOTADMISSIBLE 1.76e125

#define VECTOR(type,n) (type *) malloc(((n)+1)*sizeof(type))
