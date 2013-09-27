#include "cfunc.h"
#include "interp.h"


/* LEAVE THIS LINE ALWAYS UNCHANGED:  */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  int i, j, l, dim, nw, np, offs, *nnodes, nWeights=2,*numbs;
  double **nodes,*probs;

  DECL_MAT(x);
  DECL_MAT(Grid);
  DECL_MAT(nPerDim);
  DECL_SCALAR(iMultiLin);

  DECL_OUTPUT(NumbGrid);
  DECL_OUTPUT(Weights);

  CHECK_ARGN(4,2);


  GET_ARG_MAT(x,1);
  GET_ARG_MAT(Grid,2);
  GET_ARG_MAT(nPerDim,3);
  iMultiLin = GET_ARG_SCALAR(4);

  dim = x_nc;
  assert(dim==Grid_nc,"incompatible dimensions in interpolation");
  nodes = VECTOR(double*,dim);
  numbs = VECTOR(int,100);
  probs = VECTOR(double,100);
  nnodes = VECTOR(int,dim);
  INCR(i,1,dim){
    nodes[i] = &MATEL(Grid,1,i) - 1;
    nnodes[i] = (int) VECEL(nPerDim,i);
  }
  
  if(iMultiLin==0)
    nWeights = InterpWeightsSimplex(numbs,probs,
				    &MATEL(x,1,1) - 1,
				    nodes,
				    nnodes,
				    dim);
  else {
    nWeights = InterpWeightsTensor(numbs,probs,
				   &MATEL(x,1,1) - 1,
				   nodes,
				   nnodes,
				   dim,
				   0);
  }
  CREATE_OUTPUT(1,NumbGrid,nWeights,1,REAL);
  CREATE_OUTPUT(2,Weights,nWeights,1,REAL);
  INCR(i,1,nWeights){
    MATEL(NumbGrid,i,1) = numbs[i];
    MATEL(Weights,i,1) = probs[i];
  }
  free(numbs);
  free(probs);
  free(nodes);
  free(nnodes);
}

void piksr2(int n, double arr[], int brr[])
{
  int i,j;
  double a;
  int b;

  for (j=2;j<=n;j++) {
    a=arr[j];
    b=brr[j];
    i=j-1;
    while (i > 0 && arr[i] > a) {
      arr[i+1]=arr[i];
      brr[i+1]=brr[i];
      i--;
    }
    arr[i+1]=a;
    brr[i+1]=b;
  }
}


int InterpWeightsSimplex(int *NumbGrid,
			 double *Probs,
			 double *State,
			 double **nodes,
			 int *nPerDim,
			 int dim)
{
  int *index = VECTOR(int,dim);
  int *nPerDimCumul = VECTOR(int,dim);
  int *dorder = VECTOR(int,dim);
  double *frac = VECTOR(double,dim);
  int i,ii,nWeights;
  int findex;

  nPerDimCumul[1] = 1;
  INCR(i,2,dim){
    nPerDimCumul[i] = nPerDimCumul[i-1]*nPerDim[i-1];
  }
  INCR(i,1,dim){
    int n = nPerDim[i];
    double xi = State[i];
    dorder[i] = i;
    if(n==1){
      index[i] = 1;
      frac[i] = 1;
    } else {
      hunt_double(nodes[i],nPerDim[i],xi,&index[i]);
      ii = index[i];
      assert(ii>0 && ii<n, "wrong ii");
      index[i]++; /* for index2n later;  */
      frac[i] = (xi-nodes[i][ii]) / (nodes[i][ii+1]-nodes[i][ii]);
      if(!(frac[i]>=0 && frac[i]<=1)){
	if(frac[i]>-1e-10)
	  frac[i] = 0;
	else if(frac[i]<1+1e-10)
	  frac[i] = 1;
	else {
	  printf("xi = %e; ii = %d; frac[i] = %e\n",
		 xi,ii,frac[i]);
	  assert(0,"junk");
	}
      }
    }
  }

  piksr2(dim,frac,dorder);
  findex = index2n(index, nPerDim, dim); 
  /* printf("%d %d %d %d %d\n",index[1],index[2],nPerDim[1],nPerDim[2],findex); */
  NumbGrid[1] = findex;
  Probs[1] = 1;
  nWeights = dim+1;
  INCR(i,1,dim){
    if(frac[i]==1){/* all fracs left are =1, nothing to be done;  */
      /* printf("break: %d\n",i); */
      nWeights = i;
      break;
    }
    ii = dorder[i];
    if(nPerDim[ii]>1)
      findex = findex - nPerDimCumul[ii];
    NumbGrid[i+1] = findex;
    Probs[i+1] = 1-frac[i];
    Probs[i] -= 1-frac[i];
  }

  free(index);
  free(nPerDimCumul);
  free(dorder);
  free(frac);
  return nWeights;
}

void myassert(int cond,const char *msg)
{
  if(!cond)
    mexErrMsgTxt(msg);
}
