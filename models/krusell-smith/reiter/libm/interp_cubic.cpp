#include <mex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

#include "interp.h"
#include "simplearr.h"
#define INCR(i,a,b) for((i)=(a);(i)<=(b);(i)++)

/* contains the routines on finding the interpolation weights  */


static
void IndexMatrix(SimpleArray<int> &indx, const SimpleArray<int> &nnodes)
{
  int dim = nnodes.size();
  myassert(dim>=1,"dim wrong");
  int i,j,l,n=1;
  INCR(i,1,dim)
    n*= nnodes(i);
  indx.resize(n,dim);

  INCR(j,1,dim)
    indx(1,j) = 1;

  INCR(i,2,n){
    /* advance indices:  */
    for(j=dim;j>=1;j--){
      if(indx(i-1,j)<nnodes(j)){
	indx(i,j) = indx(i-1,j) + 1;
	INCR(l,1,j-1)
	  indx(i,l)=indx(i-1,l);
	break;
      } else 
	indx(i,j)=1;
    }
  }
}

static
void QuadInterpWeights(double *wx, double y, double nonc)
{
  double t4 =   y *  y;
  double t12 =   nonc *  nonc;
  wx[0] =  (1 / (1 + nonc) * (-nonc * y + 2 * nonc - 3 * y + t4 + 2)) / 0.2e1;
  wx[1] =  (1 / (-1 + t12) * (y - 2) * y);
  wx[2] = - ((-nonc - 1 + y) * y / (-1 + nonc)) / 0.2e1;
}
static
void QuadInterpWeights(double *wx, double *dwx, double y, double nonc,
		       double dx)
{
  double t4 =   y *  y;
  double t12 =   nonc *  nonc;
  wx[0] =  (1 / (1 + nonc) * (-nonc * y + 2 * nonc - 3 * y + t4 + 2)) / 0.2e1;
  wx[1] =  (1 / (-1 + t12) * (y - 2) * y);
  wx[2] = - ((-nonc - 1 + y) * y / (-1 + nonc)) / 0.2e1;

  double t1 = 2 * y;
  double t3 = 1 /  dx;
  dwx[0] = - (1 / (1 + nonc) * t3 * (nonc + 3 - t1)) / 0.2e1;
  dwx[1] =  (2 / (-1 + t12) * (y - 1) * t3);
  dwx[2] =  (1 / (-1 + nonc) * t3 * (nonc + 1 - t1)) / 0.2e1;
}

static inline
void QuadInterpWeights(double *wx, double *dwx, const SimpleArray<double> &X,
		       double x, int jlo)
{
  double dx = (X(jlo+2)-X(jlo))/2;
  double y = (x-X(jlo)) / dx;
  double nonc = (X(jlo+1)-X(jlo))/dx - 1;
  if(dwx)
    QuadInterpWeights(wx, dwx, y, nonc, dx);
  else
    QuadInterpWeights(wx, y, nonc);
}

/*
  The following routine prepares the 1D-interpolation at a given point x
  the interpolation is linear in the nodes, using 1 to 4 nodes
  it returns the location of the first node;
  the size of the Array weights indicates the number of nodes used

  number of nodes:
  1 if x coincides with a grid point
  2 if X.size==2 or iLinear=1 (linear interp in x)
  3 if X.size==3 or point in extreme interval (quadr. in x)
  4 otherwise (cubic in x)
*/

static
int CubicApprPars(SimpleArray<double> &weights, 
		  int *nweights,
		  double x,
		  double *X0,
		  int n,
		  int iLinear)
{
  int jlo;
  *nweights = -1;
  SimpleArray<double> X(X0+1,n);
  hunt_double(X0, n, x, &jlo);
  if(jlo==0 || jlo>=n){
    if(jlo==1 && x==X(1)){
      *nweights = 1;
      weights(1) = 1;
      return 1;
    } else {
      fprintf(stderr,"Error in interpolation: %e is not between %e and %e\n",
	      x,X(1),X(n));
      myassert(0,"not between");
    }
  }
  if(x==X(jlo)){
    *nweights = 1;
    weights(1) = 1;
    return jlo;
  } 
  if(x==X(jlo+1)){
    *nweights = 1;
    weights(1) = 1;
    return jlo+1;
  } 
  if(iLinear || X.size()==2){
    *nweights = 2;
    double frac = (x-X(jlo))/(X(jlo+1)-X(jlo));
    weights(1) = 1-frac;
    weights(2) = frac;
    return jlo;
  }
  if(jlo==1 || jlo==n-1){
    *nweights = 3;
    if(jlo==n-1)
      jlo--;
    QuadInterpWeights(&weights(1), 0, X, x, jlo);
    return jlo;
  } else {
    *nweights = 4;
    double Wlow[3],Wupp[3];
    double frac = (x-X(jlo))/(X(jlo+1)-X(jlo));
    double wlow = 1-frac;
    double wupp = frac;
    QuadInterpWeights(Wlow, 0, X, x, jlo-1);
    QuadInterpWeights(Wupp, 0, X, x, jlo);

    weights(1) = wlow*Wlow[0];
    weights(2) = wlow*Wlow[1] + wupp*Wupp[0];
    weights(3) = wlow*Wlow[2] + wupp*Wupp[1];
    weights(4) =                wupp*Wupp[2];
    return jlo-1;
  }
}


int InterpWeightsTensor(int *NumbGrid,
			double *Probs,
			double *State,
			double **nodes,
			int *nPerDim,
			int dim,
			int iLinear) 
{
  SimpleArray<SimpleArray<double> > weights(dim);
  int i;
  for(i=1;i<=dim;i++){
    //weights(i).SetNew();
    myassert(weights(i).Owns()==0,"owns!");
    weights(i).resize(4);
  }
  SimpleArray<int> istart(dim),nnodes(dim);

  INCR(i,1,dim){
    double x = State[i];
    int n = nPerDim[i];
    double *X = nodes[i];
    if(x<X[1] || x>X[n]){
      myassert(0,"out of bounds");
    }
    istart(i) = CubicApprPars(weights(i), 
			      &nnodes(i),
			      x,
			      X,n,iLinear);
    // printf("%d  %d\n",istart(i),nnodes(i));
  }
  SimpleArray<int> indxmat;
  IndexMatrix(indxmat,nnodes);
  int nn = indxmat.size(1);
  SimpleArray<int> indx(dim);

  int lauf = 0;
  int ii,k;
  double sumprob = 0;
  INCR(ii,1,nn){
    double weight = 1;
    INCR(k,1,dim){
      int ik = indxmat(ii,k);
      weight*=weights(k)(ik);
      indx(k) = istart(k) + ik - 1;
    }
    if(weight){
      lauf++;
      NumbGrid[lauf] = index2n(indx.begin()-1,nPerDim,dim);
      Probs[lauf] = weight;
      sumprob += weight;
    }
  }

  myassert(fabs(1-sumprob)<1e-8,"probs dont add up");
  return lauf;
}
