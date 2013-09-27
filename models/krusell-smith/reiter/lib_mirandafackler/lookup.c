#include "mex.h"
#include <math.h>
/* 
% LOOKUP  Performs a table lookup.
% USAGE:
%   ind=lookup(tabvals,x,endadj);
% INPUTS:
%   tabvals: a sorted vector of values
%   x: a matrix a values
%   endadj: a optional endpoint adjustment: 0, 1, 2 or 3.
% Returns a matrix of size(x) with element (i,j) equal to
%   max k: x(i,j)>=tabvals(k)
%
% Optional endpoint adjustments:
%   0: no adjustments
%   1: values of x < min(tabvals) will return 
%        length(tabvals=tabvals(1))
%   2: values of x > max(tabvals) will return 
%        m-length(tabvals=tabvals(end))
%   3: adjustments 1 and 2 will be performed

% Based on the HUNT algorithm in Press, et al, Numerical Recipes.
*/

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
   mxArray *pBa;
   double *table, *x, *ind, xi;
   int i, j, jlo, jhi, inc, n, m, p, numfirst, n1, n2;
    
   if (nrhs<2)
       mexErrMsgTxt("Two arguments must be passed");
   if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]))           
       mexErrMsgTxt("Input arguments of inproper type");
   if (mxIsSparse(prhs[0]) || mxIsSparse(prhs[1]))           
       mexErrMsgTxt("Input arguments of inproper type");
   if (nrhs>2) p=(int)*mxGetPr(prhs[2]);
   else p=0;

   n=mxGetNumberOfElements(prhs[0]);      /* # of elements in table */
   if (n==0) mexErrMsgTxt("Table must contain at least one element");
   m = mxGetNumberOfElements(prhs[1]);    /* # of elements in x */ 
   plhs[0]=mxDuplicateArray(prhs[1]);     /* mxArray for output */
   if (m==0) return;
   table=mxGetPr(prhs[0]);                /* pointer to table data  */ 
   x=mxGetPr(prhs[1]);                    /* pointer to values data */
   ind=mxGetPr(plhs[0]);                  /* pointer to output data */

   /* Lower endpoint adjustment */
   numfirst=1;
   while (numfirst<n && table[numfirst]==table[0]) numfirst++;

   /* Upper endpoint adjustment */
   if (p>=2) 
   {
     n--;
     while (table[n]==table[n-1]) n--;
   }
  
   n1=n-1;
   n2=n-2;

   /* handle 1-value lists separately */
   if (n-numfirst<1)
   {
     if (p==1 || p==3) for (i=0; i<m;i++) ind[i]=numfirst;  
     else for (i=0; i<m; i++)
     {
        if (table[0]<=x[i]) ind[i]=numfirst;  
        else ind[i]=0;
     }
   }
   else
   {
   jlo=0;
   for (i = 0; i < m; i++)
   {
     inc=1;
     xi=x[i];
     if (xi>=table[jlo])
     {
       jhi=jlo+1;
       while (xi>=table[jhi])
       {
         jlo=jhi;
         jhi+=inc;
         if (jhi>=n) { jhi=n; break; }
         else { inc += inc; }
       }
     }
     else
     {
       jhi=jlo;
       jlo--;
       while (xi<table[jlo])
       {
         jhi=jlo;
         jlo-=inc;
         if (jlo<0) { jlo=-1; break; }
         else { inc += inc; }
       }
     }
     while (jhi-jlo>1)
     {
       j=(jhi+jlo)/2;
       if (xi>=table[j]) jlo=j; 
       else jhi=j; 
     }
     ind[i]=jlo+1;
     if (jlo<0)
     {
       jlo=0; 
       if (p==1 || p==3) ind[i]=numfirst;
     }
     if (jlo==n1) jlo=n2;
   }
   }
}