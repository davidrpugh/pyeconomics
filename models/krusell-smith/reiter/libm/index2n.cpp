/* index and ngrid have unit-offset;  */
#include "cfunc.h"
#include "interp.h"

int index2n(const int *index, const int *ngrid, int dimen)
{
  int n,i;
  n = index[dimen] - 1;
  for(i=dimen-1;i>=1;i--){
    n = ngrid[i]*n + index[i]-1;
  }
  n++;
  return n;
}
