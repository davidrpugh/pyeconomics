#include "cfunc.h"
#include "interp.h"

void hunt_double(const double *xx, int n, double x, int *jlo)
{
  long jm,jhi,inc;
  assert(xx[2]>xx[1],"must be in ascending order");
  if (x <= xx[1]){
    if (x == xx[1])
      *jlo=1;
    else
      *jlo=0;
    return;
  } else if(x >= xx[n]) {
    if(x > xx[n])
      *jlo=n+1;
    else /* case x==xx[n]: give lower jlo;  */
      *jlo=n-1;
    return;
  } 
  if (*jlo <= 0 || *jlo > n) {
    *jlo=0;
    jhi=n+1;
  } else {
    inc=1;
    if (x >= xx[*jlo]) {
      if (*jlo == n) return;
      jhi=(*jlo)+1;
      while (x >= xx[jhi]) {
	*jlo=jhi;
	inc += inc;
	jhi=(*jlo)+inc;
	if (jhi > n) {
	  jhi=n+1;
	  break;
	}
      }
    } else {
      if (*jlo == 1) {
	*jlo=0;
	return;
      }
      jhi=(*jlo)--;
      while (x < xx[*jlo]) {
	jhi=(*jlo);
	inc <<= 1;
	if (inc >= jhi) {
	  *jlo=0;
	  break;
	}
	else *jlo=jhi-inc;
      }
    }
  }
  while (jhi-(*jlo) != 1) {
    jm=(jhi+(*jlo)) >> 1;
    if (x > xx[jm])
      *jlo=jm;
    else
      jhi=jm;
  }
}

