/*
  To make a matlab function from this file, type
  "mex value_c.c"
  on the command level or from matlab on tornasol.
  That's all. (It produces a file value_c.mexsol)
*/

/* LEAVE THIS LINE ALWAYS UNCHANGED:  */
#include "cfunc.h"

/* LEAVE THIS LINE ALWAYS UNCHANGED:  */

void schum_prep(double *c, double *dvc, int *ilinear, double *xa, double *xb, double *va, double *vb, double *dva, double *dvb,int n);
double schum_eval(double *deriv, double x0, double  c, double  dvc, 
		  double  xa, double  xb, double  va, double  vb, 
		  double  Dva, double  Dvb);
double PolyEval3PDD(double *deriv, double x, double  a, double  va, double  dva, 
		    double   b, double  dvb);
void Hunt(const double *xx, int n, double x, int *jlo);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  double *x0,*x,*v,*dv,x0adj,*vx,*dvx;
  double *xa,*xb,*va,*vb,*dva,*dvb;
  double *c, *dvc;
  int *ilinear;
  int i,n,nx,loc;

  DECL_MAT(m_x0);
  DECL_MAT(m_x);
  DECL_MAT(m_v);
  DECL_MAT(m_dv);
  /* DECLARE OUTPUT MATRICES: */
  DECL_OUTPUT(m_vx);
  DECL_OUTPUT(m_dvx);

  /* FUNCTION GETS 2 INPUTS AND 1 OUTPUT ARGUMENT:  */
  CHECK_ARGN(4,2);


  GET_ARG_MAT(m_x0,1);
  GET_ARG_MAT(m_x,2);
  GET_ARG_MAT(m_v,3);
  GET_ARG_MAT(m_dv,4);

  /*  printf("really me\n"); */
  n = m_x0_nr;
  nx = m_x_nr - 1;

  assert(nx+1==m_v_nr && nx+1==m_dv_nr,'dimensions dont fit');/*MORE CHECKS!! */
  
  CREATE_OUTPUT(1,m_vx,n,1,REAL);
  CREATE_OUTPUT(2,m_dvx,n,1,REAL);

  x0 = &MATEL(m_x0,1,1)-1;
  x = &MATEL(m_x,1,1)-1;
  v = &MATEL(m_v,1,1)-1;
  dv = &MATEL(m_dv,1,1)-1;
  vx = &MATEL(m_vx,1,1)-1;
  dvx = &MATEL(m_dvx,1,1)-1;

  c = VECTOR(double,nx);
  dvc = VECTOR(double,nx);
  ilinear= VECTOR(int,nx);

#if 0
  if(n<length(x)){
    
    xa = x;
    xb = x+1;
    va = v(loc);
    vb = v(loc+1);
    dva = dv(loc);
    dvb = dv(loc+1);  
    dvmean = (vb-va)./(xb-xa);

    [c,dvc,ilinear] = schum_prep(xa,xb,va,vb,dva,dvb);


    fac = (x0(ilinear)-xa(ilinear)) ./ (xb(ilinear)-xa(ilinear));
    vx(ilinear) = (1-fac).*va(ilinear) + fac.*vb(ilinear);
    dvx(ilinear) = dvmean(ilinear);

    inl = ~ilinear;
    [vx(inl),dvx(inl)] = schum_eval(x0(inl), c(inl), dvc(inl), xa(inl), xb(inl), va(inl), vb(inl), dva(inl), dvb(inl));
  }
#endif
  xa = x;
  xb = x+1;
  va = v;
  vb = v+1;
  dva = dv;
  dvb = dv+1;
  
  schum_prep(c,dvc,ilinear,xa,xb,va,vb,dva,dvb,nx);

  INCR(i,1,n){
    x0adj = x0[i];
    if(x0adj<x[1]){
      if(x0adj<x[1]-1e-12)
       mexErrMsgTxt("below bounds in schumaker");
      else
	x0adj = x[1];
    }
    if(x0adj>x[nx+1]){/*+1, since nx = length(x)-1; */
      if(x0adj > x[nx+1]+1e-12)
       mexErrMsgTxt("above bounds in schumaker");
      else
	x0adj = x[nx+1];
    }

    Hunt(x, nx+1, x0adj, &loc); /*loc <= length(x)-1; */
    if(ilinear[loc]){  /*all x0 that fall in interval with linear interpolation; */
      vx[i] = va[loc] + dvc[loc]*(x0adj-xa[loc]);
      dvx[i] = dvc[loc];
    } else {
      vx[i] = schum_eval(&dvx[i],x0adj, c[loc], dvc[loc], xa[loc], xb[loc], va[loc], vb[loc], dva[loc], dvb[loc]);
    }
  }
  free(c);
  free(dvc);
  // printf("%d  %d %d %d %d\n",ilinear[0],ilinear[1],ilinear[2],ilinear[3],ilinear[4]);
  free(ilinear);
}


double schum_eval(double *deriv, double x0, double  c, double  dvc, double  xa, double  xb, double  va, double  vb, 
		      double  Dva, double  Dvb)
{
  double v;
  if(x0<=c){  /* <= important, to handle case where c=b; */
    v  = PolyEval3PDD(deriv,x0,xa,va,Dva,c,dvc);
  } else {
    v = PolyEval3PDD(deriv,x0,xb,vb,Dvb,c,dvc);
  }
  return v;
}

double PolyEval3PDD(double *deriv, double x, double  a, double  va, double  dva, 
		    double   b, double  dvb)
{
  int i;
  double d2va,d2vad,d,v;
  d2va = (dvb-dva)/(b-a);
  d = x-a;
  d2vad = d2va*d;
  *deriv = dva + d2vad;
  v = va + d*(dva + 0.5*d2vad);
  return v;
}


void schum_prep(double *c, double *dvc, int *ilinear, double *xa, double *xb, double *va, double *vb, double *dva, double *dvb,int n)
{
  double nan,dvmean,crit,fac,dev,t,aux;
  int i,iexact,iequ,ilow,ihigh;
  double dum0=0,dum1=0;
  //  nan = 0./0.;
  nan = dum0/dum1;
  /* c = NaN*ones(n,1); */
  /* dvc = NaN*ones(n,1); */

  INCR(i,1,n){
    c[i] = dvc[i] = nan;
    dvmean = (vb[i]-va[i])/(xb[i]-xa[i]);

    /* exactly satisfies quadratic equation; */
    /* monotonicity and concavity satisfied, if in the data; */
    iexact = (dva[i]+dvb[i]) == 2*dvmean;
    if(iexact){
      c[i] = xb[i];
      dvc[i] = dvb[i];
      ilinear[i] = 0;
    } else {
      /* inotex = ~iexact; */
      iequ = (dva==dvb); /* & inotex; */
      if(iequ){
	c[i] = 0.5*(xa[i]+xb[i]);
      } else {
	c[i] = (vb[i] - dvb[i]*xb[i] - va[i] + dva[i]*xa[i]) / (dva[i] - dvb[i]);
      }
  
      crit = 1e-10*(xb[i]-xa[i]);
      ilinear[i] = (fabs(c[i]-xa[i])<crit || fabs(c[i]-xb[i])<crit); /* & inotex; */
      if(ilinear[i])
	printf("a=%f; b=%f;c=%f;\n",xa[i],xb[i],c);
      if(ilinear[i]){
	dvc[i] = dvmean; /*dvc set to mean slope if linear, can be used for interpol!! */
      } else {
	ilow = c[i]<xa[i]; /* & inl; */
	if(ilow){
	  dev = xa[i] - c[i];
	  fac = 1/(1+exp(-dev)) - 0.5;
	  c[i] = fac*xb[i] + (1-fac)*xa[i];
	}
	ihigh = c[i]>xb[i];  /*why missing inl in matlab file?? */
	if(ihigh){
	  dev = c[i]-xb[i];
	  fac = 1/(1+exp(-dev)) - 0.5;
	  c[i] = fac*xa[i] + (1-fac)*xb[i];
	}
	t = (c[i]-xa[i])/(xb[i]-xa[i]);
	aux = 2*dvmean - t*dva[i] - (1-t)*dvb[i];
	dvc[i] = aux;
      }
    }
  }
}

/* xx is unit offset!! */
void Hunt(const double *xx, int n, double x, int *jlo)
{
  int jm,jhi,inc;
  assert(xx[2]>xx[1],"must be ascending");
  if (x <= xx[1]){
    if (x == xx[1])
      *jlo=1;
    else
      *jlo=0;
    return;
  } else if(x >= xx[n]) {
    if(x > xx[n])
      *jlo=n+1;
    else /* case x==xx[n]: give lower jlo; */
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
