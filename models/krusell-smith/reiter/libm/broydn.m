% Broydn's method to solve system of nonlinear equations
% INPUT:
%     fname:  name of vector function that should be set to zero
%     xold:   starting value for parameters
%     opts:   a vector of options; this can be empty, then just default values are taken;
%             otherwise:
%       opts(1): tolerance level on function value (should be around 10^(-6))
%       opts(2): 1 if Jacobian through automatic differentiation 
%               (class deriv1) should be used, zero otherwise
%       opts(3): 1 if output on iterations is desired
%       opts(4): steplength for forward difference Jacobian, if used (should be around 10^(-5))
%     additional arguments will be passed to the function
%  OUTPUT:
%     x:      parameter vector that solves the system
%     check:  0 if ok, 1 if there is some problem (no solution found)
%
% The function fname should set the first element of its return vector to
% 1e100 in order to indicate overflow.
function [x, check, B] = broydn(fname,xold,opts,varargin);
  global broydn_init_B;
  % global Parameter broydn_init_B does not interfere with recursive calling of broydn,
  %  since only used in the beginning to initialize B, must be set anew before
  %  each call of broydn

  if(length(opts)>0)
    TOLF = opts(1);
  else  %default
    TOLF = 1e-6;
  end
  if(length(opts)>1)
    iAD = opts(2);
  else  %default
    iAD = 0;
  end
  if(length(opts)>2)
    iprint = opts(3);
  else  %default
    iprint = 0;
  end
  if(length(opts)>3)
    step_jacob = opts(4);
  else  %default
    step_jacob = 1e-5;
  end
    
  % initialize outputs:
  check = 0;
  B = broydn_init_B;
  broydn_init_B = [];  %must be set each time before calling broydn!

  xold = xold(:);
  x = xold;

  MAXITS = 2000;
  EPS = 1.0e-7;
  TOLX = EPS;
  STPMX = 100.0;
  TOLMIN = 1.0e-6;
  n=size(x,1);
  
  [f,fvec]=fmin_br(x,fname,varargin{:});
  fvcold = fvec;
  fold=f;
  fx = fvec;
  if abs(fvec(1))>=1e100
    error('overflow in function given to broydn at initial vector');
  end;
  test = max(abs(fvec));
  if (test<TOLF)  %changed from 0.01*TOLF to TOLF;
    return;
  end;
  stpmax=STPMX*max(sqrt(x'*x),n);
  if(isempty(B))
    restrt=1;
  else
    restrt=0;
  end
  since_restrt = 0;
  alam = 1;
  alam2 = 0;
  for its=1:MAXITS
    since_restrt = since_restrt+1;
    txt = sprintf('its:  %d; f:  %e; step size taken: %e, %e\n',its,getval(f),alam);
    if(iprint==1)
      disp(txt(1:end-1));  %without EOL
    elseif(iprint==2)
      appendf('broydn_out.txt',txt);
    end;
    % determine whether restart (new computation of Jacobian):
    if(fileexists('broydn_restrt.txt'))
      system('del broydn_restrt.txt');
      restrt=1;
    end
    if(restrt==1 | since_restrt>=3*n) %compute Jacobian at 3n iterations since last restart
      since_restrt = 0;
      alam = 1;  %starting value of backstepping parameter for line search;
      if(iAD==1)
	xder = deriv1(x);
	xjac = feval(fname,xder,varargin{:});
	B = full(getjacob(xjac));
	fvec = getval(fvec);
      else
	% broydn_in_jacob = 1;
	B = jacob(fname,{x,fx},step_jacob,varargin{:});
	assert(size(B,1)==size(B,2));
	% broydn_in_jacob = 0;
      end;
%	disp(sprintf('condition number in broydn is %e',cond(B)));
    elseif(its>1)
      s = x - xold;
      skip=1;
      w = (fvec-fvcold) - B*s;
      for i=1:n
	if abs(w(i)) >= EPS*(abs(fvec(i))+abs(fvcold(i)))
	  skip=0;
	else
	  w(i)=0.0;
	end;
      end;
      if skip==0
	B = B + w*s' / (s'*s);
      end;
    end;
%    if(cond(B)>1e15)
%      error('condition B')
%    end
    maxr = max(abs(B'))';
    % p = - B \ fvec;
    p = - (B./(repmat(maxr,1,size(B,1)))) \ (fvec./maxr);
    g = B'*fvec;
    xold = x;
    fvcold = fvec;
    fold=f;
    [x,f,fvec,check,alam2] = lnsrch_br(xold,fold,g,p,stpmax,fname,varargin{:});
    if(alam2==1)
      alam = 5*alam;
    else
      alam = alam*alam2*2;
    end
    alam = max(min(alam,1),1e-5);
    fx = fvec;
    test = max(abs(fvec));
    if (test < TOLF)
      check=0;
      return;
    end;
    if check==1
      if (restrt)
	return;
      else 
	test=0.0;
	den=max(f,0.5*n);
	test = max(abs(g) .* max([abs(x');ones(1,n)])') / den;
	if (test < TOLMIN)
	  return;
	else
	  restrt=1;
	end;
      end;
    else 
      restrt=0;
      test= max(abs(x-xold) ./ max([abs(x');ones(1,n)])');
      if (test < TOLX)
	return;
      end;
    end;
  end;
  check = 1;
  % warning('MXITS exceeded in broydn');
  return;
  
 

function [x,f,fvec,check,alam] = lnsrch_br(xold,fold,g,p,stpmax,varargin);
  
  n = size(xold,1);
  ALF = 1.0e-4;
  TOLX = 1.0e-7;
  check=0;
  sum = p'*p;
  sum=sqrt(sum);
  if (sum > stpmax)
    p = p * stpmax/sum;
  end;
  slope = g'*p;
  test=max(abs(p) ./ max([abs(xold');ones(1,n)])');
  % min() introduced by MR 5.4.2006:
  % otherwise, alamin can be >1 !!
  alamin=min(0.1,TOLX/test);  
  alam=1.0;
  for its=1:100000
    x=xold + alam*p;
    try
      [f,fvec]=fmin_br(x,varargin{:});
    catch
      f = 1e100;  %something outrageous;
    end;
    if (alam < alamin) 
      x = xold;
      check=1;
      return;
    else
      if f <= fold+ALF*alam*slope
	return;
      else 
	if alam == 1.0
	  tmplam = -slope/(2.0*(f-fold-slope));
	else 
	  rhs1 = f-fold-alam*slope;
	  rhs2=f2-fold2-alam2*slope;
	  a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	  b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	  if (a == 0.0)
	    tmplam = -slope/(2.0*b);
	  else 
	    disc=b*b-3.0*a*slope;
	    if (disc<0.0)
	      error('Roundoff problem in lnsrch');
	    else
	      tmplam=(-b+sqrt(disc))/(3.0*a);
	    end;
	  end;
	  if (tmplam>0.5*alam)
	    tmplam=0.5*alam;
	  end;
	end;
      end;
    end;
    alam2=alam;
    f2 = f;
    fold2=fold;
    alam=max(tmplam,0.1*alam);
  end;
  



function [f,fvec] = fmin_br(x,funcname,varargin);
  fvec = feval(funcname,x,varargin{:});
  f = 0.5*fvec'*fvec;
