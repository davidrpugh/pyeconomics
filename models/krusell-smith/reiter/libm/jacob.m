% computes Jacobian by forward differences
% Input:
% func: string, name of function
% x: point at which to take Jacobian
%    if function value at x is already known, then give x={starting point, function value}
% step: scalar, relative stepwidth;
% more arguments will be passed on to the function;
%
function jac = jacob(func,x,step,varargin);
  if(iscell(x))  %function value at x is given:
    f0 = x{2};
    x = x{1};
  else
    f0 = feval(func,x,varargin{:});
  end
  n = size(x,1);
  m = size(f0,1);
  jac = zeros(m,n);
  x0 = x;
  for i=1:n
    step2 = step*max(1,abs(x0(i)));
    x = x0;
    x(i) = x0(i) + step2;
    jac(1:m,i) = (feval(func,x,varargin{:}) - f0)/step2;
  end;
