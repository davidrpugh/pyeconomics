function z = midpoint(x)
  if(isvector(x))
    n = length(x);
    z = 0.5*(x(1:n-1)+x(2:n));
  else
    n = size(x,1);
    z = 0.5*(x(1:n-1,:)+x(2:n,:));
  end
