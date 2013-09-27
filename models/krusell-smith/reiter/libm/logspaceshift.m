% grid of n points between xa and xb,
% equidistant in log(x+xshift)
% if nargin==5:
function [grid,xshift] = logspaceshift(xa,xb,n,x2,n_at_x2);
	       
  if(nargin==5)
    %   xshift such that approximately n_at_shift points from xa to xshift
    %   solves log(x2+xshift)-log(xa+xshift) = (n_at_x2)/n * (log(xb+xshift)-log(xa+xshift))
    frac = n_at_x2 / n;
    xshift = fzero(@find_xshift,[-xa+1e-8;1000*xb],[],xa,xb,x2,frac);
  elseif(nargin==4)
    xshift = x2;
  else  %(nargin<4)
    xshift = 0;
  end;
  grid = exp(linspace(log(xa+xshift),log(xb+xshift),n))-xshift;
  % to avoid roundoff errors at end points:
  grid(1) = xa;
  grid(end) = xb;


%   solves log(x2+xshift)-log(xa+xshift) = frac * (log(xb+xshift)-log(xa+xshift))
function z = find_xshift(xshift,xa,xb,x2,frac)
  z = log(x2+xshift)-log(xa+xshift) - frac * (log(xb+xshift)-log(xa+xshift));