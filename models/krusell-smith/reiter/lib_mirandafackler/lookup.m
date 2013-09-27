% LOOKUP  Performs a table lookup.
% USAGE
%   ind=lookup(tabvals,x,endadj);
% INPUTS
%   tabvals : a sorted vector of n values
%   x       : an array of values
%   endadj  : a optional endpoint adjustment: 0, 1, 2 or 3.
% OUTPUT
%   ind : an array (same size as x) of indices from 1 to n
%
% Returns an array of size(x) with element (i,j) equal to
%   max k: x(i,j)>=tabvals(k)
%
% Optional endpoint adjustments:
%   0: no adjustments
%   1: values of x < min(tabvals) will return 
%        length(tabvals=tabvals(1))
%   2: values of x > max(tabvals) will return 
%        m-length(tabvals=tabvals(end))
%   3: adjustments 1 and 2 will be performed
%
% With endadj=3 all the indices are between 1 and n-1
% To find the nearest table value to each x use:
%   ind = lookup(tabvals,x,3);
%   ind = ind + (x-tabvals(ind) > tabvals(ind+1)-x);
%   nearest = tabvals(ind);
%
% Coded in C.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function ind=lookup(tabvals,x,endadj);

if nargin<2
  error('At least two parameters must be specified');
end
if nargin<3 endadj=0; end
if isempty(endadj) endad=0; end

n=prod(size(x));
if min(size(tabvals))>1
  error('tabvals must be a vector');
else 
  tabvals=tabvals(:);
  if any(diff(tabvals)<0)
    error('tabvals must be sorted in ascending order')
  end
end
m=length(tabvals);
if endadj>=2, m=m-length(find(tabvals==tabvals(end))); end

[temp,ind]=sort([tabvals(1:m); x(:)]);
temp=find(ind>m);
j=ind(temp)-m;
ind=reshape(temp-(1:n)',size(x));
ind(j)=ind(:);

if endadj==1 | endadj==3
  ind(ind==0)=length(find(tabvals==tabvals(1))); 
end
