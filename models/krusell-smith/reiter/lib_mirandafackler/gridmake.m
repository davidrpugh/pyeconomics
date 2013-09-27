% GRIDMAKE Forms grid points
% USAGE
%   X=gridmake(x);
%   X=gridmake(x1,x2,x3,...);
%   [X1,X2,...]=gridmake(x1,x2,x3,...);
%   X=gridmake({y11,y12},x2,{y21,y22,y23});
%
% Expands matrices into the associated grid points.
% If N is the dx2 matrix that indexes the size of the inputs
%   GRIDMAKE returns a prod(N(:,1)) by sum(N(:,2)) matrix.
%   The output can also be returned as either
%      d matrices or 
%      sum(N(:,2)) matices
% If any of the inputs are grids, they are expanded internally
% Thus
%    X=gridmake({x1,x2,x3})
%    X=gridmake(x1,x2,x3)
% and
%    x={x1,x2,x3}; X=gridmake(x{:})
% all produce the same output.
%
% Note: the grid is expanded so the first variable change most quickly.
%
% Example:
%  X=gridmake([1;2;3],[4;5])
% produces
%     1     4
%     2     4
%     3     4
%     1     5
%     2     5
%     3     5
% 
% The function performs an action similar to NDGRID, the main difference is in the
%   increased flexability in specifying the form of the inputs and outputs.
%
% Also the inputs need not be vectors.
%   X=gridmake([1;2;3],[4 6;5 7])
% produces
%     1     4     6
%     2     4     6
%     3     4     6
%     1     5     7
%     2     5     7
%     3     5     7
% 
% See also: ndgrid

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function varargout=gridmake(varargin)
m=prod(size(varargin));
n=nargout;
Z=[];
d=zeros(1,m+1);
for i=1:m
  if isa(varargin{i},'cell')
    Z=gridmake2(Z,gridmake(varargin{i}{:}));
  else
    Z=gridmake2(Z,varargin{i});
  end
  d(i+1)=size(Z,2);
end

varargout=cell(1,max(n,1));
if n<=1
  varargout{1}=Z;
elseif n==m
  for i=1:m
    varargout{i}=Z(:,d(i)+1:d(i+1));
  end
elseif n==size(Z,2)
  for i=1:n
    varargout{i}=Z(:,i);
  end
else
  error(['An improper number of outputs requested - should be 1, ' num2str(m)  ' or ' num2str(size(Z,2))])
end

% Expands gridpoints for 2 matrices
function Z=gridmake2(X1,X2)
if isempty(X1); Z=X2; return; end
m=size(X1,1);
n=size(X2,1);
ind1=(1:m)';
ind2=1:n;
Z=[X1(ind1(:,ones(n,1)),:) X2(ind2(ones(m,1),:),:)];