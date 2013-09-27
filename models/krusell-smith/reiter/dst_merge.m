% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% merges columns of D into one marginal distribution
function D2 = dst_merge(D)
  m = size(D.p,2);
  X = D.x(:,1);
  for j=2:m
    X = [X;D.x(:,j)];
  end
  D2.x = unique(sort(X));
  F = dst_distrfunc(D,D2.x);
  D2.p = sum(diff(F),2);

  