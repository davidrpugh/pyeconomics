% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% loads initial distribution for simulations
function D = loaddistr(U)

  d = load('data/pdist.txt');
  D.x = repmat([d(1,1);d(:,1)],1,2);
  D.p = [d(:,2)*(U/sum(d(:,2))) d(:,3)*((1-U)/sum(d(:,3)))];

