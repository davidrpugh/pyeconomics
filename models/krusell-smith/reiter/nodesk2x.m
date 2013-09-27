% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% generates grid of cash-on-hand to given grid of capital, using
%   interest rate R and wage wagei
function nodesX = nodesk2x(nodesK,R,wagei)
  global Params;

  nodesX = zeros(length(nodesK),Params.nInd);
  for i=1:Params.nInd
    nodesX(:,i) = R*nodesK + wagei(i);
  end
