% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% adjust unemployment rate in distribution D
function D = adjustd(D,ia,iAssert);
  global Params;
  U = Params.U(ia);
  if(nargin>2 & iAssert)
    assert(abs(U-sum(D(:,1)))<1e-5);
  end
  D(:,1) = D(:,1) * (U/sum(D(:,1)));
  D(:,2) = D(:,2) * ((1-U)/sum(D(:,2)));
