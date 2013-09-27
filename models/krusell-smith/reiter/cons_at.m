% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% computes consumption at cash on hand x, given value function grid Vcoh, for type ii
function c = cons_at(x,Vcoh,ii)

  if(nargin==2)  %do it for_ both individual states
    nInd = length(Vcoh);  %Vcoh is a cell array
    c = zeros(size(x,1),nInd);
    for ii=1:nInd
      [v,dv,c(:,ii)] = interp_grid1d(x(:,ii),Vcoh{ii}.x,Vcoh{ii}.v,Vcoh{ii}.dv);
    end
  else
    c = zeros(size(x,1),1);
    [v,dv,c] = interp_grid1d(x,Vcoh{ii}.x,Vcoh{ii}.v,Vcoh{ii}.dv);
  end
