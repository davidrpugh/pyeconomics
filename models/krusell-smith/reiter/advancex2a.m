% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% forward distribution of cash-on-hand to distribution end-of-period capital
function Daac = advancex2a(DCoh,Vcoh)
  global Params;

  a = zeros(size(DCoh.x));
  for ii=1:Params.nInd
    if(nargin>2)
      for ii=1:Params.nInd
	Vii = struct('x',Vi.x(:,ii),'v',Vi.v(:,ii),'dv',Vi.dv(:,ii));
	a(:,ii) = coh2a(DCoh.x(:,ii),Vii);
      end
    else
      a(:,ii) = coh2a(DCoh.x(:,ii),Vcoh{ii});
    end
  end;
  Daac.x = repmat(Params.DNodes,1,Params.nInd);
  Daac.p = forward_linear(Daac.x,a,DCoh.p);
