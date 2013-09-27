% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% defines equilibrium:
%  equilibrium found if xequ2 = xequ;
% Inputs:
%   xequ: vector of equilibrium variables
%   DCoh:    distribution of cash-on-hand, after income accrues
%   ia:   current aggregate productivity state (ia=1,2)
%   Vaac: value function, as function of end-of-period assets and aggregate states
% Output:
%  xequ2: equilibrium found if xequ2 = xequ;
function xequ2 = reaction(xequ,DCoh,ia, Vaac)
  global Params;

  % interp_mom needs row vector:
  VInt = interp_mom(inbounds(xequ(:))',ia,Vaac);
  %[isc,notc] = v_concave(Params.VcohAtNext);
  %assert(isc==1);
  
  for ii=1:Params.nInd
    [Params.VcohEndog{ii}.x,Params.VcohEndog{ii}.v,...
     Params.VcohEndog{ii}.dv,Params.VcohEndog{ii}.contr] = ...
	v_optimc(VInt.x,VInt.v(:,ii),VInt.dv(:,ii));
  end;


  Params.DaacTry = transmatx2a(DCoh.x,repmat(Params.DNodes,1,Params.nInd),VInt,DCoh);
  % Params.DaacTry = advancex2a(DCoh,Params.VcohEndog);
  h = make_hvec(Params.DaacTry.x);
  mom2 = mom2stat(h'*Params.DaacTry.p(:)); 
  xequ2 = mom2(Params.iMomK:end);

