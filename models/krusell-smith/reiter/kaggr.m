% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% residual file to find aggregate k in steady state
% input: guess k
% outputs:
%   z: residual, equals zero in equilibrium
%   Daac: invariant distribution end-of-period capital
%   DXstst: invariant distribution cash on hand
%   V:  steady state value function
%   VAac: steady state value function, evaluated end of period (after consumption)
function [z,Daac,DXstst,V,VAac] = kaggr(k)
  global Params;

  [V,VAac] = dynpro(k);
  nt = Params.nInd;
  nD = Params.nDistr;
  
  nodesK = Params.DNodes;
  nodesK2 = repmat(nodesK,1,nt);
  nodesX = zeros(length(nodesK),2);
  for i=1:nt
    nodesX(:,i) = Params.Rstst*nodesK + Params.YIndStst(i)*Params.wstst;
  end

  Ta2x = kron(Params.TransProbStst',speye(nD));

  %xnext = zeros(nD+1,nt);
  %for i=1:nt
  %  xnext(:,i) = coh2a(nodesX(:,i),V{i});
  %end
  %loc = lookup(nodesK,xnext,3);
  %[iTo,iFr,Val] = forwmat_h(nodesK2,xnext,loc,diff(nodesK2));
  %Tx2a = sparse(iFr,iTo,Val,nt*nD,nt*nD);


  Tx2a = transmatx2a(nodesX,nodesK2,VAac);


  K2K = Tx2a * Ta2x;

  opts.disp=0;
  [x,eval] = eigs(K2K,[],1,1+1e-10,opts);
  D = x/sum(x);
  assert(min(D)>-1e-10);
  D = max(D,0);
  D = D/sum(D);

    
  Daac.x = repmat(Params.DNodes,1,nt);
  Daac.p = reshape(D,Params.nDistr,nt);
  meanK = dst_moment(1,Daac);
  z = meanK - k;




  if(nargout>1)
    DXstst = advancea2x(Daac,Params.TransProbStst,Params.Rstst,Params.YIndStst'*Params.wstst);
  end;
  disp(sprintf('%g  ',[k meanK Params.Rnorm-Params.Rstst z]));
