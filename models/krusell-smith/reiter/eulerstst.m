% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% computes euler residuals in steady state;
% for usage, cf. findkstst.m
function res = eulerstst(Coh,Vcoh,R,w,TransProb)
  global Params;

  res = zeros(length(Coh),Params.nInd);

  c1 = cons_at(Coh,Vcoh);
  muc1 = margutil(c1);


  for ii=1:Params.nInd
    a1 = Coh(:,ii) - c1(:,ii);
    isConstrained = a1(:)<1e-6;

    muc2 = 0;
    Coh2 = repmat(R*a1,1,2) + repmat(w*Params.YIndStst(:)',size(a1,1),1);
    c2 = cons_at(Coh2,Vcoh);
    for iinext=1:Params.nInd
      prob = TransProb(ii,iinext);
      muc2 = muc2 + prob*margutil(c2(:,iinext));
    end

    c1bar = min(Coh(:,ii),invmargutil(Params.betta*R*muc2));
    res(:,ii) = (c1bar-c1(:,ii))./c1bar;
  end
  