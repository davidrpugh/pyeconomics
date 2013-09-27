% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% computes euler residuals
% for usage, cf. simul.m
function res = eulerres(Coh,Vcoh,VcohNext,R,w,ia,TransProb)
  global Params;

  res = zeros(length(Coh),Params.nInd);

  c1 = cons_at(Coh,Vcoh);
  muc1 = margutil(c1);


  for ii=1:Params.nInd
    j = (ia-1)*2 + ii;
    a1 = Coh(:,ii) - c1(:,ii);
    isConstrained = a1(:)<1e-6;

    muc2 = 0;
    sumprob = 0;
    for ianext=1:2
      Coh2 = repmat(R(ianext)*a1,1,2) + repmat(w(ianext)*Params.YInd{ianext}(:)',size(a1,1),1);
      c2 = cons_at(Coh2,VcohNext{ianext});
      for iinext=1:Params.nInd
	jnext = (ianext-1)*2 + iinext;
	prob = TransProb(j,jnext);
	sumprob = sumprob+prob;
	muc2 = muc2 + prob*R(ianext)*margutil(c2(:,iinext));
      end
    end
    assert(abs(sumprob-1)<1e-8);

    c1bar = min(Coh(:,ii),invmargutil(Params.betta*muc2));
    res(:,ii) = (c1bar-c1(:,ii))./c1bar;
  end
  