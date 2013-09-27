% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% simulates the model
function [ser,id] = dosimul(tagRes,tagSim,tSimul,iChange)
  global Params;

  if(nargin<4)
    iChange=0;
  end

  setparam;
  loadstst;

  load(resfilename('V',tagRes,Params.iter));

  if(Params.whichShocks==1)
    IA = load('agg_switch.txt');
  elseif(Params.whichShocks==2)
    IA = kron([1;2],ones(100,1));
  elseif(Params.whichShocks==3)
    IA = [2*ones(200,1);ones(200,1)];
  else  % random:
    rand('seed',424872);
    shockalld = rand(10000,2);
    T = length(shockalld);
    IA = zeros(T,1);
    if(shockalld(1)>0.5)
      IA(1) = 2;
    else
      IA(1) = 1;
    end
    for t=2:T
      IA(t) = changestaterandomly(IA(t-1),shockalld(t,1), Params.ProbAggr);
    end;
  end
  IE = load('ind_switch.txt');



  D = startdistr(IA(1));

  [ser,id] = simul(D,IA(1:tSimul),VAac,xEqu,@recu_transmom,IE);

  jedcstat(ser,id,tagSim);


  if(~isempty(tagSim))
    cmd = sprintf('ser%s=ser;id%s=id; save res/ser%s ser%s id%s;',...
		  tagSim,tagSim,tagSim,tagSim,tagSim);
    eval(cmd);
  end