% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% generates filename for results file, depending on parameters
function fname = resfilename(type,tag,iter)
  if(strcmp(type,'stst') | strcmp(type,'ststbad')) %steady state
    fname = sprintf('res/%s_%s.mat',tag,type);
  else
    if(iter==-1)
      iter = 100;
    end
    fname = sprintf('res/%s%s_%d.mat',type,tag,iter);
  end
