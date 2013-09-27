% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% Coh is cash on hand, before consumption
% Aac are assets after consumption
function [Coh,v,dv,c] = v_optimc(Aac,va,dva)
  
    c = invmargutil(dva(:));
    Coh = Aac+c; 
    v = util(c)+va(:);
    dv = dva;

    % same as, but shorter than
    % dv = margutil(c);

  