% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% solves model for one parameter constellation
function doone;
  global Params;
  setparam;
  tag = maketag(Params);
  load shocks;
  rf = resfilename('V',tag,Params.iter);
  if(~fileexists(rf))
    if(Params.iter>0) %take results from last round
      rf = resfilename('V',tag,Params.iter-1);
    elseif(strcmp(Params.approxtype,'lin')==0)
      tag2 = maketag(Params,0,'lin');
      rf = resfilename('V',tag2,0);
    end
  end
  if(fileexists(rf))
    load(rf);
    recu(VAac,xEqu);
  else
    recu;
  end;

