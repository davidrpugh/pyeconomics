% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% create vector H for moment computation:
function [hvec] = make_hvec(nodes,nMom)
  global Params;

  if(nargin<2)
    nMom = Params.nMoments;
  end

  h1a = dst_h_moment(1,nodes(:,1));
  h1b = dst_h_moment(1,nodes(:,2));
  h0 = ones(size(h1a));
  zz = zeros(size(h1a));
  hvec = [ [h0;h0] [h0;zz] [h1a;h1b] ];
  Params.iMomK = 3;
  if(nMom==2)
    if(Params.useC)
      c1 = interp1lin(Params.DNodes,Params.Cstst(:,1),midpoint(nodes(:,1)));
      c2 = interp1lin(Params.DNodes,Params.Cstst(:,2),midpoint(nodes(:,2)));
      c = [c1;c2];
      % hvec = [hvec  (c-linproject(c,hvec))];
      hvec = [hvec  c];
    else
      % variance:
      hvec = [hvec  [h1a;zz] ];
    end
  end;

