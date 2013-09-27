% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% interpolates value function in the aggregate moments
% mom: each point is a row; (for 1 point, row vector)
% We assume we interpolate Vaac
function Vint = interp_mom(mom,ia,V)
  global Params;
  nn = size(V.v);

  assert(isvector(V.x));  %we interpolate VAac;
  Vint.x = V.x;
  Vint.v = 0;
  Vint.dv = 0;
  mom = mom(:)';
  if(strcmp(V.approxmom.type,'lin'))
    [I,W] = interp_simpl(inbounds(mom),V.approxmom.nodesmat,V.approxmom.n,0);
    for i=1:length(I)
      Vint.v = Vint.v + squeeze(W(i)*V.v(I(i),:,ia,:));
      Vint.dv = Vint.dv + squeeze(W(i)*V.dv(I(i),:,ia,:));
    end;
  elseif(strcmp(V.approxmom.type,'pwcub'))
    [I,W] = interp_simpl(mom,V.approxmom.nodesmat,V.approxmom.n,1);
    for i=1:length(I)
      Vint.v = Vint.v + squeeze(W(i)*V.v(I(i),:,ia,:));
      Vint.dv = Vint.dv + squeeze(W(i)*V.dv(I(i),:,ia,:));
    end;
  elseif(strcmp(V.approxmom.type,'spli'))
    B = splinebasis(V.approxmom.breaks,3,mom);
    nn = size(Vint.v);
    Vint.v = reshape(B*V.coeff_v,nn);
    Vint.dv = reshape(B*V.coeff_dv,nn);
  elseif(strcmp(V.approxmom.type,'complpoly'))
    global Params;
    B = polbasis(mom,Params.Bcoeff);
    nn = size(Vint.v);
    nc = size(V.coeff_v);
    Vint.v  = reshape(B*reshape(V.coeff_v(:,:,ia,:),nc(1),nc(2)*nc(4)),nc([2 4]));
    Vint.dv = reshape(B*reshape(V.coeff_dv(:,:,ia,:),nc(1),nc(2)*nc(4)),nc([2 4]));
  else
    error('not yet impl');
  end;

