% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% computes proxy distribution for moments s and exogenous state ia
function D = dsf_s(s,ia)
  global Params;
  m = stat2mom(s);
  % REFERENCE DISTRIBUTION:
  D = finddref(s(Params.iMomK:end),ia);
  h = make_hvec(D.x);
  % NOTHING TO BE DONE IF MOMENT CONDITIONS SATISFIED:
  if( max(abs(h'*D.p(:)-m))<1e-8)
    return;
  end
  % OTHERWISE_ CHOOSE D CLOSEST TO REFERENCE DISTRIBUTION:
  D.p = selectdistr(m,D.p,h,[],D.p);
  % CHECK MOMENT CONDITIONS:
  dd = max(abs(h'*D.p(:)-m));
  if(dd>1e-6)
    error('DSF not found');
  end


