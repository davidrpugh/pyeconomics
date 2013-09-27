% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
function Dsel = dsf_d(D,ia)
  h = make_hvec(D.x);
  m = (h'*D.p(:));
  s = mom2stat(m);
  s(3:end) = inbounds(s(3:end));
  Dsel = dsf_s(s,ia);




  