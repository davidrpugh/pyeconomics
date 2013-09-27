% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% makes 3 adjustments to a distribution D
%  - transition of different types, with transition matrix PiTypes
%  - add Add to the grid
%  - multiply grid points by Fac
function Dnext = d_mixtypes(D,PiTypes,Add,Fac)

  assert(all(abs(sum(PiTypes')-1)<1e-10));
  [nx,nt] = size(D.p);

  if(nargin<4)
    Fac = ones(1,nt);
  end
  if(nargin<3)
    Add = zeros(1,nt);
  end

  Dnext.p = D.p * PiTypes;
  Dnext.x = D.x .* repmat(Fac,nx+1,1) + repmat(Add,nx+1,1);

