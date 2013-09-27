% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% generates a default structure for a one-dimensional value function
function g = grid1d_def_interp(varargin);
  g = grid1d_def(varargin{:});
  g.interpAboveBound = @v_above;
  g.interpBelowBound = @v_below;
  g.interpCAboveBound = @c_above;
  g.interpCBelowBound = @c_below;

