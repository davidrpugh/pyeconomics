% Inputs:  k: k-th moment
%          nodes: (n,nTypes) matrix
function H = dst_h_moment(k,nodes)
  assert(k>=0,'order of moment must be nonnegative');
  width = diff(nodes);
  % ni = (nodes(:).^(k+1)) / (k+1);
  ni = (nodes.^(k+1)) / (k+1);
  H = diff(ni) ./ max(width,1e-300);  %to avoid warning "division by zero";
  izero = width==0;
  if(any(izero))
    % H(izero) = nodes([izero;false]) .^ k;
    nn = nodes(1:end-1,:);
    H(izero) = nn([izero]) .^ k;
  end;
  
  H = H(:);

