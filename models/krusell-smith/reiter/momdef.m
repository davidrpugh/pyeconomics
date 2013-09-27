% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% generates structure with the information on the interpolation of moments
function f = momdef(nodes,type)
  if(iscell(nodes))
    [f.nodesmat,f.n] = gridmat(nodes);
    f.d = length(nodes);
  else
    if(isvector(nodes))
      f.n = length(nodes);
      f.nodesmat = nodes(:);
      f.d = 1;
    else
      [f.n,f.d] = size(nodes);
      f.nodes = nodes;
    end
  end;
  f.type = type;
  if(strcmp(type,'spli'))
    n = f.n;
    if(f.d==1)
      f.breaks = linspace(f.nodesmat(1,1),f.nodesmat(n(1),1),n(1)-2)';
    else
      f.breaks = cell(1,f.d);
      for i=1:f.d
	f.breaks{i} = linspace(f.nodesmat(1,i),f.nodesmat(n(i),i),n(i)-2)';
      end;
    end
    f.B = splinebasis(f.breaks,3,nodes);
  end;
  if(strcmp(type,'complpoly'))
    global Params;
    f.B = polbasis(nodes,Params.Bcoeff);
  end;

