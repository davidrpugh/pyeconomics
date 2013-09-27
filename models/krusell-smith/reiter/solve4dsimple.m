% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
function f = solve4dsimple(H,g,m0)
  SolutionFound = 0;

  [n,p] = size(H);

  f = zeros(n,1);
  IsZero = false(n,1);
  m = m0;
    
  fracEliminate = 0.05;

  while(1)
    lauf = 0;
    nActive = n - sum(IsZero);
    if(nActive<p)
      break;
    end
    indx = find(IsZero==false);
      

    HAct = H(indx,:);
    gAct = g(indx);
    lam = (HAct'*HAct)\(m-HAct'*gAct);

    x = gAct + HAct*lam;
      
      
    nwrong = sum(x<0);
    if(nwrong==0)
      SolutionFound = 1;
      f(indx) = x;
      break;
    else 
      nelim = floor(fracEliminate*nActive) + 1;
      if(nwrong<=nelim)
	IsZero(indx(x<0)) = true;
      else
	[unused,iSmall] = sort(x);
	IsZero(indx(iSmall(1:nelim))) = true;
      end
    end
  end

  % fprintf(1,'nact:  %d\n',nActive);
  assert(SolutionFound);