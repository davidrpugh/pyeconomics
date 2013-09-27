% min (f-g)' Omega (f-g)
% s.t. Hf=m
% 
% Omega is diagonal, the diagonal part is given by
%   the vector omega
function [f,info] = solve4dens(g,Ht,m,omega,f_start)

  if(nargin<5)
    f_start = g;
  end
  [n,p] = size(Ht);
  if(n ~= length(g) | p ~= length(m))
    error('wrong inputs in solve4dens');
  end;

  % Scaling:
  % typical element of f is 1 (cf. minbbfgs)
  % typical element of H is 1/n (since first equation sums up to 1):

  % s_h = n*sqrt(sum(Ht .* Ht)/n);
  % s_h = n*mean(abs(Ht));
  s_h = n*(abs(Ht)'*g);
  scale = reshape(1 ./ s_h,1,p);

  Ht = elelcol(Ht,scale);
  m = m.*scale';
  m = m * n; %extra scaling, to make f~1
  g = g*n; % same scaling for_ g;
  f_start = f_start*n; % same scaling for_ f_start;
  H = Ht';
  % typical element of omega is 1/n (s.t. opj.function is about 1):
  s_o = sum(omega.*omega);
  scale = sqrt(1 / n / s_o);
  omega = scale*omega;

  % Matlab's quadprog does not yet work reasonably for_ this
  % kind of problem; does not accept sparse matrices
  % f = quadprog(sparse(1:n,1:n,omega),-omega.*g,[],[],Ht',m,zeros(n,1));


  % starting value for f is g;
  imonitor = 0;
  myf = @myderivs;
  % initial penalty parameter;
  % don't set too high, creates convergence problems:
  c = 1e2;
  mu = zeros(p,1);
  xlow = zeros(n,1);
  % ff = myderivs(g, g,H,m,omega,c,mu,xlow);
  % [f,info,slack] = minbounds(g, xlow, myf, imonitor, g,H,m,omega,c,mu,xlow);
  if(0)
    [f,info,slack] = minbounds_sp(f_start, xlow, @myderivs2, imonitor, g,H,m,omega,c,mu,xlow);
    crit = 1e-8;
    for i=1:20
      if(info)
	c = c/100;
	% disp('Warning: convergence problems in solve4dens; reduce penalty parameter c');
      end
      discrep = H*f-m;
      if(i>2 & all(abs(discrep)<crit))
	% disp('solve4 success');
	break;
      end
      mu = mu + c*(discrep);
      c = min(10*c,1e6);
      [f,info,slack] = minbounds_sp(f, xlow, @myderivs2, imonitor, g,H,m,omega,c,mu,xlow);
    end;
    if(any(abs(discrep)>crit))
      disp('Warning: distribution not found: discrep is')
      disp(discrep');
      info = 1;
    else
      info = 0;
    end
  end
  f = solve4dsimple(H',g,m); info=0;
  f = f/n;  %undo extra scaling;


function [arg1,hess] = myderivs(f,g,H,m,omega,c,mu,xlow)
  % 0.5*(f-g)'Omega*(f-g) + 0.5*c*(H*f-m)'(H*f-m) + mu'*(H*f-m)

  fg = f-g;
  Ht = H';
  HH = Ht*H;
  resid = H*f - m;

  if(nargout==1) %arg1 is function value:
    arg1 = 0.5*(sum(fg.*fg.*omega) + c*(resid'*resid)) + mu'*resid;
  else
    arg1 = omega.*fg + c*(Ht*resid) + Ht*mu;
    hess = diag(omega) + c*HH;
  end;

function [arg1,I,hess,eps] = myderivs2(f,g,H,m,omega,c,mu,xlow)
  % 0.5*(f-g)'Omega*(f-g) + 0.5*c*(H*f-m)'(H*f-m) + mu'*(H*f-m)

  eps = 1e-7; %typical x-values should be around 1!;

  fg = f-g;
  resid = H*f - m;

  if(nargout==1) %arg1 is function value:
    arg1 = 0.5*(sum(fg.*fg.*omega) + c*(resid'*resid)) + mu'*resid;
  else    % arg1 is gradient!
    Ht = H';
    arg1 = omega.*fg + c*(Ht*resid) + Ht*mu;
    I = f<xlow+eps & arg1>0;
    notI = ~I;
    hess = zeros(size(f));
    HI = H(:,I);
    HnotI = H(:,notI);
    hess(notI) = - sherman_diag(arg1(notI),HnotI',c*HnotI,omega(notI));
    if(any(I))
      hess(I) = - arg1(I) ./ (omega(I) + c*diag(HI'*HI));
    end;
  end;


