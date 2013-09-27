function Dvec = selectdistr(m,D,Hvec,omega,Dstart);

  if(nargin<5)
    Dstart = D;
  end
  if(nargin<4 | isempty(omega))
    omega = ones(size(Hvec,1),1);
  end
  [Dvec,info] = solve4dens(D(:),Hvec,m,omega,Dstart(:));
  %if(info)
  %  error('distribution not found');
  %end
  kontr = Hvec'*Dvec - m;
  if(any(abs(kontr)>1e-6))
    disp(sprintf('warning: max error in selectdistr = %e',max(abs(kontr))));
  end;
  if(any(abs(kontr)>1e-2))
    disp(kontr);
    error('distribution not found');
  end;

  Dvec = reshape(Dvec,size(D));

