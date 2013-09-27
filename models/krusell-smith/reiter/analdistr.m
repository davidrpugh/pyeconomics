% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% loads the distribution from simulation and extracts information
% for reference distribution
load(Params.RefDistrib,'allD');

% RELEVANT RANGE OF DISTRIBUTION:
loadstst;
relevant = sum(DKstst.p>1e-8,2);
last_relevant = find(relevant,1,'last');


clear deltaP X indx;

T = length(allD);
indx{1} = false(T,1);
indx{2} = false(T,1);
K = zeros(T,1);
C = zeros(T,1);
P = zeros([size(allD{1}.p) T]);
h = make_hvec(repmat(Params.DNodes,1,2),2);
for t=1:T
  for j=1:2
    x = allD{t}.x(:,j);
    if(abs(x(end)-Params.DNodes(end))>1e-8)
      assert(0);
      p(:,j) = changenodesneu(x,Params.DNodes,allD{t}.p(:,j));
    else
      p(:,j) = allD{t}.p(:,j);
    end
  end
  m = h'*p(:);
  if(m(2)>0.07)
    IA(t) = 1;
    indx{1}(t)=true;
  else
    IA(t) = 2;
    indx{2}(t)=true;
  end
  P(:,:,t) = p;
  K(t) = m(3);
  C(t) = m(4);
end

  
for j=1:2
  ii{j} = find(indx{j});
  y = C(ii{j});
  X = [ones(length(ii{j}),1) K(ii{j}) K(ii{j}).^2];
  % X = [ones(length(ii{j}),1) K(ii{j})];
  CofK{j} = (X'*X)\(X'*y);
  yhat = X*CofK{j};
  e = y-yhat;
  % disp([min(e) max(e)]);
  % plot(X(:,2),y,'.',X(:,2),yhat,'.r');
end

Dmean{1} = allD{1};
Dmean{2} = allD{1};
for j=1:2
  Dmean{j}.p = squeeze(mean(P(:,:,IA==j),3));
end
save(Params.RefDistrib,'allD','CofK','Dmean');
