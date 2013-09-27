% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% solves for the steady state level of capital
global Params;

if(Params.BadStSt)  %the case_ required by DJJ
  %kstst = fzero(@kaggr,[1.0001 1.01]*Params.Knorm);
  %disp(kstst);
  kstst = 43;
else
  % solve for_ aggregate capital:
  warning off;
  kstst = fzero(@kaggr,[1.001 1.02]*Params.Knorm);
  warning on;
end

[resid,DKstst,DXstst,VCohStSt,VAacStSt] = kaggr(kstst); 

Coh = midpoint(DXstst.x);


if(Params.BadStSt)  % calculations required by DJJ
  GridK = [(0:0.01:4.999999)'; (5:0.1:100.0001)']; % DJJ
  GridX = nodesk2x(GridK,Params.Rstst,Params.wstst*Params.YIndStst);
  for i=1:Params.nInd
    xcrit(i) = VCohStSt{i}.x(1);
    kcrit(i) = (xcrit(i)-Params.wstst*Params.YIndStst(i))/Params.Rstst;  
  end
  GridC = cons_at(GridX,VCohStSt);
  Knext = GridX - GridC;
  FF = fopen(['res/cfuncstst' maketag(Params,1) '.txt'],'w');
  for i=1:size(Knext,1)
    fprintf(FF,'%12.2f  %14.8f  %14.8f  %14.8f  %14.8f\n',GridK(i),Knext(i,1),Knext(i,2),GridC(i,1),GridC(i,2));
  end
  fclose(FF);

  IE = load('ind_switch.txt');
  ne = length(IE);
  X = zeros(ne,1); 
  C = X;
  K = zeros(ne+1,1);
  K(1) = 43;
  for t=1:ne
    ie = IE(t);
    wagei = Params.YIndStst(ie)*Params.wstst;
    X(t) = Params.Rstst*K(t) + wagei;
    C(t) = cons_at(X(t),VCohStSt,ie);
    K(t+1) = X(t)-C(t);
  end
  FF = fopen(['res/simstst' maketag(Params,1) '.txt'],'w');
  for t=1:length(IE)
    fprintf(FF,'%12.0f  %14.8f  %14.8f\n',IE(t),K(t),C(t));
  end
  fclose(FF);
end

nodesKEuler = (0:0.01:100.0001)';
nodesXEuler = nodesk2x(nodesKEuler',Params.Rstst,Params.wstst*Params.YIndStst);
res = eulerstst(nodesXEuler,VCohStSt,Params.Rstst,Params.wstst,Params.TransProbStst);
name = 'stst';
if(Params.BadStSt)  %the case_ required by DJJ
  name = 'ststbad';
end
FF = fopen(['res/' name 'res' maketag(Params,1) '.txt'],'w');
assert(kcrit(2)<0);
fprintf(FF,'K at which constraint binding: %14.8f\n',kcrit(1));
ares = abs(res);
meanres = mean(ares(1:1000,:));  %only use the first 1000 points in K, as required by DJJ;
[maxres,imax] = max(ares);
fprintf(FF,'Euler residuals are in the form (c-\\bar c)/\\bar c (0.01 means 1 percent)\n');
fprintf(FF,'max EulerRes:         %14s  %14s\n',number_acc(maxres(1)),number_acc(maxres(2)));
fprintf(FF,'cap. at EulerRes:     %14.8f  %14.8f\n',nodesKEuler(imax(1)),nodesKEuler(imax(2)));
fprintf(FF,'average EulerRes:     %14s  %14s\n',number_acc(meanres(1)),number_acc(meanres(2)));
fprintf(FF,'Unempl = %14.8f; \\bar l =  %14.8f; \\tau = %14.8f\n',Params.Ustst,Params.lbar,Params.taustst);
fprintf(FF,'Wage = %14.8f; Interest rate(percent) =  %14.8f\n',Params.wstst,r_transform(Params.Rstst));
fclose(FF);
FF = fopen(['res/' name 'resdef' maketag(Params,1) '.txt'],'w');
fprintf(FF,'\\def\\maxeulerstst{%0s}\n',number_acc(max(ares(:)),1));
fprintf(FF,'\\def\\meaneulerstst{%0s}\n',number_acc(mean(mean(ares(1:1000,:))),1));
fclose(FF);

% WRITE COMPLETE SET OF EULER RESIDUALS:
FF = fopen(['res/' name 'euler' maketag(Params,1) '.txt'],'w');
for i=1:length(nodesKEuler)
  fprintf(FF,'%10.2f  %14s  %14s\n',nodesKEuler(i),number_acc(res(i,1)),number_acc(res(i,2)));
end 
fclose(FF);



assert(all(Params.DNodes==DKstst.x(:,1)));
% consumption, as a function_ of beginning-of-period assets:
Cstst = cons_at(nodesk2x(Params.DNodes,Params.Rstst,Params.wstst*Params.YIndStst),VCohStSt);
Caggrstst = sum(sum(DKstst.p.*midpoint(Cstst)));

Pstst = Params;
save(resfilename(name,maketag(Params,1)),'kstst','DKstst','DXstst','Pstst','Cstst','Caggrstst','VCohStSt','VAacStSt');
