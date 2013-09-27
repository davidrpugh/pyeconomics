% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
function [ser,idSer] = simul(Dbeg,IA,VAac,xEqu,funceval,IEmpl)
  global Params;

  tic;
  'enter simul'

  T = length(IA);

  % skip tSkip first distributions:
  tSkip = 5000;
  if(Params.iSaveDistr & T>tSkip)
    allD = cell(T-tSkip,1);
  end
  


  % number of moments and predictions
  nMP = Params.nMoments*(1+Params.nPred);
  nMoreMoms = 5;
  nVars = 36;
  ncSer = nMP + 3*nMoreMoms + nVars + 2;  %2 from h2;
  idSer.mom = 1;
  idSer.moremom = nMP+1;
  nv = nMP + 3*nMoreMoms;
  idSer.m2 = nv+1;
  nv = nv + 2;
  idSer.ia = nv+1;
  idSer.meanres = nv+2;
  idSer.maxres = nv+3;
  idSer.dD = nv+4;
  idSer.dDc = nv+5;
  idSer.meanU = nv+6;
  idSer.meanE = nv+7;
  idSer.fracConstrU = nv+8;
  idSer.fracConstrE = nv+9;
  idSer.fracConstr = nv+10;
  idSer.C = nv+11;
  idSer.Y = nv+12;
  idSer.R = nv+13;
  idSer.w = nv+14;
  idSer.L = nv+15;
  idSer.cHH = nv+16;
  idSer.kHH = nv+17;
  idSer.yHH = nv+18;
  idSer.ie = nv+19;
  idSer.cstst = nv+20;
  idSer.Klag = nv+21;
  idSer.Inv = nv+22;
  idSer.capU = nv+23;
  idSer.capE = nv+24;
  offs = nv+24;
  idSer.p05 = offs+1;
  idSer.p10 = offs+2;
  idSer.p20 = offs+3;
  idSer.p30 = offs+4;
  offs = offs+4;
  idSer.p05U = offs+1;
  idSer.p10U = offs+2;
  idSer.p20U = offs+3;
  idSer.p30U = offs+4;
  offs = offs+4;
  idSer.p05E = offs+1;
  idSer.p10E = offs+2;
  idSer.p20E = offs+3;
  idSer.p30E = offs+4;
  
  ser = zeros(T,ncSer);
  ser(:,idSer.ia) = IA(1:T);

  % initial individual wealth
  kHHlag = 43;

	
  Params.eulerstats = [];
  meanres = 0;
  maxres = 0;


  for t=1:T
    if(mod(t,500)==0)
      fprintf(1,'%d  ',t);
    end
    ia = IA(t);

    if(t>1)
      % UPDATE DISTRIBUTION:
      Dbeg = d_mixtypes(Daac,Params.ProbInd{IA(t-1),ia});
    end
    Dbeg.p = adjustd(Dbeg.p,ia,1);  %adjust each period, other wise the error in U is of order 1e-6;

    % Klag = dst_moment(1,Dbeg);  %aggregate K same for_ Dend and Dbeg!
    [Klag,CapUE,MeanCapUE] = dst_moment(1,Dbeg);  %changed 17.1.08
    [DCoh,wage,R] = dbeg2dcoh(Dbeg,ia);
    [Y,wage2,R2] = prodfunc(Klag,Params.Lvec(ia),Params.zAggr(ia));

    if(t==1)  %initialize xequ
      xequ = mean(xEqu(:,:,ia));
      xequ = xequ(:);
    end

    % next values:
    xequ = solverecu(xequ,DCoh,[],[],ia,VAac,1);
    Daac = Params.DaacTry;  %was set in last call to reaction.m;
    if(Params.useDsel)
      D2 = dsf_d(Daac,ia); 
      errorDSF = max(cumsum(Daac.p(:))-cumsum(D2.p(:)));
      Daac = D2;
    else
      errorDSF = 0;
    end

    % compute Euler residuals in period t=1, and then every Params.iEulerres periods:
    if(Params.iEulerres>0 & mod(t,Params.iEulerres)==1)
      P2 = Params;  % save Params;
      for ianext=1:2
	[DNext{ianext},VCohNext{ianext},wNext(ianext),RNext(ianext)] = forwardd(VAac,xequ,Daac,ia,ianext);
      end
      Params = P2;  % restore Params;
      P2 = [];
      % Params.Eulerres always contains the last sample of residuals:
      Params.Eulerres = eulerres(midpoint(DCoh.x),Params.VcohEndog,VCohNext,RNext,wNext,ia,Params.TransProb);
      % here we collect summary statistics:
      maxres = max(abs(Params.Eulerres(:)));
      % meanres = sqrt(sum(Params.Eulerres(:).^2 .* DCoh.p(:) ));
      meanres = mean(mean(abs(Params.Eulerres(1:1000,:)))); %only use the first 1000 points in K, as required by DJJ;
      Params.eulerstats = [Params.eulerstats;[maxres meanres]];
      % disp(sprintf('max and mean euler:  %e  %e',maxres,meanres));
    end

    % compute current moments:
    m = mom2stat(make_hvec(Daac.x)'*Daac.p(:));
    meanK = dst_moment(1,Daac);
    Invest = meanK - (1-Params.delta)*Klag;
    Xtotal = dst_moment(1,DCoh);
    L = Params.Lvec(ia);
    z = Params.zAggr(ia);
    Ytotal = Xtotal - (1-Params.delta)*Klag;  %includes capital income
    assert(max(abs([Ytotal-Y R2-R wage2-wage]))<5e-10);
    xHH = R*kHHlag + Params.YInd{ia}(IEmpl(t))*wage;
    yHH = xHH - kHHlag;
    moms = m(Params.iMomK:end)';
    cHH = cons_at(xHH,Params.VcohEndog,IEmpl(t));
    U = sum(Dbeg.p(:,1));
    assert(abs(U-Params.U(ia))<1e-8);
    % fprintf(1,'t=%d; U = %f, error U:  %e\n',t,U,U-Params.U(ia));

    % COMPUTE HIGHER MOMENTS:
    for j=1:nMoreMoms
      %changed to Dbeg on 27.2.08:
      [moreMoms(j),dummy,condMoms(j,:)] = dst_moment(j,Dbeg,0);  %non-central moments
      if(j>1)
	moreMoms(j) = moreMoms(j).^(1/j) / moreMoms(1);
	condMoms(j,:) = condMoms(j,:).^(1/j) ./ condMoms(1,:);
      end
    end
    Dmarg = dst_merge(Dbeg);  %changed to Dbeg on 17.1.08
    if(0)
      probSmall = dst_distrfunc(Dmarg,[0.05;0.1;0.2;0.3]);
      probSmallUE = dst_distrfunc(Dbeg,[0.05;0.1;0.2;0.3],1); %changed to Dbeg on 17.1.08
    else  %compute percentiles:
      probSmall = dst_quantile(Dmarg,[0.05;0.1;0.2;0.3]);
      probSmallUE = dst_quantile(Dbeg,[0.05;0.1;0.2;0.3]); %changed to Dbeg on 17.1.08
    end
    cSS = 0;  %fill with other variables


    % save distribution and moments:
    if(Params.iSaveDistr & t>tSkip)
      allD{t-tSkip} = Daac;
    end


    % put in time series:
    ser(t,1:Params.nMoments) = moms;
    ser(t,idSer.moremom:idSer.moremom+3*nMoreMoms-1) = [moreMoms condMoms(:)'];
    mom2 = zeros(4,1);


    
    % statistics of Daac:
    C = Xtotal - meanK;
    kHH = xHH - cHH;
    h2 = make_hvec(Daac.x,2);
    % mh2 = mom2stat(h2'*Daac.p(:),2); causes problems!!!
    mh2 = h2'*Daac.p(:);
    ser(t,idSer.m2:idSer.m2+1) = mh2(3:4);
    fracConstr = [Daac.p(1,1)/U Daac.p(1,2)/(1-U) sum(Daac.p(1,:))];
    ser(t,nv+2:nv+nVars) = [meanres maxres errorDSF 0 MeanCapUE  fracConstr C Ytotal R wage L cHH kHHlag yHH IEmpl(t) cSS Klag Invest CapUE probSmall' probSmallUE(:)'];


    kHHlag = kHH;
  end;
  fprintf(1,'\n',t);

  
  if(Params.iMakeForec)
    save(resfilename('S',maketag(Params,0),Params.iter),'ser','idSer');
  else
    save lastsimul ser idSer;
  end

  if(Params.iSaveDistr & T>tSkip)
    tag = maketag(Params,0);
    % iter = Params.iter+1, since distribution counts for_ next round:
    fn = resfilename('D',tag,Params.iter+1);
    save(fn,'allD');
    Params.RefDistrib = fn;
  end

  toc;