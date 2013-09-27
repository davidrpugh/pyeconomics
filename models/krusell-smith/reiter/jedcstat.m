% Matlab file for_ paper "Solving the Incomplete Markets Model With Aggregate Uncertainty by Backward Induction"
% Michael Reiter, Institute for Advanced Studies, August 2008
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% computes statistics from the JEDC comparison paper
function jedcstat(ser,id,tag)
  global Params;

  IA = ser(:,id.ia);
  indxGood = IA==2;
  indxBad = IA==1;
  FF = fopen(sprintf('res/aggrstat_%s.txt',tag),'w');
  fprintf(FF,'correlation of individual and aggregate consumption: %0.6f\n',...
	  corrc(ser(:,id.C),ser(:,id.cHH)));
  fprintf(FF,'correlation of individual consumption and aggregate income: %0.6f\n',...
	  corrc(ser(:,id.Y),ser(:,id.cHH)));
  fprintf(FF,'correlation of individual consumption and aggregate capital stock: %0.6f\n',...
	  corrc(ser(:,id.Klag),ser(:,id.cHH)));
  fprintf(FF,'correlation of individual consumption and individual income: %0.6f\n',...
	  corrc(ser(:,id.yHH),ser(:,id.cHH)));
  fprintf(FF,'correlation of individual consumption and individual capital stock: %0.6f\n',...
	  corrc(ser(:,id.kHH),ser(:,id.cHH)));
  fprintf(FF,'standard deviation of individual consumption: %0.6f\n',...
	  std(ser(:,id.cHH)));
  fprintf(FF,'standard deviation of individual capital: %0.6f\n',...
	  std(ser(:,id.kHH)));
  fprintf(FF,'autocorrelation of individual consumption, lag 1,2,3: %0.6f  %0.6f  %0.6f\n',...
	  ac(ser(:,id.cHH),1),ac(ser(:,id.cHH),2),ac(ser(:,id.cHH),3));
  fprintf(FF,'autocorrelation of individual capital, lag 1,2,3: %0.6f  %0.6f  %0.6f\n',...
	  ac(ser(:,id.kHH),1),ac(ser(:,id.kHH),2),ac(ser(:,id.kHH),3));
  fprintf(FF,'autocorrelation of consumption growth, lag 1,2,3: %0.6f  %0.6f  %0.6f\n',...
	  acd(ser(:,id.cHH),1),acd(ser(:,id.cHH),2),acd(ser(:,id.cHH),3));

  if(Params.whichShocks<2)
    fprintf(FF,'fraction of times an agent is at the constraint: %18.10e\n',...
	    mean(ser(:,id.fracConstr)));
    fprintf(FF,'fraction of times an agent is at the constraint in the good aggregate state: %18.10e\n',...
	    mean(ser(indxGood,id.fracConstr)));
    fprintf(FF,'fraction of times an agent is at the constraint in the bad aggregate state: %18.10e\n',...
	    mean(ser(indxBad,id.fracConstr)));
    fprintf(FF,'average values of the 5th and 10th percentile\n');
    fprintf(FF,'%30s: %18.10e  %18.10e\n','unconditional',mean(ser(:,id.p05)),mean(ser(:,id.p10)));
    fprintf(FF,'%30s: %18.10e  %18.10e\n','unconditional employed',mean(ser(:,id.p05E)),mean(ser(:,id.p10E)));
    fprintf(FF,'%30s: %18.10e  %18.10e\n','unconditional unempl.',mean(ser(:,id.p05U)),mean(ser(:,id.p10U)));
    fprintf(FF,'%30s: %18.10e  %18.10e\n','good state',mean(ser(indxGood,id.p05)),mean(ser(indxGood,id.p10)));
    fprintf(FF,'%30s: %18.10e  %18.10e\n','good state employed',mean(ser(indxGood,id.p05E)),mean(ser(indxGood,id.p10E)));
    fprintf(FF,'%30s: %18.10e  %18.10e\n','good state unempl.',mean(ser(indxGood,id.p05U)),mean(ser(indxGood,id.p10U)));
    fprintf(FF,'%30s: %18.10e  %18.10e\n','bad state',mean(ser(indxBad,id.p05)),mean(ser(indxBad,id.p10)));
    fprintf(FF,'%30s: %18.10e  %18.10e\n','bad state employed',mean(ser(indxBad,id.p05E)),mean(ser(indxBad,id.p10E)));
    fprintf(FF,'%30s: %18.10e  %18.10e\n','bad state unempl.',mean(ser(indxBad,id.p05U)),mean(ser(indxBad,id.p10U)));
  end
  
  nMoreMoms = 5;  % must be consistent with simul.m!!
  fprintf(FF,'AVERAGES AND STANDARD DEVIATIONS OF MOMENTS:\n%5s  %18s  %18s  %18s  %18s  %18s  %18s\n',...
	  'Moment','Ave Population','Std Population','Ave Unemployed','Std Unemployed','Ave Employed','Std Employed');
  for j=1:nMoreMoms
    jm = id.moremom+j-1;
    fprintf(FF,'%d:      %18.10e  %18.10e  %18.10e  %18.10e  %18.10e  %18.10e\n',...
	    j,...
	    mean(ser(:,jm)),std(ser(:,jm)),...
	    mean(ser(:,jm+nMoreMoms)),std(ser(:,jm+nMoreMoms)),...
	    mean(ser(:,jm+2*nMoreMoms)),std(ser(:,jm+2*nMoreMoms)));
  end


  fprintf(FF,'Aggregate time series, original series, undetrended:\n');
  fprintf(FF,'%22s%12s%12s%12s%12s%12s\n','','ave','std','corr(1)','corr(2)','corr(3)');
  bcstat(FF,r_transform(ser(:,id.R)),'rental rate (perc.)',0);
  bcstat(FF,ser(:,id.w),'wage rate',0);
  bcstat(FF,ser(:,id.Y),'aggr. inc.',0);
  bcstat(FF,ser(:,id.C),'aggr. cons.',0);
  bcstat(FF,ser(:,id.Inv),'aggr. investm.',0);

  fprintf(FF,'%22s%12s%12s%12s%12s%12s%12s%12s\n','Corr. of Y(t) with',...
	  'x(t-3)','x(t-2)','x(t-1)','x(t)','x(t+1)','x(t+2)','x(t+3)');
  corrY(FF,r_transform(ser(:,id.R)),'rental rate (perc.)',ser(:,id.Y),0);
  corrY(FF,ser(:,id.w),'wage rate',ser(:,id.Y),0);
  corrY(FF,ser(:,id.Y),'aggr. inc.',ser(:,id.Y),0);
  corrY(FF,ser(:,id.C),'aggr. cons.',ser(:,id.Y),0);
  corrY(FF,ser(:,id.Inv),'aggr. investm.',ser(:,id.Y),0);

%  fprintf(FF,'Aggregate time series, logged and HPfiltered(1600):\n');
%  fprintf(FF,'%22s%12s%12s%12s%12s%12s\n','','ave','std','corr(1)','corr(2)','corr(3)');
%  bcstat(FF,r_transform(ser(:,id.R)),'rental rate (perc.)',1);
%  bcstat(FF,ser(:,id.w),'wage rate',1);
%  bcstat(FF,ser(:,id.Y),'aggr. inc.',1);
%  bcstat(FF,ser(:,id.C),'aggr. cons.',1);
%  bcstat(FF,ser(:,id.Inv),'aggr. investm.',1);

%  fprintf(FF,'%22s%12s%12s%12s%12s%12s%12s%12s\n','Corr. of Y(t) with',...
%	  'x(t-3)','x(t-2)','x(t-1)','x(t)','x(t+1)','x(t+2)','x(t+3)');
%  corrY(FF,r_transform(ser(:,id.R)),'rental rate (perc.)',ser(:,id.Y),1);
%  corrY(FF,ser(:,id.w),'wage rate',ser(:,id.Y),1);
%  corrY(FF,ser(:,id.Y),'aggr. inc.',ser(:,id.Y),1);
%  corrY(FF,ser(:,id.C),'aggr. cons.',ser(:,id.Y),1);
%  corrY(FF,ser(:,id.Inv),'aggr. investm.',ser(:,id.Y),1);


  fclose(FF);

  FF = fopen(sprintf('res/ts_%s.txt',tag),'w');
  for t=1:size(ser,1)
    % aggregate state (just to check we are doing the same thing)
    % individual state (again just to check)
    % begining-of-period capital of our one individual (43 in first period/row)
    % individual consumption choice
    % for_ the beginning-of-period distribution
    %   mean capital stock of the unemployed
    %   mean capital stock of the employed
    %   rental and wage rate as implied by these two means
    %   5th percentile of the unemployed
    %   10th percentile of the unemployed
    %   5th percentile of the employed
    %   10th percentile of the employed

    fprintf(FF,'%d  %d  %18.10e  %18.10e  %18.10e  %18.10e  %18.10e  %18.10e  %18.10e  %18.10e  %18.10e  %18.10e\n',...
	    ser(t,id.ia),ser(t,id.ie),ser(t,id.kHH),ser(t,id.cHH),...
	    ser(t,id.meanU),ser(t,id.meanE),r_transform(ser(t,id.R)),ser(t,id.w),...
	    ser(t,id.p05U),ser(t,id.p10U),ser(t,id.p05E),ser(t,id.p10E));
  end
  fclose(FF);


function c = corrc(x,y)
  m = corrcoef(x,y);
  c = m(1,2);

function c = corrlag(x,y,ilag)
  n = length(x);
  c = corrc(x(1+ilag:n),y(1:n-ilag));

function c = ac(x,ilag)
  c = corrlag(x,x,ilag);

function c = acd(x,ilag)
  n = length(x);
  y = x(2:n)./x(1:n-1) - 1;
  c = ac(y,ilag);

function bcstat(FF,x,name,doHP)
  if(doHP)
    x = log(x);
    x = x - hpfilter(x,1600);
  end
 fprintf(FF,'%20s:  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f\n',...
	 name,mean(x),std(x),ac(x,1),ac(x,2),ac(x,3));
 %Average and standard deviation 
 %Autocorrelation of both prices (up to three lags)

function corrY(FF,x,name,y,doHP)
  % • Standard table of auto and cross correlations (with output) at leads and lags (up to three)
  if(doHP)
    x = log(x);
    y = log(y);
    x = x - hpfilter(x,1600);
    y = y - hpfilter(y,1600);
  end
  fprintf(FF,'%20s:  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f\n',...
	  name,corrlag(y,x,3),corrlag(y,x,2),corrlag(y,x,1),corrlag(y,x,0),...
	  corrlag(x,y,1),corrlag(x,y,2),corrlag(x,y,3));


