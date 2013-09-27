% The program for the article "Solving the incomplete markets model with
% aggregate uncertainty using the Krusell-Smith algorithm" from the special 
% JEDC issue edited by Den Haan, Judd and Juillard (2008)  
%
% Written by Lilia Maliar, Serguei Maliar and Fernando Valli (2008)
%
% "TEST.m" (should be run after "MAIN.m"; it uses "Inputs_for_test" and
% "Solution_to_model" for computing the statistics reported in Den Haan's
% 2008, comparison article)  
%__________________________________________________________________________
clc;
clear all;
%__________________________________________________________________________

load Solution_to_model; % contains a solution to the model computed by MAIN.m
load Inputs_for_test;   % contains three objects provided by Den Haan, Judd 
                        % and Juillard, 2008): 
                        % 1) "kdist" is initial distribution of capital on a 
                        % 1000-point grid for unemployed and employed agents; 
                        % 2) "agshock" is a 10,000-period realization of 
                        % aggregate shock;
                        % 3) "idshock" is a 10,000-period realization of 
                        % idiosyncratic shock for one agent;
%__________________________________________________________________________
%
% 1. PARAMETERS
%__________________________________________________________________________

T=10000; % simulation length
J=1001;  % number of grid points for non-stochastic simulation (is equal to 
         % 1001 and not to 1000 because we include explicitly point k=0)
ur=zeros(nstates_ag,1); % vector of unemployment rates in two agg. states
ur(1)=ur_b; ur(2)=ur_g; 
%__________________________________________________________________________
%
% 2. SERIES
%__________________________________________________________________________

% Time series of aggregate variables

km_ag=zeros(T,1);       % aggregate capital (mean of capital distribution)
c_ag=zeros(T-1,1);      % aggregate consumption
income_ag=zeros(T-1,1); % aggregate income
wealth_ag=zeros(T-1,1); % aggregate wealth
i_ag=zeros(T-1,1);      % aggregate net investment

prod_ag=(agshock==1)*a(1)+(agshock==2)*a(2);   % aggregate productivity
labor_ag=(agshock==1)*er_b+(agshock==2)*er_g;  % aggregate employment (L in 
                                               % the paper)

irate=zeros(T-1,1); % interest rate
wage=zeros(T-1,1);  % wage

% Time series of individual variables

k_ind=zeros(T,1);        % individual capital
k_ind(1,1)=43;           % initial individual capital 
c_ind=zeros(T-1,1);      % individual consumption
income_ind=zeros(T-1,1); % individual income

% Beginning-of-period capital distributions

kdistu=zeros(T,J); % beginning-of-period capital distributions in all 
                   % periods for the unemployed
kdiste=zeros(T,J); % for the employed

% Initial-period capital distributions

% kdist has dimensionality J*2 with J=1001; in the last point J, the probability 
% is chosen to normalize the sum of probabilities to 1

kdistu(1,1:J)=[kdist(1:end-1,1)' 1-sum(kdist(1:end-1,1))]; % a raw vector of 
                       % the initial capital distribution for the unemployed
kdiste(1,1:J)=[kdist(1:end-1,2)' 1-sum(kdist(1:end-1,2))]; % for the employed

% Grid for the capital distribution

kvalues_min=0;   % minimum grid value
kvalues_max=100; % maximum grid value
kvalues=linspace(kvalues_min,kvalues_max,J)'; % grid for capital distribution
                                              % with J grid points
% _________________________________________________________________________
%
% 3. TRANSITION PROBABILITIES 
%__________________________________________________________________________

% prob_ag(i,j) is the probability of tomorrow's agg. shock (i=1,2) given 
% today's agg. shock (j=1,2)

prob_ag=zeros(2,2);
prob_ag(1,1)=prob(1,1)+prob(1,2); prob_ag(2,1)=1-prob_ag(1,1);  
prob_ag(2,2)=prob(3,3)+prob(3,4); prob_ag(1,2)=1-prob_ag(2,2);

% p_xy_zw is the probability of idiosyncratic shock epsilon'=w conditional 
% on aggregate shocks s'=y, s=x and idiosyncratic shock epsilon=z 

p_bb_uu = prob(1,1)/prob_ag(1,1); p_bb_ue=1-p_bb_uu;
p_bb_ee = prob(2,2)/prob_ag(1,1); p_bb_eu=1-p_bb_ee;
p_bg_uu = prob(1,3)/prob_ag(2,1); p_bg_ue=1-p_bg_uu;
p_bg_ee = prob(2,4)/prob_ag(2,1); p_bg_eu=1-p_bg_ee;
p_gb_uu = prob(3,1)/prob_ag(1,2); p_gb_ue=1-p_gb_uu;
p_gb_ee = prob(4,2)/prob_ag(1,2); p_gb_eu=1-p_gb_ee;
p_gg_uu = prob(3,3)/prob_ag(2,2); p_gg_ue=1-p_gg_uu;
p_gg_ee = prob(4,4)/prob_ag(2,2); p_gg_eu=1-p_gg_ee;

% Transition probabilities from one idiosyncratic shock, epsilon, to another, 
% epsilon', given that agg. shocks are s and s'

p_bb = [[p_bb_uu p_bb_ue]*ur_b; [p_bb_eu p_bb_ee]*er_b];
p_bg = [[p_bg_uu p_bg_ue]*ur_b; [p_bg_eu p_bg_ee]*er_b];
p_gb = [[p_gb_uu p_gb_ue]*ur_g; [p_gb_eu p_gb_ee]*er_g];
p_gg = [[p_gg_uu p_gg_ue]*ur_g; [p_gg_eu p_gg_ee]*er_g];

%__________________________________________________________________________
%
% 4. NON-STOCHASTIC SIMULATION
%__________________________________________________________________________

for t=1:T
    
   t
   
   % Aggregate capital
   
   km_ag(t)=kdistu(t,:)*kvalues*(1-labor_ag(t))+kdiste(t,:)*kvalues*labor_ag(t);
      % aggregate capital=capital of the unemployed + capital of the employed;
      % (1-labor_ag(t)) is the share of unemployed people in the economy
   
   % Prices
   
   irate(t)=alpha*prod_ag(t)*(km_ag(t)/labor_ag(t)/l_bar)^(alpha-1);
      % interest rate 
   wage(t)=(1-alpha)*prod_ag(t)*(km_ag(t)/labor_ag(t)/l_bar)^alpha;
      % wage
  %________________________________________________________________________
  %
  % Individual capital function, k'
  %________________________________________________________________________
   
  kprimet(:,1)=interpn(k,km,kprime(:,:,agshock(t),1),kvalues,km_ag(t)*ones(J,1),'linear'); 
      % interpolate the capital function k' (computed in "MAIN") of the 
      % unemployed agent in kvalues for the given agg. capital km_ag(t)
      
  kprimet(:,2)=interpn(k,km,kprime(:,:,agshock(t),2),kvalues,km_ag(t)*ones(J,1),'linear');
      % the same for the employed agent
      
   kprimet=kprimet.*(kprimet>=kvalues_min).*(kprimet<=kvalues_max)+kvalues_min*(kprimet<kvalues_min)+kvalues_max*(kprimet>kvalues_max);
      % restrict individual capital k' to be within [kvalues_min, kvalues_max]
      % range
   %_______________________________________________________________________
   %
   % Inverse of the individual capital function, k'(x) 
   %_______________________________________________________________________
   
   % To invert k'(x), we treat k' as an argument and x as a value of the 
   % function in the argument k'
   % Note that an inverse of k' is not well defined for those points of the 
   % grid for which the value of k' is the same (for example, for unemployed 
   % agents, we have k'=0 for several small grid values) 
   % Therefore, when inverting, we remove all but one grid points for which
   % the values of k' are repeated
   
   index_min=zeros(1, 2); % this variable will count how many times k'=kvalues_min 
      % in the employed and the unemployed states separately
      
   index_min=sum(kprimet==kvalues_min); % count how many times k'=kvalues_min
   
   first=index_min.*(index_min>0)+1*(index_min==0); % if index_min>0, consider  
      % k' starting from the (index_min)-th grid point; otherwise, consider  
      % k' starting from the first grid point 
      
   index_max=zeros(1, 2); % this variable will count how many times k'=kvalues_max 
      % in the employed and the unemployed states 
      
   index_max=sum(kprimet==kvalues_max); % count how many times k'=kvalues_max
   
   last=J-((index_max-1).*(index_max>0)+0*(index_max==0)); % if index_max>0, 
      %consider k' until the grid point (J-(index_max-1)); otherwise, 
      % consider k' until the last grid point, which is J 

   xt(:,1)=interp1(kprimet(first(1):last(1),1),kvalues(first(1):last(1),1), kvalues,'linear'); 
      % find x(k') in the unemployed state (state 1) by interpolation (see 
      % condition (10) in the paper)
      
   xt(:,2)=interp1(kprimet(first(2):last(2),2),kvalues(first(2):last(2),1), kvalues,'linear');
      % find x(k') in the employed state (state 2) by interpolation

   xt=xt.*(xt>=kvalues_min).*(xt<=kvalues_max)+kvalues_min*(xt<kvalues_min)+kvalues_max*(xt>kvalues_max);
      % restrict k to be in [kvalues_min, kvalues_max] 
      
   % Notice that some values of xt at the beginning and in the end of the 
   % grid will be NaN. This is because there are no values of kprimet 
   % (i.e., k') that correspond to some values of xt (i.e., x). 
   % For example, to have kprimet(xt)=0 for an employed agent xt must be 
   % negative, which is not allowed, so we get NaN. These NaN values of xt 
   % create a problem when computing the end-of-period capital distribution 
   % for terminal (but not for initial) grid-value points. To deal with this 
   % problem, we set xt in the end of the grid to kvalues_max whenever they 
   % are NaN. 

   % unemployed (consider xt(:,1))
   j=0; 
   while isnan(xt(J-j,1))==1; % consider xt=NaN from the end of the grid (j=0)
      xt(J-j,1)=kvalues_max;  % when xt=NaN, set xt=kvalues_max
      j=j+1;
   end

   % employed (consider xt(:,2))
   j=0; 
   while isnan(xt(J-j,2))==1;
      xt(J-j,2)=kvalues_max;
      j=j+1;
   end
%__________________________________________________________________________
%
% End-of-period cumulative capital distribution 
%__________________________________________________________________________

  Fu=zeros(J,1); % the end-of-period capital distribution for the unemployed; 
                 % note that we do not store the cumulative density over time   
  Fe=zeros(J,1); % the end-of-period capital distribution for the employed;

  for j=1:J; % we have J values of xt(j,1) and J values of xt(j,2) which are 
             % inverse of kprimet(kvalues) for the unemployed and employed
      
     % unemployed
              for i=1:J; % consider points on the grid (kvalues)
              if kvalues(i,1)<=xt(j,1); % if a grid point i considered is 
                  % <= xt(j,1)
                 Fu(j)=Fu(j)+kdistu(t,i); % then compute the cumulative 
                 % distribution Fu by adding the probability between the 
                 % points kvalues(i-1,1) and kvalues(i,1)                  
              end
              if kvalues(i,1)>xt(j,1); % if a grid point i considered is 
                  % >xt(j,1)
                 Fu(j)=Fu(j)+(xt(j,1)-kvalues(i-1,1))/(kvalues(i,1)-kvalues(i-1,1))*kdistu(t,i); break
                  % then compute the cumulative distribution Fe by adding 
                  % the probability to be the point between points kvalues(i-1,1) 
                  % and xt(j,1) and then break
              end     
       end
       
       % employed
              for i=1:J; 
              if kvalues(i,1)<=xt(j,2);
                  Fe(j)=Fe(j)+kdiste(t,i); 
               end
               if kvalues(i,1)>xt(j,2);
                 Fe(j)=Fe(j)+(xt(j,2)-kvalues(i-1,1))/(kvalues(i,1)-kvalues(i-1,1))*kdiste(t,i); break  
               end
         end
  end
%__________________________________________________________________________
%
% Next period's beginning-of-period distribution
%__________________________________________________________________________
 
if t < T % we do not compute next period's beginning-of-period distributions 
         % for t=T (i.e., kdistu(T+1,:) and kdiste(T+1,:)) as we do not have 
         % agshock(T+1) 
     
%  Mass of agents in different idiosyncratic states conditional on agg.
%  states

   if (agshock(t)==1)&(agshock(t+1)==1);g=p_bb; end % g is a 2*2 matrix; 
   % for example, g(1,2) is the mass of agents who were unemployed at t and 
   % employed at t+1, conditional on being in a bad agg. state at both t 
   % and t+1
   if (agshock(t)==1)&(agshock(t+1)==2);g=p_bg; end
   if (agshock(t)==2)&(agshock(t+1)==1);g=p_gb; end
   if (agshock(t)==2)&(agshock(t+1)==2);g=p_gg; end

% Next period's beginning-of-period distribution (see formulas in Den Haan, 
% Judd and Juillard, 2008)    
   
   Pu=(g(1,1)*Fu+g(2,1)*Fe)/(g(1,1)+g(2,1)); % for the unemployed
   Pe=(g(1,2)*Fu+g(2,2)*Fe)/(g(1,2)+g(2,2)); % for the employed
   
  % unemployed 
   kdistu(t+1,1)=Pu(1,1); % probability of having k=kvalues_min at t+1
   kdistu(t+1,2:J-1)=Pu(2:J-1,1)'-Pu(1:J-2,1)'; % probabilities of different 
                                                % grid points kvalues
   kdistu(t+1,J)=1-sum(kdistu(t+1,1:J-1));      % probability of k=kvalues_max 
                                                % is set so that "kdistu" is
                                                % normalized to one
     
  % employed 
   kdiste(t+1,1)=Pe(1,1); 
   kdiste(t+1,2:J-1)=Pe(2:J-1,1)'-Pe(1:J-2,1)';
   kdiste(t+1,J)=1-sum(kdiste(t+1,1:J-1));
 end
end; % end of the NON-STOCHASTIC SIMULATION
%__________________________________________________________________________
%
% 5. QUANTITIES
%__________________________________________________________________________

for t=1:T-1; 
    
% Individual capital

k_ind(t+1,1)=interpn(k,km,kprime(:,:,agshock(t),idshock(t,1)),k_ind(t,1),km_ag(t),'linear');

k_ind(t+1,1)=k_ind(t+1,1).*(k_ind(t+1,1)>=k_min).*(k_ind(t+1,1)<=k_max)+k_min*(k_ind(t+1,1)<k_min)+k_max*(k_ind(t+1,1)>k_max);
           % restrict k_ind to be in [k_min,k_max] range

% Individual income

income_ind(t,1)=k_ind(t,1)*irate(t)+(idshock(t,1)-1).*l_bar*wage(t)+mu*(2-idshock(t,1)).*wage(t)-...
ur(agshock(t))/(1-ur(agshock(t)))*mu*(idshock(t,1)-1).*wage(t);
 
% Individual consumption

c_ind(t,1)=k_ind(t,1)*(1-delta+irate(t))+(idshock(t,1)-1).*l_bar*wage(t)+mu*(2-idshock(t,1)).*wage(t)-...
ur(agshock(t))/(1-ur(agshock(t)))*mu*(idshock(t,1)-1).*wage(t)-k_ind(t+1,1);

end 

% Output
output=prod_ag(1:T-1,1).*km_ag(1:T-1,1).^alpha.*(labor_ag(1:T-1,1)*l_bar).^(1-alpha);

% Aggregate consumption
c_ag=km_ag(1:T-1,1)*(1-delta)+output(1:T-1,1)-km_ag(2:T,1);

% Aggregate income
income_ag=output(1:T-1,1);

% Aggregate investment
i_ag=km_ag(2:T)-(1-delta)*km_ag(1:T-1);
%__________________________________________________________________________
%
% 6. STATISTICS
%__________________________________________________________________________
% 
% Capital distribution
%__________________________________________________________________________
 
kdist=zeros(T,J);  % capital distribution for both employed and unemployed

 for t=1:T

 kdist(t,:)=kdistu(t,:)*(1-labor_ag(t))+kdiste(t,:)*labor_ag(t);

 % Cumulative distribution 
 
 cdist=cumsum(kdistu(t,:)*(1-labor_ag(t))+kdiste(t,:)*labor_ag(t));
                             % both employed and unemployed
 cdistu=cumsum(kdistu(t,:)); % unemployed
 cdiste=cumsum(kdiste(t,:)); % employed
 
 % 5 and 10 percentiles employed and unemployed
  
 per5(t,1)=interp1(cdist(1:J-100),kvalues(1:J-100),0.05);
                                        % both 
 per5u(t,1)=interp1(cdistu(1:J-100),kvalues(1:J-100),0.05);
                                        % unemployed
 per5e(t,1)=interp1(cdiste(1:J-100),kvalues(1:J-100),0.05);
                                        % employed
 per10(t,1)=interp1(cdist(1:J-100),kvalues(1:J-100),0.1);
 per10u(t,1)=interp1(cdistu(1:J-100),kvalues(1:J-100),0.1);
 per10e(t,1)=interp1(cdiste(1:J-100),kvalues(1:J-100),0.1);
  
 % 1-st and normalized 2,3,4 and 5th moments of capital distribution of 
 % employed, unemployed and both
 
 m1(t,1)=km_ag(t); % both
 m1u(t,1)=kdistu(t,:)*kvalues; % unemployed
 m1e(t,1)=kdiste(t,:)*kvalues; % employed

 m2(t,1)=(kdist(t,:)*kvalues.^2)^(1/2)/m1(t);
 m2u(t,1)=(kdistu(t,:)*kvalues.^2)^(1/2)/m1u(t);
 m2e(t,1)=(kdiste(t,:)*kvalues.^2)^(1/2)/m1e(t);
  
 m3(t,1)=(kdist(t,:)*kvalues.^3)^(1/3)/m1(t);
 m3u(t,1)=(kdistu(t,:)*kvalues.^3)^(1/3)/m1u(t);
 m3e(t,1)=(kdiste(t,:)*kvalues.^3)^(1/3)/m1e(t);

 m4(t,1)=(kdist(t,:)*kvalues.^4)^(1/4)/m1(t);
 m4u(t,1)=(kdistu(t,:)*kvalues.^4)^(1/4)/m1u(t);
 m4e(t,1)=(kdiste(t,:)*kvalues.^4)^(1/4)/m1e(t);

 m5(t,1)=(kdist(t,:)*kvalues.^5)^(1/5)/m1(t);
 m5u(t,1)=(kdistu(t,:)*kvalues.^5)^(1/5)/m1u(t);
 m5e(t,1)=(kdiste(t,:)*kvalues.^5)^(1/5)/m1e(t); 
    
 end

% Table 9 in Den Haan's (2008) comparison paper

Table9=cell(22,2);

Table9(1,1)= {'sample means'};
Table9(1,2)= {' '};

% Sample means

Table9(2,1)= {'1st unemp'};Table9(2,2)={mean(m1u)};
Table9(3,1)= {'2nd unemp'};Table9(3,2)={mean(m2u)};
Table9(4,1)= {'3rd unemp'};Table9(4,2)={mean(m3u)};
Table9(5,1)= {'4th unemp'};Table9(5,2)={mean(m4u)};
Table9(6,1)= {'5th unemp'};Table9(6,2)={mean(m5u)};
Table9(7,1)= {'1st emp'};Table9(7,2)={mean(m1e)};
Table9(8,1)= {'2nd emp'};Table9(8,2)={mean(m2e)};
Table9(9,1)= {'3rd emp'};Table9(9,2)={mean(m3e)};
Table9(10,1)= {'4th emp'};Table9(10,2)={mean(m4e)};
Table9(11,1)= {'5th emp'};Table9(11,2)={mean(m5e)};

Table9(12,1)= {'sample std'};
Table9(12,2)= {' '};

% Sample standard deviations

Table9(13,1)= {'1st unemp'};Table9(13,2)={std(m1u)};
Table9(14,1)= {'2nd unemp'};Table9(14,2)={std(m2u)};
Table9(15,1)= {'3rd unemp'};Table9(15,2)={std(m3u)};
Table9(16,1)= {'4th unemp'};Table9(16,2)={std(m4u)};
Table9(17,1)= {'5th unemp'};Table9(17,2)={std(m5u)};
Table9(18,1)= {'1st emp'};Table9(18,2)={std(m1e)};
Table9(19,1)= {'2nd emp'};Table9(19,2)={std(m2e)};
Table9(20,1)= {'3rd emp'};Table9(20,2)={std(m3e)};
Table9(21,1)= {'4th emp'};Table9(21,2)={std(m4e)};
Table9(22,1)= {'5th emp'};Table9(22,2)={std(m5e)};

% Table 11 in Den Haan's (2008) comparison paper

Table11=cell(14,2);

SeriesT11=[agshock per5u per10u per5e per10e]; 
series_sort=sortrows(SeriesT11); % sort series by the aggregate shock
Tb=sum(series_sort(:,1)==1);     % number of bad realizations of agg. shock
mean_all=mean(series_sort(1:T,:));     % unconditional means
mean_bad=mean(series_sort(1:Tb,:));    % bad-state means
mean_good=mean(series_sort(Tb+1:T,:)); % good-state means

Table11(1,1)= {'unemployed'};
Table11(1,2)= {' '};
Table11(2,1)= {'5%'};       Table11(2,2)={mean_all(2)};
Table11(3,1)= {'5%, good'}; Table11(3,2)={mean_good(2)};
Table11(4,1)= {'5%, bad'};  Table11(4,2)={mean_bad(2)};
Table11(5,1)= {'10%'};      Table11(5,2)={mean_all(3)};
Table11(6,1)= {'10%, good'};Table11(6,2)={mean_good(3)};
Table11(7,1)= {'10%, bad'}; Table11(7,2)={mean_bad(3)};

Table11(8,1)= {'employed'};
Table11(8,2)= {' '};
Table11(9,1)= {'5%'};        Table11(9,2)={mean_all(4)};
Table11(10,1)= {'5%, good'}; Table11(10,2)={mean_good(4)};
Table11(11,1)= {'5%, bad'};  Table11(11,2)={mean_bad(4)};
Table11(12,1)= {'10%'};      Table11(12,2)={mean_all(5)};
Table11(13,1)= {'10%, good'};Table11(13,2)={mean_good(5)};
Table11(14,1)= {'10%, bad'}; Table11(14,2)={mean_bad(5)};

%__________________________________________________________________________
%
% Properties of individual policy rules 
%__________________________________________________________________________

% Table 12 in Den Haan's (2008) comparison paper

Table12=cell(19,2);
Table12(1,1)= {'correlations'};
Table12(1,2)= {' '};

% Correlation between individual and aggregate consumption

Table12(2,1)= {'c_ind,c_ag'};
aux1=corrcoef(c_ind,c_ag);
Table12(2,2)={aux1(1,2)};
   
% Correlation between individual consumption and aggregate income

Table12(3,1)= {'c_ind,income_ag'};
aux1=corrcoef(c_ind,income_ag);
Table12(3,2)={aux1(1,2)};
   
% Correlation between individual consumption and aggregate capital

Table12(4,1)= {'c_ind,km_ag'};
aux1=corrcoef(c_ind,km_ag(1:T-1,1));
Table12(4,2)={aux1(1,2)};
   
% Correlation between individual consumption and individual capital

Table12(5,1)= {'c_ind,km_ind'};
aux1=corrcoef(c_ind,k_ind(1:T-1,1));
Table12(5,2)={aux1(1,2)};

Table12(6,1)= {'autocorrelations'};
Table12(6,2)= {' '};
   
% Autocorrelation of individual consumption (up to three lags)

Table12(7,1)= {'c_ind(t),c_ind(t-1)'};
aux1=corrcoef(c_ind(2:T-1,1),c_ind(1:T-2,1));
Table12(7,2)={aux1(1,2)}; % lag=1
Table12(8,1)= {'c_ind(t),c_ind(t-2)'};
aux1=corrcoef(c_ind(3:T-1,1),c_ind(1:T-3,1));
Table12(8,2)={aux1(1,2)}; % lag=2
Table12(9,1)= {'c_ind(t),c_ind(t-3)'};
aux1=corrcoef(c_ind(4:T-1,1),c_ind(1:T-4,1));
Table12(9,2)={aux1(1,2)}; % lag=3
   
% Autocorrelation of individual capital (up to three lags)   

Table12(10,1)= {'k_ind(t),k_ind(t-1)'};
aux1=corrcoef(k_ind(2:T-1,1),k_ind(1:T-2,1));
Table12(10,2)={aux1(1,2)}; % lag=1
Table12(11,1)= {'k_ind(t),k_ind(t-1)'};
aux1=corrcoef(k_ind(3:T-1,1),k_ind(1:T-3,1));
Table12(11,2)={aux1(1,2)}; % lag=2
Table12(12,1)= {'k_ind(t),k_ind(t-1)'};
aux1=corrcoef(k_ind(4:T-1,1),k_ind(1:T-4,1));
Table12(12,2)={aux1(1,2)}; % lag=3
   
% Autocorrelation of individual consumption growth

Table12(13,1)= {'c_gr(t),c_gr(t-1)'};
c_gr=c_ind(2:T-1,1)./c_ind(1:T-2,1);
aux1=corrcoef(c_gr(2:T-2,1),c_gr(1:T-3,1));
Table12(13,2)={aux1(1,2)};   

Table12(14,1)= {'means'};
Table12(14,2)= {' ' };

% Mean of individual consumption

Table12(15,1)= {'c_ind'};
Table12(15,2)={mean(c_ind(1:T-1,1))};   

% Mean of individual capital

Table12(16,1)= {'k_ind'};
Table12(16,2)={mean(k_ind(1:T-1,1))};

Table12(17,1)= {'standard deviations'};
Table12(17,2)= {' '};

% Standard deviation of individual consumption

Table12(18,1)= {'c_ind'};
Table12(18,2)={std(c_ind(1:T-1,1))};   

% Standard deviation of individual capital

Table12(19,1)= {'k_ind'};
Table12(19,2)={std(k_ind(1:T-1,1))};
%__________________________________________________________________________
%
% 7. DYNAMIC EULER EQUATION ACCURACY TEST
%__________________________________________________________________________

k_tilda(1,1) = k_ind(1,1); % initially, k_tilda (alternative series) is 
                           % equal to individual capital found before

   for t=1:T-1
       
      t
       
      k_hat(t+1,1)=interpn(k,km,kprime(:,:,agshock(t),idshock(t,1)),k_tilda(t,1),km_ag(t),'linear');
      % k_hat is a temporary variable found by interpolation of the individual
      % capital function k', i.e., k_hat(t+1)=k'(k_tilda, epsilon)
              
     % Future (at t+1) aggregate capital found from the ALM 
     
     if agshock(t)==1 % bad state
      kmalm_hat(t+1)=exp(B(1)+B(2)*log(km_ag(t)));
     else % good state
      kmalm_hat(t+1)=exp(B(3)+B(4)*log(km_ag(t)));
     end

     % Future (at t+1) interest rate
     
      irate_b=alpha*a(1)*(kmalm_hat(t+1)/er_b/l_bar)^(alpha-1); % bad state
      irate_g=alpha*a(2)*(kmalm_hat(t+1)/er_g/l_bar)^(alpha-1); % good state
      
     % Future (at t+1) wage 
     
      wage_b=(1-alpha)*a(1)*(kmalm_hat(t+1)/er_b/l_bar)^(alpha); % bad state
      wage_g=(1-alpha)*a(2)*(kmalm_hat(t+1)/er_g/l_bar)^(alpha); % good state

      % Future value of the temporary variable k_hat (i.e., k'') in four 
      % possible future states
      
      k_hat_bu(t+2,1)=interpn(k,km,a2,epsilon2,kprime,k_hat(t+1,1),kmalm_hat(t+1),1,1,'linear');
      k_hat_be(t+2,1)=interpn(k,km,a2,epsilon2,kprime,k_hat(t+1,1),kmalm_hat(t+1),1,2,'linear');
      k_hat_gu(t+2,1)=interpn(k,km,a2,epsilon2,kprime,k_hat(t+1,1),kmalm_hat(t+1),2,1,'linear');      
      k_hat_ge(t+2,1)=interpn(k,km,a2,epsilon2,kprime,k_hat(t+1,1),kmalm_hat(t+1),2,2,'linear');      
     
      % Future (at t+1) marginal utility, c(t+1)^(-gamma), with c(t+1)
      % found from the budget constraint
      
      mut_bu(t+1,1)=((1-delta+irate_b)*k_hat(t+1,1)+wage_b*l_bar*epsilon_u+mu*wage_b*(1-epsilon_u)-mu*wage_b*epsilon_u*ur(1)/(1-ur(1))-k_hat_bu(t+2,1))^(-gamma);
      mut_be(t+1,1)=((1-delta+irate_b)*k_hat(t+1,1)+wage_b*l_bar*epsilon_e+mu*wage_b*(1-epsilon_e)-mu*wage_b*epsilon_e*ur(1)/(1-ur(1))-k_hat_be(t+2,1))^(-gamma);
      mut_gu(t+1,1)=((1-delta+irate_g)*k_hat(t+1,1)+wage_g*l_bar*epsilon_u+mu*wage_g*(1-epsilon_u)-mu*wage_g*epsilon_u*ur(2)/(1-ur(2))-k_hat_gu(t+2,1))^(-gamma);
      mut_ge(t+1,1)=((1-delta+irate_g)*k_hat(t+1,1)+wage_g*l_bar*epsilon_e+mu*wage_g*(1-epsilon_e)-mu*wage_g*epsilon_e*ur(2)/(1-ur(2))-k_hat_ge(t+2,1))^(-gamma);
       
      % Conditional expectation (on the right-hand-side of the Euler equation)
      
      curr_state=(agshock(t)-1)*2+idshock(t,1); % supplementary variable for 
                                                % computing the transition 
                                                % probabilities from "prob" 

      E_tilda(t,1)=mut_bu(t+1,1)*(1-delta+irate_b)*prob(curr_state,1)+mut_be(t+1,1)*(1-delta+irate_b)*prob(curr_state,2)+...
          mut_gu(t+1,1)*(1-delta+irate_g)*prob(curr_state,3)+ mut_ge(t+1,1)*(1-delta+irate_g)*prob(curr_state,4);
      
      % Current consumption found from the Euler equation 
      
      c_tilda(t,1)=(beta*E_tilda(t,1))^(-1/gamma); 
      
      % Future (at t+1) capital restored from the budget constraint 
      
      k_tilda(t+1,1)=k_tilda(t,1)*(1-delta+irate(t))+(idshock(t,1)-1).*l_bar*wage(t)+mu*(2-idshock(t,1)).*wage(t)-...
      ur(agshock(t))/(1-ur(agshock(t)))*mu*(idshock(t,1)-1).*wage(t)-c_tilda(t,1);
  
      % Check that k_tilda is not negative 
      
      k_tilda(t+1,1)=k_tilda(t+1,1)*(k_tilda(t+1,1)>0); % set k_tilda=0 
                                                        % if it is negative
  
      % Re-compute c_tilda whenever k_tilda was initially negative 
      
      c_tilda(t,1)=k_tilda(t,1)*(1-delta+irate(t))+(idshock(t,1)-1).*l_bar*wage(t)+mu*(2-idshock(t,1)).*wage(t)-...
      ur(agshock(t))/(1-ur(agshock(t)))*mu*(idshock(t,1)-1).*wage(t)-k_tilda(t+1,1);

   end
   
% Dynamic Euler equation test errors

errk=abs((k_ind(1:T,1)-k_tilda(1:T,1))./mean(k_ind(1:T,1)))*100;
errc=abs((c_ind(1:T-1,1)-c_tilda(1:T-1,1))./c_ind(1:T-1,1))*100;

% Table 14 in Den Haan's (2008) comparison paper

Table14(1,1)= {'capital (scl)'};
Table14(1,2)= {' '};
Table14(2,1)= {'Average'}; Table14(2,2)={mean(errk)};
Table14(3,1)= {'Maximum'}; Table14(3,2)={max(errk)};

Table14(4,1)= {'consumption'};
Table14(4,2)= {' '};
Table14(5,1)= {'Average'}; Table14(5,2)={mean(errc)};
Table14(6,1)= {'Maximum'}; Table14(6,2)={max(errc)};

%__________________________________________________________________________
%
% 8. ACCURACY OF AGGREGATE LAW OF MOTION (ALM)
%__________________________________________________________________________

% Mean of  the capital distribution computed from the ALM

kmalm=zeros(T,1); % this variable represents km computed from the ALM
kmalm(1)=km_ag(1);
for ii=1:T-1 
   if agshock(ii)==1
      kmalm(ii+1)=exp(B(1)+B(2)*log(kmalm(ii)));
   else
      kmalm(ii+1)=exp(B(3)+B(4)*log(kmalm(ii)));
   end
   
end

% Table 15 in Den Haan's (2008) comparison paper

Table15(1,1)= {'mean of K'};
Table15(1,2)= {' '};
Table15(2,1)= {'in panel'}; Table15(2,2)={mean(km_ag)};
Table15(3,1)= {'from alm'}; Table15(3,2)={mean(kmalm)};
Table15(4,1)= {'% differ'}; Table15(4,2)={abs(mean(km_ag)-mean(kmalm))/mean(km_ag)*100};

Table15(5,1)= {'std of K'};
Table15(5,2)= {' '};
Table15(6,1)= {'in panel'}; Table15(6,2)={std(km_ag)};
Table15(7,1)= {'from alm'}; Table15(7,2)={std(kmalm)};
Table15(8,1)= {'% differ'}; Table15(8,2)={abs(std(km_ag)-std(kmalm))/std(km_ag)*100};


% Table 16 in Den Haan's (2008) comparison paper

km_error=abs((km_ag-kmalm)./km_ag)*100; % percentage error between km 
   % computed from the capital distribution (km_ag) and km computed from 
   % the ALM (kmalm)
Table16(1,1)= {'average'};
Table16(1,2)= {' '};
Table16(2,1)= {'total'};  Table16(2,2)={mean(km_error)};
Table16(3,1)= {'maximal'};
Table16(3,2)= {' '};
Table16(4,1)= {'total'};  Table16(4,2)={max(km_error)};

%__________________________________________________________________________
%
% 9. EULER EQUATION ERRORS ON SIMULATED TIME PATH
%__________________________________________________________________________

c_eeq=zeros(T-1,1); % individual consumption computed from the Euler 
                    % equation using true individual and aggregate capital 

   for t=1:T-1
       
      t

      curr_state=(agshock(t)-1)*2+idshock(t,1); % current state
      
      % Future (at t+1) interest rate 
      
      irate_b=alpha*a(1)*(km_ag(t+1)/er_b/l_bar)^(alpha-1);
      irate_g=alpha*a(2)*(km_ag(t+1)/er_g/l_bar)^(alpha-1);
      
      % Future (at t+1) wage
      
      wage_b=(1-alpha)*a(1)*(km_ag(t+1)/er_b/l_bar)^(alpha);
      wage_g=(1-alpha)*a(2)*(km_ag(t+1)/er_g/l_bar)^(alpha);
      
      % Future (at t+2) value of capital (i.e., k'') in four possible
      % states (temporary variable)
      
      k_hat1_bu=interpn(k,km,a2,epsilon2,kprime,k_ind(t+1,1),km_ag(t+1),1,1,'linear');
      k_hat1_be=interpn(k,km,a2,epsilon2,kprime,k_ind(t+1,1),km_ag(t+1),1,2,'linear');
      k_hat1_gu=interpn(k,km,a2,epsilon2,kprime,k_ind(t+1,1),km_ag(t+1),2,1,'linear');      
      k_hat1_ge=interpn(k,km,a2,epsilon2,kprime,k_ind(t+1,1),km_ag(t+1),2,2,'linear');      
     
      % Current consumption computed from the Euler equation
      
      curr_state=(agshock(t)-1)*2+idshock(t,1); % supplementary variable for 
                                                % computing the transition 
                                                % probabilities from "prob" 
      
      c_eeq(t,1)=(beta*(((1-delta+irate_b)*k_ind(t+1,1)+wage_b*l_bar*epsilon_u+mu*wage_b*(1-epsilon_u)-mu*wage_b*ur_b/(1-ur_b)*epsilon_u-k_hat1_bu)^(-gamma)*(1-delta+irate_b)*prob(curr_state,1)...
         +((1-delta+irate_b)*k_ind(t+1,1)+wage_b*l_bar*epsilon_e+mu*wage_b*(1-epsilon_e)-mu*wage_b*ur(1)/(1-ur(1))*epsilon_e-k_hat1_be)^(-gamma)*(1-delta+irate_b)*prob(curr_state,2)...
         +((1-delta+irate_g)*k_ind(t+1,1)+wage_g*l_bar*epsilon_u+mu*wage_g*(1-epsilon_u)-mu*wage_g*ur(2)/(1-ur(2))*epsilon_u-k_hat1_gu)^(-gamma)*(1-delta+irate_g)*prob(curr_state,3)...
         +((1-delta+irate_g)*k_ind(t+1,1)+wage_g*l_bar*epsilon_e+mu*wage_g*(1-epsilon_e)-mu*wage_g*ur(2)/(1-ur(2))*epsilon_e-k_hat1_ge)^(-gamma)*(1-delta+irate_g)*prob(curr_state,4)))^(-1/gamma);
          
      % Future (at j+1) capital (i.e., k') restored from the budget constraint 
      
      k_eeq=k_ind(t,1)*(1-delta+irate(t))+(idshock(t,1)-1).*l_bar*wage(t)+mu*(2-idshock(t,1)).*wage(t)-...
      ur(agshock(t))/(1-ur(agshock(t)))*mu*(idshock(t,1)-1).*wage(t)-c_eeq(t,1);
  
      % Check that k_eeq is not negative 
      
      k_eeq=k_eeq*(k_eeq>0); % set k_eeq=0 if it is negative
  
      % Re-compute c_eeq whenever k_eeq was initially negative 
      
      c_eeq(t,1)=k_ind(t,1)*(1-delta+irate(t))+(idshock(t,1)-1).*l_bar*wage(t)+mu*(2-idshock(t,1)).*wage(t)-...
      ur(agshock(t))/(1-ur(agshock(t)))*mu*(idshock(t,1)-1).*wage(t)-k_eeq;

   end
   
% Table 1 in Maliar, Maliar and Valli's (2008)

c_error=abs((c_ind(1:T-1,1)-c_eeq)./c_eeq)*100; % percentage difference 
                                                % between c_eeq and c_ind
Table1(1,1)= {'average'};Table1(1,2)={mean(c_error')};
Table1(2,1)= {'maximal'};Table1(2,2)={max(c_error')};
%__________________________________________________________________________
%
% 10. RESULTS
%__________________________________________________________________________

% Results corresponding to Tables 9,11,12,14,15,16 in Den Haan (2008) and
% Table 1 in Maliar, Maliar and Valli (2008)

disp('Table 9. Means and standard deviations of cross-sectional moments'); Table9
disp('Table 11. Means of the 5th and 10th percentile'); Table11
disp('Table 12. Properties of individual policy rules'); Table12
disp('Table 14. Percentage errors from dynamic Euler accuracy test'); Table14
disp('Table 15. Moments of Kt in panel and according to aggregate law of motion'); Table15
disp('Table 16. Accuracy aggregate policy rule'); Table16
disp('Table 1.  Euler equation errors'); Table1
%__________________________________________________________________________
%
% SAVE RESULTS IN FILE "Results_of_test"
%__________________________________________________________________________
%
save Results_of_test;