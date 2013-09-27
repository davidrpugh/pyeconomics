% The program for the article "Solving the incomplete markets model with
% aggregate uncertainty using the Krusell-Smith algorithm" from the special 
% JEDC issue edited by Den Haan, Judd and Juillard (2008)  
%
% Written by Lilia Maliar, Serguei Maliar and Fernando Valli (2008)

function [kmts,kcross]  = AGGREGATE_NS(l_bar,alpha,prob,ur_b,ur_g,T,J,kvalues_min,kvalues_max,ngridk,ngridkm,nstates_ag,nstates_id,idshock,agshock,km_max,km_min,kprime,km,k,epsilon2,ndiscard,k_min,k_max,kcross,a,a2);

ur=zeros(nstates_ag,1); ur(1)=ur_b; ur(2)=ur_g; % vector of unemployment 
% rates in two aggregate states

er_b=(1-ur_b);   % employment rate in a bad aggregate state
er_g=(1-ur_g);   % employment rate in a good aggregate state

kmts=zeros(T,1);         % a time series of the mean of capital distribution 

% Time series of aggregate variables

prod_ag=(agshock==1)*a(1)+(agshock==2)*a(2);   % aggregate productivity
labor_ag=(agshock==1)*er_b+(agshock==2)*er_g;  % aggregate employment (L in 
                                               % the paper)

irate=zeros(T-1,1); % interest rate
wage=zeros(T-1,1);  % wage

% Beginning-of-period capital distributions

kdistu=zeros(T,J); % beginning-of-period capital distributions in all 
                   % periods for the unemployed
kdiste=zeros(T,J); % for the employed

% Initial-period capital distributions

% kdist has dimensionality J*2 with J=1001; in the last point J, the probability 
% is chosen to normalize the sum of probabilities to 1

kdistu(1,1:J)=[kcross(1,1:J-1) 1-sum(kcross(1,1:J-1))]; % a raw vector of 
                       % the initial capital distribution for the unemployed
kdiste(1,1:J)=[kcross(2,1:J-1) 1-sum(kcross(2,1:J-1))]; % for the employed

% Grid for the capital distribution

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
   
   % Aggregate capital
   
   kmts(t)=kdistu(t,:)*kvalues*(1-labor_ag(t))+kdiste(t,:)*kvalues*labor_ag(t);
      % aggregate capital=capital of the unemployed + capital of the employed;
      % (1-labor_ag(t)) is the share of unemployed people in the economy
   kmts(t)=kmts(t)*(kmts(t)>=km_min)*(kmts(t)<=km_max)+km_min*(kmts(t)<km_min)+km_max*(kmts(t)>km_max); % restrict kmts to be within [km_min, km_max]

   % Prices
   
   irate(t)=alpha*prod_ag(t)*(kmts(t)/labor_ag(t)/l_bar)^(alpha-1);
      % interest rate 
   wage(t)=(1-alpha)*prod_ag(t)*(kmts(t)/labor_ag(t)/l_bar)^alpha;
      % wage
  %________________________________________________________________________
  %
  % Individual capital function, k'
  %________________________________________________________________________
   
  kprimet(:,1)=interpn(k,km,kprime(:,:,agshock(t),1),kvalues,kmts(t)*ones(J,1),'linear'); 
      % interpolate the capital function k' (computed in "MAIN") of the 
      % unemployed agent in kvalues for the given agg. capital km_ag(t)
      
  kprimet(:,2)=interpn(k,km,kprime(:,:,agshock(t),2),kvalues,kmts(t)*ones(J,1),'linear');
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
% Aggregate capital in T
%__________________________________________________________________________

kmts(T)=kdistu(T,:)*kvalues*(1-labor_ag(T))+kdiste(T,:)*kvalues*labor_ag(T);
                          % aggregate capital in T
kmts(T)=kmts(T)*(kmts(T)>=km_min)*(kmts(T)<=km_max)+km_min*(kmts(T)<km_min)+km_max*(kmts(T)>km_max); % restrict kmts to be within [km_min, km_max]
%__________________________________________________________________________
%
% Terminal distribution of capital
%__________________________________________________________________________

kcross(1,:)=kdistu(T,:); % distribution for unemployed
kcross(2,:)=kdiste(T,:); % distribution for employed