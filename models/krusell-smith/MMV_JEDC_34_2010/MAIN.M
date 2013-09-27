% The program for the article "Solving the incomplete markets model with
% aggregate uncertainty using the Krusell-Smith algorithm" from the special 
% JEDC issue edited by Den Haan, Judd and Juillard (2008)  
%
% Written by Lilia Maliar, Serguei Maliar and Fernando Valli (2008)
% 
% The program includes the following files: 
%
% 1. "MAIN.m" (computes a solution and stores the results in "Solution")
% 2. "SHOCKS.m"      (a subroutine of MAIN.m; generates the shocks)
% 3. "INDIVIDUAL.m"  (a subroutine of MAIN.m; computes a solution to the 
%                     individual problem)
% 4. "AGGREGATE_ST.m"   (a subroutine of MAIN.m; performs the stochastic 
%                        simulation)
% 5. "AGGREGATE_NS.m"   (a subroutine of MAIN.m; performs the non-stochastic 
%                        simulation)
% 6. "Inputs_for_test" (contains initial distribution of capital and 
%    10,000-period realizations of aggregate shock and idiosyncratic shock 
%    for one agent provided by Den Haan, Judd and Juillard, 2008) 
% 7. "TEST.m" (should be run after "MAIN.m"; it uses "Inputs_for_test" and
%    "Solution_to_model" for computing the statistics reported in Den Haan's
%    2008, comparison article)  
%
% See the web page of the authors for the updated versions of the program 
% __________________________________________________________________________
clc;
clear all;
%__________________________________________________________________________
%
% Simulation method (stochastic and non-stochastic) 
%__________________________________________________________________________

Method=2; % Method=1 stands for the stochastic simulation; for the non-
          % stochastic simulation, set Method=2
          
% Parameters for Method 1 only

N=10000;          % number of agents for the stochastic simulation

% Parameters for Method 2 only

J=1000;           % number of grid points for the non-stochastic simulation
kvalues_min=0;    % minimum grid value for the non-stochastic simulation
kvalues_max=100;  % maximum grid value for the non-stochastic simulation
%__________________________________________________________________________
%
% Parameters
%__________________________________________________________________________

beta=0.99;       % discount factor
gamma=1;         % utility-function parameter
alpha=0.36;      % share of capital in the production function
delta=0.025;     % depreciation rate
delta_a=0.01;    % (1-delta_a) is the productivity level in a bad state, 
                 % and (1+delta_a) is the productivity level in a good state
mu = 0.15;       % unemployment benefits as a share of wage
l_bar=1/0.9;     % time endowment; normalizes labor supply to 1 in a bad state
T=1100;          % simulation length
ndiscard=100;    % number of periods to discard

nstates_id=2;    % number of states for the idiosyncratic shock
nstates_ag=2;    % number of states for the aggregate shock

epsilon_u=0;     % idiosyncratic shock if the agent is unemployed
epsilon_e=1;     % idiosyncratic shock if the agent is employed

ur_b=0.1;        % unemployment rate in a bad aggregate state
er_b=(1-ur_b);   % employment rate in a bad aggregate state
ur_g=0.04;       % unemployment rate in a good aggregate state
er_g=(1-ur_g);   % employment rate in a good aggregate state

% Matrix of transition probabilities in Den Haan, Judd, Juillard (2008)

prob=[0.525 0.35 0.03125 0.09375  
   0.038889 0.836111 0.002083 0.122917
   0.09375 0.03125 0.291667 0.583333
   0.009115 0.115885 0.024306 0.850694];

kss=((1/beta-(1-delta))/alpha)^(1/(alpha-1)); % steady-state capital in a 
        % deterministic model with employment rate of 0.9 (i.e., l_bar*L=1 
        % where L is aggregate labor in the paper) 
%__________________________________________________________________________
%
% Generation of shocks
%__________________________________________________________________________

% Method 1 uses both idiosyncartic and aggregate shocks, while Method 2 uses
% only aggregate shock

[idshock,agshock]  = SHOCKS(prob,T,N,ur_b);
%__________________________________________________________________________
%
% Grids
%__________________________________________________________________________

% Grid for capital in the individual problem

k_min=0;                   % minimum grid-value of capital
k_max=1000;                % maximum grid-value of capital
ngridk=100;                % number of grid points
x=linspace(0,0.5,ngridk)'; % generate a grid of ngridk points on [0,0.5] 
                           % interval  
y=x.^7/max(x.^7);          % polynomial distribution of grid points, formula 
                           % (7) in the paper
k=k_min+(k_max-k_min)*y;   % transformation of grid points from [0,0.5] 
                           % interval to [k_min,k_max] interval

% Grid for mean of capital

km_min=30;                           % minimum grid-value of the mean of 
                                     % capital distribution, km 
km_max=50;                           % maximum grid value of km
ngridkm=4;                           % number of grid points for km 
km=linspace(km_min,km_max,ngridkm)'; % generate a grid of ngridkm points on 
                                     % [km_min,km_max] interval 
%__________________________________________________________________________
%
% Parameters of idiosyncratic and aggregate shocks
%__________________________________________________________________________

epsilon=zeros(nstates_id,1);  % vector of possible idiosyncratic states
epsilon(1)=epsilon_u; epsilon(2)=epsilon_e; % the unemployed and employed 
                              % states are 0 and 1, respectively 
epsilon2=zeros(nstates_id,1); % vector of possible idiosyncratic states 
epsilon2(1)=1; epsilon2(2)=2; % the unemployed and employed states are 
                              % 1 and 2, respectively  

a=zeros(nstates_ag,1);        % vector of possible aggregate states
a(1)=1-delta_a; a(2)=1+delta_a; % bad and good aggregate states are 1-delta_a 
                              % and 1+delta_a, respectively
a2=zeros(nstates_ag,1);       % vector of possible aggregate states
a2(1)=1; a2(2)=2;             % bad and good aggregate states are 1 and 2, 
                              % respectively
%__________________________________________________________________________
%
% Initial conditions 
%__________________________________________________________________________

kprime=zeros(ngridk,ngridkm,nstates_ag,nstates_id); % next-period individual 
   % capital (k') depends on four state variables: individual k, aggregate k, 
   % aggregate shock, idiosyncratic shock 

% Initial capital function

for i=1:ngridkm
   for j=1:nstates_ag
      for h=1:nstates_id
         kprime(:,i,j,h)=0.9*k;
      end
   end
end

% Initial distribution of capital is chosen so that aggregate capital is
% near the steady state value, kss

if Method==1;    

    % inititial distribution of capital for the stochastic simulation 

    kcross=zeros(1,N)+kss;  % initial capital of all agents is equal to kss

else % i.e. Method 2
    
    % inititial density on the grid for the non-stochastic simulation

    kcross=zeros(2,J);       % density function for the employed and 
                             % unemployed agents, defined in J grid points 
    igrid=(kvalues_max-kvalues_min)/J; % interval between grid points
    jss=round(kss/igrid);    % # of the grid point nearest to kss
    kcross(1:2,jss)=1;       % all density is concentrated in the point jjs;
                             % density is zero in all other grid points
end

% Initial vector of coefficients B of the ALM (in the paper, it is b) 

% (the ALM in a bad state is ln(km')=B(1)+B(2)*ln(km) and the ALM in a good
% state is ln(km')=B(3)+B(4)*ln(km))

B=[0 1 0 1];
%__________________________________________________________________________
%
% Convergence parameters
%__________________________________________________________________________

dif_B=10^10;   % difference between coefficients B of the the ALM on 
               % successive iterations; initially, set to a large number
criter_k=1e-8; % convergence criterion for the individual capital function
criter_B=1e-8; % convergence criterion for the coefficients B in the ALM
update_k=0.7;  % updating parameter for the individual capital function
update_B=0.3;  % updating parameter for the coefficients B in the ALM
%__________________________________________________________________________
%
% SOLVING THE MODEL
%__________________________________________________________________________

iteration=0      % initial iteration  
init_time=clock; % initialize the time clock

i=1; while dif_B>criter_B % perform iterations until the difference between 
                          % coefficients is less or equal than criter_B
                          
[kprime,c]  = INDIVIDUAL(prob,ur_b,ur_g,ngridk,ngridkm,nstates_ag,nstates_id,k,km,er_b,er_g,a,epsilon,l_bar,alpha,delta,gamma,beta,mu,km_max,km_min,kprime,B,criter_k,k_min,k_max,update_k);
            % compututing a solution to the individual problem

if Method==1; 
    [kmts,kcross1]  = AGGREGATE_ST(T,idshock,agshock,km_max,km_min,kprime,km,k,epsilon2,k_min,k_max,kcross,a2);
else % i.e., Method 2 
    [kmts,kcross1]  = AGGREGATE_NS(l_bar,alpha,prob,ur_b,ur_g,T,J,kvalues_min,kvalues_max,ngridk,ngridkm,nstates_ag,nstates_id,idshock,agshock,km_max,km_min,kprime,km,k,epsilon2,ndiscard,k_min,k_max,kcross,a,a2);
end

% Time series for the ALM regression 

ibad=0;           % count how many times the aggregate shock was bad
igood=0;          % count how many times the aggregate shock was good
xbad=0;  ybad=0;  % regression-variables for a bad state
xgood=0; ygood=0; % regression-variables for a good state
for i=ndiscard+1:T-1
   if agshock(i)==1
      ibad=ibad+1;
      xbad(ibad,1)=log(kmts(i));
      ybad(ibad,1)=log(kmts(i+1));
   else
      igood=igood+1;
      xgood(igood,1)=log(kmts(i));
      ygood(igood,1)=log(kmts(i+1));
   end
end

[B1(1:2),s2,s3,s4,s5]=regress(ybad,[ones(ibad,1) xbad]);R2bad=s5(1); 
    % run the OLS regression ln(km')=B(1)+B(2)*ln(km) for a bad agg. state 
    % and compute R^2 (which is the first statistic in s5)
[B1(3:4),s2,s3,s4,s5]=regress(ygood,[ones(igood,1) xgood]);R2good=s5(1);
    % make the OLS regression ln(km')=B(3)+B(4)*ln(km) for a good agg. state 
    % and compute R^2 (which is the first statistic in s5)

dif_B=norm(B-B1) % compute the difference between the initial and obtained 
                 % vector of coefficients

% To ensure that initial capital distribution comes from the ergodic set,
% we use the terminal distribution of the current iteration as initial 
% distribution for a subsequent iteration. When the solution is sufficiently 
% accurate, dif_B<(criter_B*100), we stop such an updating and hold the 
% distribution "kcross" fixed for the rest of iterations. ·

if dif_B>(criter_B*100)
    kcross=kcross1; % the new capital distribution  replaces the old one
end

B=B1*update_B+B*(1-update_B); % update the vector of the ALM coefficients 
                 % according to the rule (9) in the paper
iteration=iteration+1

end

end_time=clock;               % end the time clock
et=etime(end_time,init_time); % compute time in seconds that has elapsed 
                              % between init_time and end_time
disp('Elapsed Time (in seconds):'); et
disp('Iterations');          iteration
format long g; 
disp('R^2 bad aggregate shock:'); R2bad(1)
disp('R^2 good aggregare shock:'); R2good(1)
format; 
%__________________________________________________________________________
%
% FIGURE OF THE AGGREGATE TIME SERIES SOLUTION
%__________________________________________________________________________

kmalm=zeros(T,1);  % represents aggregate capital computed from the ALM
kmalm(1)=kmts(1);  % in the first period km computed from the ALM (kmalm) 
                   % is equal km computed from the cross-sectional capital 
                   % distribution (kmts)
                                   
for t=1:T-1       % compute kmalm for t=2:T
   if agshock(t)==1
      kmalm(t+1)=exp(B(1)+B(2)*log(kmalm(t)));
   else
      kmalm(t+1)=exp(B(3)+B(4)*log(kmalm(t)));
   end
   
end

Tts=1:1:T;
axis([min(Tts) max(Tts) min(kmts)*0.99 max(kmts)*1.01]); axis manual; hold on; 
plot (Tts,kmts(1:T,1),'-',Tts,kmalm(1:T,1),'--'),xlabel('Time'), ylabel('Aggregate capital series'), title('Figure 1. Accuracy of the aggregate law of motion.')
legend('implied by individual policy rule', 'aggregate law of motion')
%__________________________________________________________________________
%
% SAVE RESULTS IN FILE "Solution_to_model"
%__________________________________________________________________________
%
save Solution_to_model;