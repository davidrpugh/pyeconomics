% RBC Model with variable capacity utilization.
%
% Depreciation function: delta*u^omega
% Calibrate omega and delta to generate k/y=10 and u_ss=0.9
% Calibrate variance shock to generate same volatility output 
% than basic RBC. 
% We can reduce standard deviation of shock by 39%.
%
% Jesus Fernandez-Villaverde
% Philadelphia, March 3, 2005

%----------------------------------------------------------------
% 0. Housekeeping
%----------------------------------------------------------------

close all;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var y c k i u l y_l z;
varexo e;

// same warning concerning beta and gamma
parameters beta psi alpha delta omega rho sigma;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

alpha   = 0.33;
beta    = 0.99;
psi     = 1.75;
omega   = 1.45;
delta   = 0.0265;
rho     = 0.95;  
sigma   = (0.0043/(1-alpha));

%----------------------------------------------------------------
% 3. Model   
%----------------------------------------------------------------

model; 
  (1/c) = beta*(1/c(+1))*(1+alpha*(u(+1)^alpha)*(k^(alpha-1))*(exp(z(+1))*l(+1))^(1-alpha)-delta*(u(+1)^omega));
  psi*c/(1-l) = (1-alpha)*((u*k(-1))^alpha)*(exp(z)^(1-alpha))*(l^(-alpha));
  y/k(-1) = omega*delta*(u^omega)/alpha;
  c+i = y;
  y = ((u*k(-1))^alpha)*(exp(z)*l)^(1-alpha);
  y_l = y/l;
  i = k-(1-(delta*(u^omega)))*k(-1);
  z = rho*z(-1)+e;
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

initval;
  k = 10;
  c = 0.76;
  l = 0.31;
  y = 1;
  i = 0.24;
  u = 0.8;
  z = 0; 
  e = 0;
end;

shocks;
var e = sigma^2;
end;

steady;

stoch_simul(hp_filter = 1600, order = 1);

%----------------------------------------------------------------
% 5. Some Results
%----------------------------------------------------------------

statistic1 = 100*sqrt(diag(oo_.var(1:7,1:7)))./oo_.mean(1:7);
dyntable('Relative standard deviations in %',strvcat('VARIABLE','REL. S.D.'),M_.endo_names(1:7,:),statistic1,10,8,4);
