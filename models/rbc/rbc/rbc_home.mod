% RBC Model with Household Production. 
%
% Benhabib, Rogerson, and Wright Model
%
% I borrow their calibration.
%
% Jesus Fernandez-Villaverde
% Philadelphia, March 4th, 2005.

%----------------------------------------------------------------
% 0. Housekeeping
%----------------------------------------------------------------

close all;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var c c_m c_h k k_m k_h l_h l_m i y_m z_m z_h;
varexo e_m e_h;

parameters beta psi alpha theta a eta delta rho_1 rho_2 sigmae_m sigmae_h gamma;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

% Parameters taken directly from BRW
beta  = 0.99;         
alpha = 0.36;
theta = 0.8;         
eta   = 0.08;
delta = 0.025;
rho_1 = 0.95;
rho_2 = 0.95;
sigmae_m = 0.007;
sigmae_h = 0.007;
gamma = 0.66;

% Parameters taken to match BRW observations for hours
a     = 0.33707;
psi   = 0.58756;

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model; 
  c = a*(c_m^theta)+(1-a)*(c_h^theta);
  a*(c_m^(theta-1))/c = beta*(a*(c_m(+1)^(theta-1))/c(+1))*(1+alpha*exp(z_m(+1))*(k_m(+1)^(alpha-1))*(l_m(+1)^(1-alpha))-delta);
  psi/(1-l_m-l_h) = (a*(c_m^(theta-1))/c)*(1-alpha)*exp(z_m)*(k_m^alpha)*(l_m^(-alpha));
  ((1-alpha)*a)/((1-eta)*(1-a))*(c_m^(theta-1))*exp(z_m)*(k_m^alpha)*(l_m^(-alpha)) = (c_h^theta)/l_h;
  ((alpha*a)/(eta*(1-a)))*(c_m^(theta-1))*exp(z_m)*(k_m^(alpha-1))*(l_m^(1-alpha)) = (c_h^theta)/k_h;
  c_m+k = exp(z_m)*(k_m^alpha)*(l_m^(1-alpha))+(1-delta)*k(-1);
  y_m = exp(z_m)*(k_m^alpha)*(l_m^(1-alpha));
  c_h = exp(z_h)*(k_h^eta)*(l_h^(1-eta));
  i = k-(1-delta)*k(-1);
  k(-1) = k_m+k_h;
  z_m = rho_1*z_m(-1)+e_m;
  z_h = rho_2*z_h(-1)+e_h;
end;

%----------------------------------------------------------------
% 4. Computation 
%----------------------------------------------------------------

initval;
  c    = 0.5;
  c_m  = 0.5;
  c_h  = 0.5;
  k    = 12;
  k_m  = 10;
  k_h  = 2;
  i    = 1.2;
  l_m  = 0.33;
  l_h  = 0.28;
  y_m  = 1;
  z_m  = 0;
  z_h  = 0;
  e_m  = 0;
  e_h  = 0;
end;

shocks;
  var e_m = sigmae_m^2;
  var e_h = sigmae_h^2;
  var e_m,e_h = gamma*sigmae_m*sigmae_h;
end;

steady;

stoch_simul(hp_filter = 1600, order = 1);

%----------------------------------------------------------------
% 5. Some Results
%----------------------------------------------------------------

statistic1 = 100*sqrt(diag(oo_.var(1:10,1:10)))./oo_.mean(1:10);
dyntable('Relative standard deviations in %',strvcat('VARIABLE','REL. S.D.'),M_.endo_names(1:10,:),statistic1,10,8,4);
