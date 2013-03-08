% RBC Model with investment specific technological change,
% capital capacity utilization, and indivisible labor. 
%
% Jesus Fernandez-Villaverde
% Philadelphia, March 4, 2005

%----------------------------------------------------------------
% 0. Housekeeping
%----------------------------------------------------------------

close all;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var y c k i l y_l u z1 z2;
varexo e1 e2;

parameters beta psi delta omega alpha rho1 rho2 sigma1 sigma2;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

alpha   = 0.33;
beta    = 0.99;
psi     = 1.75;
omega   = 1.45;
delta   = 0.0265;
rho1    = 0.95; 
rho2    = 0.64; 
sigma1  = (0.0022/(1-alpha));
sigma2  = 0.0022;

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model; 
  (1/c) = beta*(1/c(+1))*(alpha*(u(+1)^alpha)*(k^(alpha-1))*(exp(z1(+1))*l(+1))^(1-alpha)+(1-delta*(u(+1)^omega))*exp(z2-z2(+1)));
  psi*c = (1-alpha)*((u*k(-1))^alpha)*(exp(z1)^(1-alpha))*(l^(-alpha));
  y/k(-1) = omega*delta*(u^omega)/alpha;
  c+i = y;
  y = ((u*k(-1))^alpha)*(exp(z1)*l)^(1-alpha);
  y_l = y/l;
  k = (1-(delta*(u^omega)))*k(-1)+exp(z2)*i;  
  z1 = rho1*z1(-1)+e1;
  z2 = rho2*z2(-1)+e2;
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
  z1 = 0; 
  z2 = 0;
end;

shocks;
var e1 = sigma1^2;
var e2 = sigma2^2;
end;

steady;

stoch_simul(hp_filter = 1600, order = 1);

%----------------------------------------------------------------
% 5. Some Results
%----------------------------------------------------------------

statistic1 = 100*sqrt(diag(oo_.var(1:6,1:6)))./oo_.mean(1:6);

dyntable('Relative standard deviations in %',strvcat('VARIABLE','REL. S.D.'),M_.endo_names(1:6,:),statistic1,10,8,4);
