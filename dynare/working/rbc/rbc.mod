% Basic RBC Model 
%
% Jesus Fernandez-Villaverde
% Philadelphia, March 3, 2005

%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------

close all;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var y c k i l y_l z;
varexo e;

parameters beta psi delta alpha rho;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

alpha   = 0.33;
beta    = 0.99;
delta   = 0.023;
psi     = 1.75;
rho     = 0.95;  
sigma   = (0.007/(1-alpha));

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model; 
  (1/c) = beta*(1/c(+1))*(1+alpha*(k^(alpha-1))*(exp(z(+1))*l(+1))^(1-alpha)-delta);
  psi*c/(1-l) = (1-alpha)*(k(-1)^alpha)*(exp(z)^(1-alpha))*(l^(-alpha));
  c+i = y;
  y = (k(-1)^alpha)*(exp(z)*l)^(1-alpha);
  i = k-(1-delta)*k(-1);
  y_l = y/l;
  z = rho*z(-1)+e;
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

initval;
  k = 9;
  c = 0.76;
  l = 0.3;
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

statistic1 = 100*sqrt(diag(oo_.var(1:6,1:6)))./oo_.mean(1:6);
dyntable('Relative standard deviations in %',strvcat('VARIABLE','REL. S.D.'),M_.endo_names(1:6,:),statistic1,10,8,4);
