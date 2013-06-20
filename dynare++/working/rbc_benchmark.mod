////////// Declare variables //////////

/* ///// Endogenous variables /////
   
   k: capital
   c: consumption
   z: productivity
   y: output
   i: investment
   w: real wage
   r: net interest rate
   l: labor supply
   check1: zero profit condition
*/
var k, c, z, y, i, w, r, l; //check1 

///// Exogenous variables /////

// eps: productivity shock 
varexo eps;

////////// Declare parameters //////////
parameters beta, theta, omega, alpha, delta, rho, sigma, kss, lss;

// discount factor
beta  = 0.9896;

// coefficient of relative risk aversion
theta = 2.0;

// weights utility from consumption and leisure
omega = 0.357;

// capital's share of income
alpha = 0.40;

// depreciation rate of capital
delta = 0.0196;

// persistence of productivity process
rho   = 0.95;

// standard deviation of productivity shocks
sigma = 0.007;

// steady state capital (computed for theta=1 and delta=1)
kss =  ((omega * (1 - alpha)) / (omega * (1 - alpha) + (1 - omega) * (1 - alpha * beta))) * 
       (alpha / ((1 / beta) - 1 + delta))^(1 / (1 - alpha));

// steady state labor supply (computed for theta=1 and delta=1)
lss = (omega * (1 - alpha)) / (omega * (1 - alpha) + (1 - omega) * (1 - alpha * beta));

////////// Model equations //////////

model;
// production
y = exp(z) * k(-1)^alpha * l^(1 - alpha);

// real wage
w = (1 - alpha) * (y / l);

// net marginal product of capital
r = alpha * (y / k(-1)) - delta;

// resource constraint
y = c + i;

// consumption Euler equation
c^(omega - 1) * (1 - l)^(1 - omega) * (c^omega * (1 - l)^(1 - omega))^(-theta) = 
beta * c(+1)^(omega - 1) * (1 - l(+1))^(1 - omega) * (c(+1)^omega * (1 - l(+1))^(1 - omega))^(-theta) * (1 + r(+1));
 
// intra-temporal consumption/labor trade-off
omega * (1 - l) * w = (1 - omega) * c;

// equation of motion for capital
k = (1 - delta) * k(-1) + i;

// productivity process
z = rho * z(-1) + eps;

// check that zero profit condition holds
//check1 = y - w * l - (r + delta) * k(-1);
end;

////////// Initial values for computing steady state //////////

initval;
k = kss;
l = lss;
c = kss^alpha * lss^(1 - alpha) - delta * kss;
z = 0.0;
y = kss^alpha * lss^(1 - alpha);
i = delta * kss;
w = (1 - alpha) *kss^alpha * lss^(-alpha);
r = (1 / beta) - 1;
check1 = 0.0;
end;

////////// Variance covariance matrix //////////
vcov = [0.007];
