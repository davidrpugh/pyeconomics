////////// RBC model with inelastic labor supply //////////

/* Implementation of Romer's RBC model without a household labor 
supply decision. The level of technology is driven by the dollowing 
stochastic process:

    A(t) = (1 + g) * (z(t) / z(t-1)) * A(t-1) 

where the technology shocks z(t) follow an AR(1) process in logs:

    ln z(t) = rho * ln(z(t-1)) + eps(t) 

with eps ~ N(0, sigma).  This specification assumes that technology 
fluctuates, persistantly, around an underlying trend growth rate of g.

*/

////////// Declare variables //////////

///// Endogenous variables /////

// y: output per effective worker
// k: capital per effective worker
// c: consumption per effective worker
// i: investment per effective worker
// z: productivity process
// check1: zero profit condition
var y, k, c, i, z, check1; 

///// Exogenous variables /////

// productivity shock
varexo eps;

////////// Declare parameters //////////
parameters beta, theta, alpha, delta, n, g, rho, sigma;

// discount factor
beta = 0.9896;

// coefficient of relative risk aversion
theta = 2.0;

// capital's share of income/output
alpha = 0.4;

// depreciation rate
delta = 0.0196;

// technology growth (average quarterly growth of Solow residual for U.S. from PWT)
g = 0.0031875459362321018;

// population growth (average quarterly growth of labor force for U.S. from PWT)
n = 0.0037722357644148748;

// persistence of productivity shock
rho = 0.95;

// standard deviation of productivity shock
sigma = 0.01;

////////// The model //////////
model;

// production
y = ((z(-1) * k(-1)) / ((1 + g) * (1 + n) * z))^alpha;

// aggregate resource constraint
y = c + i;

// equation of motion for capital
k = (1 - delta) * ((z(-1) * k(-1)) / ((1 + g) * (1 + n) * z)) + i;

// consumption Euler equation
1 = beta * (1 + g)^(-theta) * ((c(+1) * z(+1)) / (c * z))^(-theta) *
    (1 + alpha * ((z * k) / ((1 + g) * (1 + n) * z(+1)))^(alpha - 1) - delta);

// productivity process
z = z(-1)^rho * exp(eps);

// check that zero profit condition holds
check1 = y -
         (1 - alpha) * ((z(-1) * k(-1)) / ((1 + g) * (1 + n) * z))^alpha -
         (alpha * ((z(-1) * k(-1)) / ((1 + g) * (1 + n) * z))^(alpha - 1) - delta + delta) *
         ((z(-1) * k(-1)) / ((1 + g) * (1 + n) * z));
end;

////////// Find the steady state //////////

// Use known analyic solution for steady state values
initval;
val = alpha * beta * (1 + g)^(-theta) / (1 - beta * (1 + g)^(-theta) * (1 - delta)); // local var
k = (1 + g) * (1 + n) * val^(1 / (1 - alpha));
y = (k / ((1 + g) * (1 + n)))^alpha;
z = 1.0;
c = (1 - delta) * (k / ((1 + g) * (1 + n))) + (k / ((1 + g) * (1 + n)))^alpha - k;
i = y - c; 
check1 = 0.0;
end;

////////// Variance covariance matrix //////////
vcov = [0.007];