////////// Declare variables //////////

///// Endogenous variables /////

/* List of variables:
   
   k: capital
   c: consumption
   z: productivity
   y: output
   i: investment
   w: real wage
   r: net interest rate

*/
var k, c, z, y, i, w, r;

///// Exogenous variables /////

// eps: productivity shock 
varexo eps;

////////// Declare parameters //////////
parameters beta, theta, alpha, delta, rho, sigma;

// discount factor
beta  = 0.9896;

// coefficient of relative risk aversion
theta = 2.0;

// capital's share of income
alpha = 0.40;

// depreciation rate of capital
delta = 0.0196;

// persistence of productivity process
rho   = 0.95;

// standard deviation of productivity shocks
sigma = 0.007;

////////// Model equations //////////

model;
// production
y = exp(z) * k(-1)^alpha;

// real wage
w = (1 - alpha) * y;

// net marginal product of capital
r = alpha * exp(z) * k(-1)^(alpha - 1) - delta;

// resource constraint
y = c + i;

// consumption Euler equation
c^(-theta) = beta * c(+1)^(-theta) * (1 + r(+1));

// equation of motion for capital
k = (1 - delta) * k(-1) + i;

// productivity process
z = rho * z(-1) + eps;
end;

////////// Initial values for computing steady state //////////
initval;
k = (alpha * beta / (1 - beta * (1 - delta)))^(1 / (1 - alpha));
c = 1 - alpha * beta * delta / (1 - beta * (1 - delta));
z = 0.0;
y = (alpha * beta / (1 - beta * (1 - delta)))^(alpha / (1 - alpha));
i = (alpha * beta / (1 - beta * (1 - delta)))^(alpha / (1 - alpha)) + 
    (alpha * beta * delta / (1 - beta * (1 - delta))) - 1;
w = (1 - alpha) * (alpha * beta / (1 - beta * (1 - delta)))^(alpha / (1 - alpha));
r = alpha * (alpha * beta / (1 - beta * (1 - delta)))^alpha - delta;
end;

////////// Variance covariance matrix //////////
vcov = [0.007];
