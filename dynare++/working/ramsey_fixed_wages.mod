////////// Declare variables //////////

///// Endogenous variables /////

/* List of variables:
   
   k: capital
   c: consumption
   l: labor 
   z: productivity
   y: output
   i: investment
   r: net interest rate
   check1: zero profit condition

*/
var k, c, l, z, y, i, r; //check1;

///// Exogenous variables /////

// eps: productivity shock 
varexo eps;

////////// Declare parameters //////////
parameters beta, theta, alpha, delta, w, rho, sigma;

// discount factor
beta  = 0.9896;

// coefficient of relative risk aversion
theta = 2.0;

// capital's share of income
alpha = 0.40;

// depreciation rate of capital
delta = 0.0196;

// real wage is assumed fixed!
w = (1 - alpha) * ((alpha * beta) / (1 - beta * (1 - delta)))^(alpha / (1 - alpha)); 

// persistence of productivity process
rho   = 0.95;

// standard deviation of productivity shocks
sigma = 0.007;

////////// Model equations //////////

model;
// production
y = exp(z) * k(-1)^alpha * l^(1 - alpha);

// With fixed wages, equilibrium is determined by firm's labor demand
w = (1 - alpha) * (y / l);

// net marginal product of capital
r = alpha * exp(z) * k(-1)^(alpha - 1) * l^(1 - alpha) - delta;

// resource constraint
y = c + i;

// consumption Euler equation
c^(-theta) = beta * c(+1)^(-theta) * (1 + r(+1));

// equation of motion for capital
k = (1 - delta) * k(-1) + i;

// productivity process
z = rho * z(-1) + eps;

// check zero profit condition holds
//check1 = y - w - (r + delta) * k(-1);
end;

////////// Initial values for computing steady state //////////
initval;
k = (alpha * beta / (1 - beta * (1 - delta)))^(1 / (1 - alpha));
c = 1 - alpha * beta * delta / (1 - beta * (1 - delta));
l = 1.0;
z = 0.0;
y = (alpha * beta / (1 - beta * (1 - delta)))^(alpha / (1 - alpha));
i = (alpha * beta / (1 - beta * (1 - delta)))^(alpha / (1 - alpha)) + 
    (alpha * beta * delta / (1 - beta * (1 - delta))) - 1;
r = alpha * (alpha * beta / (1 - beta * (1 - delta)))^alpha - delta;
//check1 = 0.0;
end;

////////// Variance covariance matrix //////////
vcov = [0.007];
