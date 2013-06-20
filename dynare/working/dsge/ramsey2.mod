////////// Ramsey model (i.e., RBC model with inelastic labor supply!) //////////

////////// Declare parameters //////////

///// Household parameters /////
parameters beta, theta;

// discount factor
beta = 0.9896;

// coefficient of relative risk aversion
theta = 2.0;

///// Firm parameters /////
parameters alpha, delta;

// capital's share of income/output
alpha = 0.4;

// depreciation rate
delta = 0.0196;

///// Growth rates /////
parameters g, n;

// technology growth (average quarterly growth of Solow residual for U.S. from PWT)
g = 0.0031875459362321018;

// population growth (average quarterly growth of labor force for U.S. from PWT)
n = 0.0037722357644148748;

///// Stochastic processes /////
parameters rho, sigma;

// persistence of productivity shock
rho = 0.95;

// standard deviation of productivity shock
sigma = 0.01;

////////// Declare variables //////////

///// Endogenous variables /////

// y: output per effective worker
// k: capital per effective worker
// c: consumption per effective worker
// i: investement per effective worker
// r: net interest rate 
// w: real wage
// z: productivity
// check1: zero profit condition
var y, k, c, i, r, w, z, check1; 

///// Exogenous variables /////

// productivity shock
varexo eps;

////////// The model //////////
model;
// aggregate resource constraint
y = c + i;

// production function 
y = ((z(-1) * k(-1)) / ((1 + g) * (1 + n) * z))^alpha;

// net interest rate
r = alpha * ((z(-1) * k(-1)) / ((1 + g) * (1 + n) * z))^(alpha - 1) - delta;

// real wage 
w = (1 - alpha) * y;

// equation of motion for capital
k = (1 - delta) * ((z(-1) * k(-1)) / ((1 + g) * (1 + n) * z)) + i;

// consumption Euler equation
1 = beta * (1 + g)^(-theta) * ((c(+1) * z(+1)) / (c * z))^(-theta) * (1 + r(+1));

// productivity process
z = z(-1)^rho * exp(eps);

// check that zero profit condition holds
check1 = y - w - (r + delta) * ((z(-1) * k(-1)) / ((1 + g) * (1 + n) * z));
end;

////////// Find the steady state //////////

kss = 1.0;
yss = (kss / ((1 + g) * (1 + n)))^alpha;
rss = alpha * (kss / ((1 + g) * (1 + n)))^(alpha - 1) - delta;
wss = (1 - alpha) * yss;
iss = (1 - delta) * (kss / ((1 + g) * (1 + n))) - kss;
css = yss - iss;
zss = 1.0;
check1ss = 0.0;

initval;
k = 1.0;
y = yss;
r = rss;
w = wss;
i = iss;
z = zss;
check1 = check1ss;
end;

steady;

////////// Stability analysis //////////
check;

////////// Shocks //////////
shocks;
var eps; stderr sigma;
end;

////////// Simulation //////////
stoch_simul(order=1, irf=40);