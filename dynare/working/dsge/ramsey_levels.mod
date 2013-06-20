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

// k: capital per effective worker
// c: consumption per effective worker
// z: productivity process
// check1: zero profit condition
var k, c, z, check1; 

///// Exogenous variables /////

// productivity shock
varexo eps;

////////// The model //////////
model;

// equation of motion for capital
k = (1 - delta) * ((z(-1) * k(-1)) / ((1 + g) * (1 + n) * z)) +
    ((z(-1) * k(-1)) / ((1 + g) * (1 + n) * z))^alpha - c;

// consumption Euler equation
1 = beta * (1 + g)^(-theta) * ((c(+1) * z(+1)) / (c * z))^(-theta) *
    (1 + alpha * ((z * k) / ((1 + g) * (1 + n) * z(+1)))^(alpha - 1) - delta);

// productivity process
z = z(-1)^rho * exp(eps);

// check that zero profit condition holds
check1 = ((z(-1) * k(-1)) / ((1 + g) * (1 + n) * z))^alpha -
         (1 - alpha) * ((z(-1) * k(-1)) / ((1 + g) * (1 + n) * z))^alpha -
         (alpha * ((z(-1) * k(-1)) / ((1 + g) * (1 + n) * z))^(alpha - 1) - delta + delta) *
         ((z(-1) * k(-1)) / ((1 + g) * (1 + n) * z));
end;

////////// Find the steady state //////////

// Use known analyic solution for steady state values
steady_state_model;

// capital per effective person
val = alpha * beta * (1 + g)^(-theta) / (1 - beta * (1 + g)^(-theta) * (1 - delta)); // local var
k = (1 + g) * (1 + n) * val^(1 / (1 - alpha));

// productivity
z = 1.0;

// consumption per effective person
c = (1 - delta) * (k / ((1 + g) * (1 + n))) + (k / ((1 + g) * (1 + n)))^alpha - k;
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