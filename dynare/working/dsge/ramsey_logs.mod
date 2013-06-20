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

// lk: log-capital per effective worker
// lc: log-consumption per effective worker
// lz: log-productivity process
var lk, lc, lz; 

///// Exogenous variables /////

// productivity shock
varexo eps;

////////// The model //////////
model;

    // equation of motion for capital
    exp(lk) = (1 - delta) * (exp(lz(-1) - lz + lk(-1)) / ((1 + g) * (1 + n))) +
    	      (exp(lz(-1) - lz + lk(-1)) / ((1 + g) * (1 + n)))^alpha - exp(lc);

    // consumption Euler equation
    1 = beta * (1 + g)^(-theta) * exp(-theta * (lc(+1) + lz(+1) - lc - lz)) *
        (1 + alpha * (exp(lz - lz(+1) + lk) / ((1 + g) * (1 + n)))^(alpha - 1) - delta);

    // productivity process
    exp(lz) = exp(rho * lz(-1) + eps);

end;

////////// Find the steady state //////////

// Use fsolve
initval;
lk = 1.0;
lc = 0.5;
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