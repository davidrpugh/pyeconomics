
////////// Declare parameters //////////

///// Household parameters /////
parameters b, beta, gamma, eta, sigma, theta, chi;

// b: coefficient of relative risk aversion (real money balances)
b = 10.0;

// beta: discount factor
beta = 0.99;

// gamma: utility weight (real money balances)
gamma = 1.0;

// eta: coefficient of relative risk aversion (labor supply)
eta = 1.0;

// sigma: coefficient of relative risk aversion (consumption)
sigma = 1.0;

// theta: elasticity of demand for individual goods (theta > 1)
theta = 1.1;

// chi: utility weight (labor supply)
chi = 1.0;

///// Firm parameters /////
parameters omega nu;

// omega: fraction of firms that can NOT change prices
omega = 0.8;

// nu: government production subsidy
nu = 1.0;

///// Policy rule parameters /////

////// Parameters for stochastic processes /////

// rho:       persistence of technology shock
// sigmaA:    standard deviation of technology shock
parameters rho sigmaA;

rho    = 0.95;
sigmaA = 0.007;
 
////////// Declare variables //////////

///// Endogenous variables /////

// C:     Aggregate consumption (flow)
// Y:     Output (flow)
// N:     Aggregate labor (flow)
// M:     Money supply (stock)
// P:     Nominal price level
// pStar: Optimal price set by firms that can change prices
// W:     Nominal wage
// F:     Denominator of RHS of firm's price equation
// K:     Numerator of RHS of firm's price equation
var MC, Y, C, N, M, P, pStar, W, F, K;

///// Exogenous variables /////

// R: Nominal interest rate 
// z: Productivity shock
varexo R, z;

////////// The model //////////

model;
// Household budget constraint (modified equation 15)
C + (M / P) = (W / P) * N + (M(-1) / P);

// consumption Euler equation (equation 21)
C^(-sigma) = beta * (1 + R) * (P / P(+1)) * C(+1)^(-sigma);

// intra-temporal FOC describing consumption-real money balances tradeoff (equation 22)
gamma * (M / P)^(-b) = C^(-sigma) * (R / (1 + R));

// intra-temporal FOC describing consumption-labor supply tradeoff (equation 23) 
chi * N^eta = C^(-sigma) * (W / P);

// aggregate production (equation 24)
Y = z * N;

// real marginal costs of production (equation 26)
MC = nu * ((W / P) / z)

// evolution of the price level (equation 38)
P^(1 - theta) = (1 - omega) * pStar^(1 - theta) + omega * P(-1)^(1 - theta);

// inflation (equation 39)
(pStar / P) = (theta / (theta - 1)) * (K / F);

// Numerator of RHS of equation 39 from lecture notes
K = C^(1 - sigma) * MC + omega * beta * E(pi(+1)^theta * K(+1));

// Denominator of the RHS of equation 39 from lecture notes
F = C^(1 - sigma) + omega * beta * E(pi(+1)^(theta - 1) * F(+1));

// stochastic process driving productivity shocks
ln(z) = rho * ln(z(-1)) + eps_z;

// Government policy rules
R = delta * pi + v;

// Government policy shock
v = rho_v * v(-1) + eps_v
end;

