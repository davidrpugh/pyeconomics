////////// Diamond (1965): Basic OLG model with production //////////

///// Declare variables /////
var y r w c c1 c2 s k check1 check2;

///// Declare parameter values /////
parameters alpha beta delta rho theta n g T;

// length of time period
T = 30;

// discount rate
rho = 0.03;

// discount factor
beta = (1 - rho)^T;

// population growth rate
n = (1 + 0.02)^T - 1;

// technology growth rate
g = (1 + 0.02)^T - 1;

// capital's share of output
alpha = 0.33;

// depreciation rate
delta = 1 - (1 - 0.10)^T;

// inverse elasticity of intertemporal substitution
theta = 1.0;

///// Steady state values /////

// capital per effective worker
kss  = 1.0;

// output per effective worker
yss  = (1 / ((1 + g) * (1 + n)))^alpha * kss^alpha;

// real wage
wss  = (1 - alpha) * yss;

// net interest rate
rss  = alpha * (1 / ((1 + g) * (1 + n)))^(alpha - 1) * kss^(alpha - 1) - delta; 

// savings
sss  = wss / (1 + beta^(- 1 / theta) * (1 + rss - delta)^((theta - 1) / theta));

// consumption of the young cohort
c1ss = wss - sss;

// consumption of the old cohort
c2ss = (1 / (1 + g)) * (1 + rss) * sss;

// consumption per effective worker
css  = c1ss + (1 / (1 + n)) * c2ss;  

///// The model /////
model;

// Cobb-Douglas production technology
y = (1 / ((1 + g) * (1 + n)))^alpha * k(-1)^alpha;

// capital is paid its net marginal product
r = alpha * (1 / ((1 + n) * (1 + g)))^(alpha - 1) *  k(-1)^(alpha - 1) - delta;

// labor is paid its marginal product
w = (1 - alpha) * y;

// period 1 flow of funds constraint (young spend their wages)
c1 + s = w;

// period 2 flow of funds constraint (old spend their savings)
c2 = (1 / (1 + g)) * (1 + r) * s;

// savings function
s = w / (1 + beta^(-1 / theta) * (1 + r(+1))^((theta - 1) / theta));

// aggregate consumption is the sum of consumption of both cohorts
c = c1 + (1 / (1 + n)) * c2;

// equation of motion for capital per effective worker
k = ((1 - delta) / ((1 + n) * (1 + g))) * k(-1) + s;

// check that zero profit condition is satisfied
check1 = y - (1 / ((1 + g) * (1 + n))) * (r + delta) * k(-1) - w;

// check that consumption Euler equation holds
check2 = c1^(-theta) - beta * (1 + g)^(-theta) * (1 + r(+1)) * c2^(-theta);

end;

///// Declare initial values /////

initval;
k  = kss;
c1 = c1ss;
c2 = c2ss;
r  = rss;
w  = wss;
s  = sss;
c  = css;
y  = yss;
end; 

///// Compute steady state values /////
steady;

///// Computes the eigenvalues of the linearized model /////
check;

///// Deterministic simulation /////
simul(periods = 100, datafile=diamond_1965_initvals);