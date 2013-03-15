////////// Basic OLG model with production //////////

///// Declare variables /////
var c c1 c2 k w r s y i check1;

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
alpha = 1 / 3;

// depreciation rate
delta = 1 - (1 - 0.10)^T;

// inverse elasticity of intertemporal substitution
theta = 1.0;

///// Steady state values /////

// capital per effective worker
kss  = 1.0;
yss  = kss^alpha;
wss  = (1 - alpha) * yss;
rss  = alpha * kss^(alpha - 1); 
sss  = wss / (1 + beta^(- 1 / theta) * (1 + rss - delta)^((theta - 1) / theta));
c1ss = wss - sss;
c2ss = (1 + n) * (1 + g) * (1 + rss - delta) * kss;
css  = c1ss + (1 / ((1 + n) * (1 + g))) * c2ss;  
iss  = yss - css;

///// The model /////
model;

// accounting identity
y = c + i;

// aggregate consumption is the sum of consumption of both cohorts
c = c1 + (1 / ((1 + n) * (1 + g))) * c2;

// zero-profit condition
y = w + r* k(-1);

// equation of motion for capital per effective worker
k = (1 / ((1 + n) * (1 + g))) * ((1 - delta) * k(-1) + i);

// capital is paid its marginal product
r = alpha * k(-1)^(alpha - 1);

// labor is paid its marginal product
w = (1 - alpha) * y;

// period 1 flow of funds constraint (young spend their wages)
c1 = w - s;

// period 2 flow of funds constraint (old spend their savings)
c2 = (1 + n) * (1 + g) * (1 + r - delta) * k(-1);

// savings function
s = w / (1 + beta^(-1 / theta) * (1 + r(+1) - delta)^((theta - 1) / theta));

// check that investment equals savings
check1 = k - (1 / ((1 + g) * (1 + n))) * s;
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

///// Check local stability of steady state /////
check;