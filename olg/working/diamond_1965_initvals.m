////////// Initial values for Diamond (1965) OLG model //////////

// Capital per effective worker
k  = 1.0;

// output per effective worker
y  = (1 / ((1 + g) * (1 + n)))^alpha * k^alpha;

// real wage
w  = (1 - alpha) * y;

// net interest rate
r  = alpha * (1 / ((1 + g) * (1 + n)))^(alpha - 1) * k^(alpha - 1) - delta; 

// savings
s  = w / (1 + beta^(- 1 / theta) * (1 + r - delta)^((theta - 1) / theta));

// consumption of the young cohort
c1 = w - s;

// consumption of the old cohort
c2 = (1 / (1 + g)) * (1 + r) * s;

// consumption per effective worker
c  = c1 + (1 / (1 + n)) * c2;  
