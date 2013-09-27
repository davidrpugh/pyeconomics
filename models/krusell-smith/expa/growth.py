import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fminbound
from scipy.interpolate import UnivariateSpline

########## Optimal growth model ##########

def u(k, kplus, z=0.0):
    """Agent has CRRA preferences."""
    c = (1 - delta) * k + np.exp(z) * k**alpha - kplus
    
    if theta != 1:
        return (c**(1 - theta)) / (1 - theta) 
    else:
        return np.log(c)
    
def Gamma(k, z=0.0):
    """The upper bound on the correspondence of feasible controls given current 
    state.

    """ 
    return (1 - delta) * k + np.exp(z) * k**alpha

def k_star():
    """Steady state value of capital."""
    return (alpha * beta / (1 - beta * (1 - delta)))**(1 / (1 - alpha))

##### Analytic solution #####

# parameters (analytic solution requires theta = delta = 1.0) 
theta = 1.0
beta  = 0.9896 
delta = 1.0
alpha = 0.4
rho   = 0.95
sigma = 0.007

# coefficients used to define the analytic value/policy functions
A = ((alpha * beta / ((1 - beta) * (1 - alpha * beta))) * np.log(alpha * beta) + 
     (1 / (1 - beta)) * np.log(1 - alpha * beta))
B = alpha / (1 - alpha * beta)
C = 1 / ((1 - alpha * beta) * (1 - beta * rho))

def analytic_v(k, z=0.0):
    """Analytic solution for the value function with logarithmic preferences and
    full deprectiation.

    Inputs:
    
        k: Current value of capital.
        z: Current value of the productivity shock.
    
    Returns: 

        v(k, z): Value of being in state k, z.

    """
    return A + B * np.log(k) + C * z

def analytic_investmentPolicy(k, z=0.0):
    """Analytic solution for the optimal policy for choosing next period's value
    of the state variable k'.

    Inputs:
    
        k: Current value of capital.
        z: Current value of the productivity shock.
    
    Returns: 

        k'(k, z): Next period's value of capital.

    """
    return (alpha * beta) * np.exp(z) * k**alpha

def analytic_consumptionPolicy(k, z=0.0):
    """Analytic solution for the optimal consumption policy function.

    Inputs:
    
        k: Current value of capital.
        z: Current value of the productivity shock.
    
    Returns: 

        c(k, z): Current period's consumption.

    """
    return (1 - alpha * beta) * np.exp(z) * k**alpha

########## Bellman and Greedy operators ##########

def optimize(v, lower, upper):
    """Wraps fminbound to find the value the minimizes -v on the closed interval
    [lower, upper] using Brent's method.

    Arguments:

        v:     A callable function representing the current value function 
               iterate.
        lower: Lower bound on the correspondence of feasible controls.
        upper: Upper bound on the correspondence of feasbile controls.

    Returns: A list containing...

        pol:    Value of the control that maximizes the function v.
        
    """  
    # find the minimal value
    pol = fminbound(lambda x: -v(x), lower, upper)
    
    return pol

    
def naiveUnivariateSplineBellman(v, u, discount, Gamma, order):
    """Naive continuous Bellman operator for the neoclassical optimal growth 
    model with inelastic labor supply. Operator uses no prior knowledge about 
    the shape of either the value or policy functions. 

    Value functions are approximated using a B-Spline representation computed 
    using the FITPACK routine CURFIT.

    Arguments:

        v:        Instance of UnivariateSpline class.
        u:        Utility function. Takes as arguments both the current and
                  future values of the state variables, k.
        discount: Effective discount factor.
        Gamma:    Correspondence of feasible values for k' given k.
        order:    The order of the B-spline representation (must satisfy
                  0 <= k <= 5).

    Returns:

        Tv: An instance of class UnivariateSpline providing the new B-spline 
            representation of the value function.

    """
    # extract the knots
    grid = v.get_knots()
    
    # containers for new value and policy functions
    vals = np.empty(grid.shape)
    
    # loop over each state and...
    for i, k in enumerate(grid):
        # current value function
        tmp_v = lambda kplus: (1 - discount) * u(k, kplus) + discount * v(kplus)
        # compute the maximizer of tmp_v
        pol = optimize(tmp_v, 0, Gamma(k, 0.0))
        # store the new value and policy
        vals[i] = tmp_v(pol)
    
    # interpolate!
    Tv = UnivariateSpline(grid, vals, k=order, s=0)
    
    return Tv

def naiveUnivariateSplineGreedy(v, u, discount, Gamma, order):
    """Computes the optimal, greedy policy for value functions approximated using 
    B-splines. Operator uses no prior knowledge about the shape of the true 
    policy functions. 

    Policy functions are approximated using a B-Spline representation computed 
    using the FITPACK routine CURFIT.

    Arguments:

        v:        Instance of UnivariateSpline class.
        u:        Utility function. Takes as arguments both the current and
                  future values of the state variables, k.
        discount: Effective discount factor.
        Gamma:    Correspondence of feasible values for k' given k.
        order:    The order of the B-spline representation (must satisfy
                  0 <= k <= 5).

    Returns:

        greedy_policy: An instance of class UnivariateSpline providing the 
                       new B-spline representation of the policy function.

    """
    # extract the knots
    grid = v.get_knots()
    
    # containers for new value and policy functions
    pols = np.zeros(grid.shape)
    
    # loop over each state and...
    for i, k in enumerate(grid):
        # current value function
        tmp_v = lambda kplus: (1 - discount) * u(k, kplus) + discount * v(kplus)
        # compute the maximizer and tmp_v(maximizer)
        pol = optimize(tmp_v, 0, Gamma(k, 0.0))
        # store the new value and policy
        pols[i] = pol
    
    # interpolate!
    greedy_policy = UnivariateSpline(grid, pols, k=order, s=0)
    
    return greedy_policy

##### Value function iteration #####
def solve_continuousValueIteration(init_v, T, tol, pts, mesg=False, **kwargs):
    """Basic implementation of Value Iteration Algorithm.

    Arguments:

        init_v:   Initial guess of the true value function. A good initial
                  guess can save a substantial amount of computational time.
        T:        A pre-defined Bellman operator.
        tol:      Convergence criterion. Algorithm will terminate when 
                  the supremum norm distance between successive value 
                  function iterates is less than tol.
        pts:      Grid of points over which to compare value function iterates.
        mesg:     Should messages be printed detailing convergence progress?
                  Default is False.

    Returns: 

        final_out: A list containing the optimal value and policy functions.

    """

    # keep track of number of iterations
    n_iter = 0

    ##### Value iteration algorithm #####
    current_v = init_v

    while True:
        next_v = T(current_v, **kwargs)
        # supremum norm convergence criterion
        change = np.max(np.abs(next_v(pts) - current_v(pts)))
        n_iter += 1
    
        # check for convergence
        if change < tol:
            if mesg == True:
                print "After", n_iter, "iterations, the final change is", change
            final_v      = next_v
            final_change = change
            break
    
        # print progress every 50 iterations
        elif n_iter % 1 == 0 and mesg == True:
            print "After", n_iter, "iterations, the change is", change
        
        current_v = next_v
    
    return [final_v, final_change]

########## Results ##########

# create a simple grid of values for capital
kMin = (1 - 0.90) * k_star()
kMax = (1 + 0.90) * k_star()
Nk   = 10
Gk   = np.linspace(kMin, kMax, Nk)

# grid for plotting has more points!
plot_grid = np.linspace(kMin, kMax, Nk * 10)

# value of a policy zero net investment policy
init_vals = u(Gk, Gk)
initial_v2 = UnivariateSpline(Gk, init_vals, k=1, s=0)

##### Plot of the initial guess for the value function #####
plt.plot(plot_grid, initial_v2(plot_grid), label='Linear')
plt.xlabel('Capital, $k$', fontsize=15)
plt.ylabel('$V_0(k)$', rotation='horizontal', fontsize=15)
plt.title('Initial value functions', fontsize=15)
plt.legend(loc='best', frameon=False)
plt.show()

##### Plot of the Bellman operator #####

# index of current capital 
i = 0

# Bellman operator using Linear approximation
tmp_v2 = lambda kplus: (1 - beta) * u(Gk[i], kplus) + beta * initial_v2(kplus)
neg_tmp_v2 = lambda kplus: -tmp_v2(kplus)

plt.plot(plot_grid, tmp_v2(plot_grid), label='Linear')
plt.xlabel('Capital, $k$', fontsize=15)
plt.ylabel('$TV_0(k)$', rotation='horizontal', fontsize=15)
plt.title('T not preserving monotonicity with StepFun?', fontsize=15)
plt.legend(loc='best', frameon=False)
plt.show()

##### Value function converges when grid is sufficently large #####

# create a simple grid of values for capital
kMin = (1 - 0.90) * k_star()
kMax = (1 + 0.90) * k_star()
Nk   = 1000
Gk   = np.linspace(kMin, kMax, Nk)

# grid for plotting has more points!
plot_grid = np.linspace(kMin, kMax, Nk * 10)

# value of a policy zero net investment policy
init_vals = u(Gk, Gk)
initial_v2 = UnivariateSpline(Gk, init_vals, k=1, s=0)

# compute the optimal value and policy functions using linear B-spline
kwargs = {'u':u, 'discount': beta, 'Gamma':Gamma, 'order':1}
final_v2, E_n = solve_continuousValueIteration(initial_v2, naiveUnivariateSplineBellman, 
                                               tol=0.01 * (1 - beta), pts = Gk, 
                                               mesg=True, **kwargs)

# plot the numeric and analytic solutions
plt.plot(plot_grid, final_v2(plot_grid) / (1 - beta), label='Linear')
plt.plot(plot_grid, analytic_v(plot_grid), 'r-', label='Analytic')
plt.xlabel('Capital, $k$', fontsize=15)
plt.ylabel('$V(k)$', rotation='horizontal', fontsize=15)
plt.title('Optimal value function using RPD', fontsize=15)
plt.legend(loc='best', frameon=False)
plt.show()

# How far it the numeric solution from the analytic solution? Not far...
print 'Approximation error for v(k) using linear B-spline:', np.max(np.abs(final_v2(Gk) / (1 - beta) - 
                                                                           analytic_v(Gk)))

# Compute the greedy policy (here is where things seem to go pear shaped for RPD)
kwargs = {'u':u, 'discount': beta, 'Gamma':Gamma, 'order':1}
final_kPol2 = naiveUnivariateSplineGreedy(final_v2, **kwargs)

plt.plot(plot_grid, final_kPol2(plot_grid), label='Linear')
plt.plot(plot_grid, analytic_investmentPolicy(plot_grid), 'r-', label='Analyic')
plt.xlabel('Capital, $k$', fontsize=15)
plt.ylabel("$k'(k)$", rotation='horizontal', fontsize=15)
plt.title("RPD has trouble approximating the optimal policy...")
plt.legend(loc='best', frameon=False)
plt.show()

# How far it the numeric solution from the analytic solution? 
print "Approximation error for k'(k) using linear B-spline:", np.max(np.abs(final_kPol2(Gk) - 
                                                                            analytic_investmentPolicy(Gk)))


final_cPol = lambda k: (1 - delta) * k + k**alpha - final_kPol2(k)

plt.plot(plot_grid, final_cPol(plot_grid), label='Linear')
plt.plot(plot_grid, analytic_consumptionPolicy(plot_grid), 'r-', label='Analyic')
plt.xlabel('Capital, $k$', fontsize=15)
plt.ylabel("$c(k)$", rotation='horizontal', fontsize=15)
plt.title("Consumption policy function")
plt.legend(loc='best', frameon=False)
plt.show()

plt.figure()
#plt.plot(plot_grid, (1 / beta) * (final_cPol(0.1) / final_cPol(plot_grid))**(-theta) - 1 + delta, label='supply, 0.1')
#plt.plot(plot_grid, (1 / beta) * (final_cPol(0.15) / final_cPol(plot_grid))**(-theta) - 1 + delta, label='supply, 0.15')
#plt.plot(plot_grid, (1 / beta) * (final_cPol(0.12) / final_cPol(plot_grid))**(-theta) - 1 + delta, label='supply, 0.20')
plt.plot(plot_grid, (1 / beta) * (final_cPol(k_star()) / final_cPol(plot_grid))**(-theta) - 1 + delta, label='supply, $k*$')
plt.plot(plot_grid, alpha * plot_grid**(alpha - 1), label='demand')
#plt.axhline((1 / beta) - 1, 'k--')
plt.xlabel('Capital, $k$', fontsize=15)
plt.ylabel("$r(k)$", rotation='horizontal', fontsize=15)
plt.title("Market for capital")
plt.legend(loc='best', frameon=False)
plt.show()


