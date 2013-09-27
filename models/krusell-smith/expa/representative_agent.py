import numpy as np

import matplotlib.pyplot as plt
from scipy import interpolate, optimize 

##### define parameters #####

# discount factor
beta = 0.99 

# coefficient of relative risk aversion
gamma = 1.0

# depreciation rate
delta = 0.025

# capital's share
alpha = 0.36

# fraction of wage paid as unemployment insurance
mu = 0.15

# time endowment
l_bar = 1 / 0.9

# aggregate shock deviation
delta_a = 0.01

##### Define some useful functions #####
def u(k, a, kplus):
    """Household has Constant Relative Risk Aversion (CRRA) preferences."""
    c = (1 - delta) * k + a * k**alpha - kplus
    if gamma != 1.0:
        return (c**(1 - gamma) - 1) / (1 - gamma)
    else:
        return np.log(c)
    
def Gamma(k, a):
    """Feasible bounds on the choice variable (i.e., consumption)."""
    return (0, a * k**alpha + (1 - delta) * k)
    
def k_star():
    """Steady state value of capital."""
    return (alpha * beta / (1 - beta * (1 - delta)))**(1 / (1 - alpha))

########## Bellman and Greedy operators ##########

def maximize(v, bounds):
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
    # extract the bounds
    lower, upper = bounds
    
    # find the minimal value
    pol = optimize.fminbound(lambda x: -v(x), lower, upper)
    
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
        tmp_v = lambda kplus: (1 - discount) * u(k, 1.0, kplus) + discount * v(kplus)
        # compute the maximizer of tmp_v
        pol = maximize(tmp_v, Gamma(k, 1.0))
        # store the new value and policy
        vals[i] = tmp_v(pol)
    
    # interpolate!
    Tv = interpolate.UnivariateSpline(grid, vals, k=order, s=0)
    
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
        tmp_v = lambda kplus: (1 - discount) * u(k, 1.0, kplus) + discount * v(kplus)
        # compute the maximizer and tmp_v(maximizer)
        pol = maximize(tmp_v, Gamma(k, 1.0))
        # store the new value and policy
        pols[i] = pol
    
    # interpolate!
    greedy_policy = interpolate.UnivariateSpline(grid, pols, k=order, s=0)
    
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
        elif n_iter % 10 == 0 and mesg == True:
            print "After", n_iter, "iterations, the change is", change
        
        current_v = next_v
    
    return [final_v, final_change]

########## Results ##########

# create a simple grid of values for capital
kMin = (1 - 0.90) * k_star()
kMax = (1 + 0.90) * k_star()
Nk   = 1000
Gk   = np.linspace(kMin, kMax, Nk)

# grid for plotting has more points!
plot_grid = np.linspace(kMin, kMax, Nk * 10)

# value of a policy zero net investment policy
init_vals = u(Gk, 1.0, Gk)
initial_v = interpolate.UnivariateSpline(Gk, init_vals, k=1, s=0)

##### Plot of the initial guess for the value function #####
plt.plot(plot_grid, initial_v(plot_grid), label='Linear')
plt.xlabel('Capital, $k$', fontsize=15)
plt.ylabel('$V_0(k)$', rotation='horizontal', fontsize=15)
plt.title('Initial value functions', fontsize=15)
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
init_vals = u(Gk, 1.0, Gk)
initial_v = interpolate.UnivariateSpline(Gk, init_vals, k=1, s=0)

# compute the optimal value and policy functions using linear B-spline
kwargs = {'u':u, 'discount': beta, 'Gamma':Gamma, 'order':1}
final_v, E_n = solve_continuousValueIteration(initial_v, 
                                              naiveUnivariateSplineBellman, 
                                              tol=0.01 * (1 - beta), pts = Gk, 
                                              mesg=True, **kwargs)

# plot the numeric and analytic solutions
plt.plot(plot_grid, final_v(plot_grid) / (1 - beta), label='Linear')
plt.xlabel('Capital, $k$', fontsize=15)
plt.ylabel('$V(k)$', rotation='horizontal', fontsize=15)
plt.title('Optimal value function', fontsize=15)
plt.legend(loc='best', frameon=False)
plt.show()

# Compute the greedy policy
kwargs = {'u':u, 'discount': beta, 'Gamma':Gamma, 'order':1}
final_kPol = naiveUnivariateSplineGreedy(final_v, **kwargs)

plt.plot(plot_grid, final_kPol(plot_grid), label='Linear')
plt.plot(plot_grid, plot_grid, 'k--')
plt.xlabel('Capital, $k$', fontsize=15)
plt.ylabel("$k'(k)$", rotation='horizontal', fontsize=15)
plt.title("Approximate optimal policy...")
plt.legend(loc='best', frameon=False)
plt.show()
