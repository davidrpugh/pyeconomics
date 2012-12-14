# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Welcome to the first lab of the inaguaral SGPE Computational Macroeconomics labs!  Your first lab will cover the Solow model from Chapter 1 of Romer's Advanced Macroeconomics.  Although Romer's treatment is excellent, I highly recommend reading Solow's original journal article entitled ['A Contribution to the Theory of Economic Growth.'](http://www.csus.edu/indiv/o/onure/econ200A/Readings/Solow.pdf)
# 
# We begin, as per usual, with some standard import statements.  All of these modules, save **sympy** and **mpmath**, you will have seen before.  We will be using [SymPy](http://sympy.org/en/index.html) to do symbolic differentiation in Python and [mpmath](http://code.google.com/p/mpmath/) to do the numerical differentiation. 

# <codecell>

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import optimize
import sympy
import mpmath

# <markdowncell>

# Task 1: Cobb - Douglas Production Function
# ------------------------------------------
# 
# How well do you understand returns to scale? Marginal productivity? Learn about the production functions using 3D plotting!  In this task, we will learn a bit of producer theory by using [matplotlib](http://matplotlib.sourceforge.net/) to plot the Cobb-Douglas production function used in the Solow, Ramsey, and Real Business Cycle (RBC) models that you will be learning in Macro I in 3D.

# <codecell>

def F(K, L, A=1):
    """
    
    Classic Cobb-Douglas production function.  Output is a function of the capital 
    stock, K, and labor, L.  Note that technology, A, is assumed to be labor-
    augmenting.
    
    Parameters, which are assumed to be already defined, are capital's share, 
    alpha, and labor's share beta. Typically, one assumes constant returns to scale 
    in which case beta = 1 - alpha.
    
    """
    return K**alpha * (A * L)**beta

# <codecell>

# Check to make sure the the function is working properly...
alpha = 1 / 3.
beta = 1 - alpha

print 'F(1, 1) =', F(K=1, L=1, A=100)
print '100^(2 / 3) =', 100**(2 / 3.)

# <markdowncell>

# Let's start by plotting the production frontier. Note that because our production function has two inputs, capital and labor, the production frontier will be a surface in 3D!

# <codecell>

# create a new Figure object 
fig = plt.figure(figsize=(8,6))

# create a 3D Axes object
ax = fig.gca(projection='3d', elev=30, azim=310)

# create a grid of (x,y) values which we will pass to function
capital = np.linspace(0, 50, 100)
labor = np.linspace(0, 50, 100)
K, L = np.meshgrid(capital, labor)

# Choose parameter values
alpha = 1 / 3.
beta = 1 - alpha
A = 1

# we will actually plot output
output = F(K, L, A)

# note the use of the new plot command!
production_frontier = ax.plot_surface(K, L, output, rstride=1, cstride=1, cmap=mpl.cm.hot, 
                                      linewidth=0, vmin=0, vmax=np.max(output), 
                                      antialiased=False)

# axes, labels, title, colorbar etc.
ax.set_xlim(0, 50)
ax.set_ylim(0, 50)
ax.set_xlabel(r'Capital ($K_{t}$)')
ax.set_ylabel(r'Labor ($L_{t}$)')
ax.set_title(r'$F(K,\ L)$ for $\alpha=%.2f, A=%.2f$' %(alpha, A), fontsize=20)
fig.colorbar(production_frontier, shrink=0.75, aspect=10)

# save the plot!
plt.savefig('Graphics/Production-frontier.png')

# display the plot!
plt.show()

# <markdowncell>

# Compare marginal product of capital with marginal product of labor...

# <codecell>

fig = plt.figure(figsize=(12, 8))

## first subplot will be be for capital ##
ax1 = fig.add_subplot(121)

# in this plot we fix L at some value, L_bar
L_bar = 10
ax1.plot(capital, F(capital, L_bar), 'r-', label=r'$F(K, \bar{L}, A)$')

# add labels, legend, title
ax1.legend(loc='best', frameon=False)
ax1.set_xlabel('Capital, $K_{t}$')
ax1.set_ylabel('Output, $Y_{t}$') 
ax1.set_ylim(0, 30)

## first subplot will be be for labor ##
ax2 = fig.add_subplot(122, sharey=ax1)

# in this plot we fix K at some value, K_bar
K_bar = 10
ax2.plot(labor, F(K_bar, labor), 'r-', label=r'$F(\bar{K}, L, A)$')

# add labels, legend, title
ax2.legend(loc='best', frameon=False)
ax2.set_xlabel('Labor, $L_{t}$')

# add a title to the figure
fig.text(0.5, 0.95, 'Comparing marginal products of $K$ and $L$', ha='center', fontsize=20, weight='bold')

# save the plot!
plt.savefig('Graphics/Comparing-marginal-products.png')

# display
plt.show()

# <markdowncell>

# We can also inspect the isoquants of our production function by creating a contour plot!

# <codecell>

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)

# Choose parameter values
alpha, A = 0.33, 1

# we will actually plot output
output = F(K, L, A)

# create the contour plot
im = ax.imshow(output, interpolation='gaussian', origin='lower', cmap=mpl.cm.hot, 
               vmin=0, vmax=np.max(output), extent=(0, 50, 0, 50))

# demarcate the contours...
CS = ax.contour(K, L, output, np.linspace(0, np.max(output), 10), colors=np.repeat('k', 10), 
                 linewidths=1, linestyles='solid')
ax.clabel(CS, inline=1, fmt='%1.2f')

# axes, labels, title, colorbar etc.
ax.set_xlim(0, 50)
ax.set_ylim(0, 50)
ax.set_xlabel(r'Capital, $K_{t}$')
ax.set_ylabel(r'Labor, $L_{t}$')
ax.set_title(r'$F(K,\ L)$ for $\alpha=%.2f, A=%.2f$' %(alpha, A), fontsize=20)
fig.colorbar(production_frontier, shrink=0.75, aspect=10)

# save the figure!
plt.savefig('Graphics/Contour-plot-for-production-frontier.png')

# display the plot!
plt.show()

# <markdowncell>

# Task 2: Coding the Solow model in Python  
# ----------------------------------------
# 
# In this task you will learn how to program a discrete time version of the Solow model in Python. Start by defining the of exogenous parameters of the Solow Model:
#     
# 1. s: savings rate 
# 2. A0: initial level of technology
# 3. g: growth rate of technology 
# 4. L0: initial level of population
# 5. n: population growth rate
# 6. delta: rate of capital depreciation
# 7. alpha: capital share of output
# 
# Model parameters are set to correspond to quarterly data (i.e., g = 0.005 or roughly 0.5% per quarter translates to 4 * 0.005 = 0.02 or 2% growth per year!)

# <codecell>

A0 = 1
g = 0.005
L0 = 100
n = 0.0025
s = 0.25
delta = 0.025
alpha = 0.33

# <markdowncell>

# Next, we need to define the key equations of the Solow model.

# <codecell>

def f(k):
    """
    
    In your lectures you will have been introduced to the intensive form 
    of the production function.  The intensive form of the production 
    function expresses output per effective worker, y, as a function of 
    the capital stock per effective worker, k.
    
    Inputs:
        1) k: capital per effective worker
        
    Returns:
        1) y: output per effective worker
        
    """
    return k**alpha
    
def capital(k):
    """
    
    Function that takes current periods capital stock per effective worker, k, 
    and returns next period's capital stock per effective worker.
    
    Inputs:
        1) k: current period's capital per effective worker 
        
    Returns:
        1) k: next period's capital per effective worker
    
    """ 
    return np.exp(-(g + n)) * (s * f(k) + np.exp(-delta) * k)

# <markdowncell>

# The equation of motion for capital stock is THE key equation of the Solow model.  Try different values for k and see what next period's k will be!

# <codecell>

print capital(2)

# <codecell>

# remember, Python function will also work accept NumPy arrays as arguments!
k = np.linspace(1, 10, num=10)
capital(k)

# <markdowncell>

# Task 3: Finding the Steady State Value of k 
# -------------------------------------------
# 
# We will do this two ways.  First, because the model is so simple, we can easily write down an analytic expression for the steady-state value in terms of the structural parameters of the model.  However, other models we will work with this year are not so simple, and for those models we will need to use numerical methods to find steady-state values.  The second part of this task simply introduces you to the tools necessary to solve for the steady-state value of capital per effective worker numerically.

# <codecell>

def analytic_k_star(): 
    """
    
    The steady-state level of capital stock per effective worker, k_star, 
    in the Solow model is a function of the 5 exogenous parameters!
    
    N.B.: The function takes no arguments because parameters are defined above!
    
    """
    return (s / (np.exp(g + n) - np.exp(-delta)))**(1 / (1 - alpha))

# <codecell>

# Display k_star for our chosen parameter values
print 'k* =', analytic_k_star()

# <codecell>

# In steady-state, next period's capital stock should equal this period's capital stock.  Does it?
print 'k* =', analytic_k_star()
print 'capital(k*) =', capital(analytic_k_star())

# <markdowncell>

# Solving for the steady-state value of k numerically. We start by defining an implicit function that describes the behavior of k near steady-state.  The utility of this approach is that, mathematically, it expresses the problem of finding the steady-state value of k in terms of finding the root of this implicit function using the scipy.optimize.fsolve() function.

# <codecell>

def implicit_capital(k):
    """
    
    Implicit function for k near steady state.  In steady-state, this function should return zero!
    
    """
    return k - capital(k)

# <codecell>

# Check that our function returns zero when evaluated at the steady-state...
print implicit_capital(analytic_k_star())

# <markdowncell>

# We will us the function scipy.optimze.fsolve() to find the root of our implicit function.  Since this is our first time using fsolve(), lets examine its help menu...

# <codecell>

optimize.fsolve?

# <codecell>

# Find the root of the implicit function using the scipy function fsolve()
numeric_k_star = optimize.fsolve(func=implicit_capital, x0=(5))

# <codecell>

# Compare the analytic and numerical solutions...note that they are identical to 10 decimal places!
print 'Analytical value for k*:', analytic_k_star() 
print 'Numerical value for k*:', numeric_k_star

# <markdowncell>

# Task 4: Assessing the Stability of the Solow Model 
# --------------------------------------------------
# 
# The stability of the Solow model is most easily seen graphically (see next task).  The Solow model is globally stable because the phase diagram cuts the $45^o$ line from above!  Mathematically, the derivative of the function capital() with respect to $k$ is less than 1 when evaluated at $k^{*}$.  Unfortunately, not all models are this simple, and this task introduces you to some of the techniques that will be used later in the course to assess dynamic stability.
# 
# As mentioned above, assessing the stability of the Solow model requires find the derivative of the equation of motion for $k$ with respect to $k$, and then evaluating this derivative at $k^{*}$.  To differentiate the equation of motion for $k$ we can either:
# 
# 1. Differentiate the equation of motion for k by hand, define a Python function for this derivative, and then evaluate this derivative at $k^{*}$.
# 2. Use SymPy to symbolically differentiate the equation of motion and evaluate the resulting derivative at $k^{*}$.
# 3. Use mpmath to numerically differentiate the equation of motion for us and evaluate the resulting derivative at $k^{*}$. 
# 
# We are going to do this in all three ways and then compare the results to make sure that they are the same.  Future labs will focus entirely on method 3 as it is the method used most often in practice to assess stability in more complicated models.

# <codecell>

# Method 1: By-hand differentiation
def analytic_capital_k(k):
    """
    
    Analytic solution for the derivative of the equation of motion for k with respect to k
    
    """
    return np.exp(-(n + g + delta)) + s * alpha * np.exp(-(n + g)) * k**(alpha - 1)

# <codecell>

# Evaluate the derivative at k_star
print 'Analytic derivative of capital() w.r.t. k evaluated at k*:', analytic_capital_k(numeric_k_star)

# <codecell>

# check with formula from slides
print 'Analytic derivative evaluated at k* (formula from slides):', alpha+(1-alpha)*np.exp(-(n + g + delta))

# <codecell>

%load_ext sympyprinting

# Method 2: Symbolic differentiation using SymPy. First, need to define k as a symbolic variable 
k = sympy.var('k')

# now we can differentiate! Display the derivative
sympy.diff(capital(k), k).simplify()

# <codecell>

# assign resulting symbolic expression to a variable
symbolic_capital_k = sympy.diff(capital(k), k)

# <codecell>

# Evaluate the symbolic derivative at k_star and compare to the analytic solution.
print 'Symbolic derivative of capital() w.r.t. k evaluated at k*:', symbolic_capital_k.evalf(n=12, subs={k:numeric_k_star})
print 'Analytic derivative of capital() w.r.t. k evaluated at k*:', analytic_capital_k(numeric_k_star)

# <codecell>

# Method 3: Numerical differentiation using Python (numerical differentiation and evaluation all in one step!)
numeric_capital_k = mpmath.diff(f=capital, x=(numeric_k_star), n=(1))

# print results
print 'Numeric derivative of capital() w.r.t. k evaluated at k*: ', numeric_capital_k
print 'Symbolic derivative of capital() w.r.t. k evaluated at k*:', symbolic_capital_k.evalf(n=12, subs={k:numeric_k_star})
print 'Analytic derivative of capital() w.r.t. k evaluated at k*:', analytic_capital_k(numeric_k_star)

# <markdowncell>

# Task 5: Graphical analysis of the Solow model using Matplotlib
# --------------------------------------------------------------
# 
# In this task you will learn how to recreate some of the basic diagrams used to analyze the Solow model using the Python library matplotlib.  First, we will create the standard Solow diagram; second, we will create a phase plot for the Solow model.  Both of these diagrams should be very familiar to you from both the textbook and your lectures.

# <codecell>

# Code for generating the standard Solow diagram

"""
# Can un-comment (i.e., remove the triple quotations above and below) to change the parameters to alter the graph
s = 0.5
g = 0
n = 0
delta = 1
alpha = 0.33
analytic_k_star()
"""

"""
# If you change the parameters, remember to revert to original parameters!
s = 0.25
g = 0.005  
n = 0.0025 
delta = 0.025 
alpha = 0.33
analytic_k_star()
"""

# Create a grid of points for plotting. Grid will be a vector of gridsize points from 0 to gridmax
gridmax, gridsize = 1.5 * analytic_k_star(), 500
grid = np.linspace(0, gridmax, gridsize)

# Create new Figure and Axes objects
fig = plt.figure()
ax = fig.add_subplot(111)

# plot output per effective worker
ax.plot(grid, f(grid), '-', color='r', label=r'$y$')
# plot savings/investment 
ax.plot(grid, s * f(grid), '-', color='g', label=r'$i_{act}$')
# plot breakeven investment
ax.plot(grid, (np.exp(g + n) - np.exp(-delta)) * grid, 'b-', label=r'$i_{br}$')

# demarcate steady-state values of y,i,c,k
ax.vlines(analytic_k_star(), 0, f(analytic_k_star()), color='k', linestyle='dashed')
ax.hlines([s * f(analytic_k_star()), f(analytic_k_star())], 0, analytic_k_star(), color='k', linestyle='dashed')

# set the axes
ax.set_xlim(0, 1.5 * analytic_k_star())
ax.set_xticks([analytic_k_star()])
ax.set_xticklabels([r'$k^*$'])

ax.set_yticks([s * f(analytic_k_star()), f(analytic_k_star())])
ax.set_yticklabels([r'$i^*$', r'$y^*$'])

# Don't forget to label your axes!
ax_xlab = ax.set_xlabel('Capital')
ax_xlab.set_family('serif')

ax_ylab = ax.set_ylabel('Output, investment, consumption')
ax_ylab.set_family('serif')

# Add a title to the plot
ax_title = ax.set_title('Deterministic Solow Model')
ax_title.set_family('serif')

# Add a legend
ax.legend(loc='right', frameon=False)

plt.savefig('Graphics/Solow-Model-Diagram.png')
plt.show()

# <codecell>

"""

Phase diagram for the deterministic Solow Model

"""

"""
# Can un-comment (i.e., remove the triple quotations above and below) to alter the graph
s = 0.5
g = 0
n = 0
delta = 1
alpha = 0.33
k_star()
"""

"""
# If you change the parameters, remember to revert to original parameters!
s = 0.25
g = 0.005  
n = 0.0025 
delta = 0.025 
alpha = 0.33
k_star()
"""

# Create a grid of points for plotting
gridmax, gridsize = 1.5 * analytic_k_star(), 500
grid = np.linspace(0, gridmax, gridsize)

# Create a new figure
fig = plt.figure()
ax = fig.add_subplot(111)

# Create the phase plot
ax.plot(grid, capital(grid), 'g-', label=r'$capital(k_{t})$')
ax.plot(grid, grid, '--', color='k')
ax.plot(analytic_k_star(), analytic_k_star(), marker='.', markersize=10, color='k')
ax.set_xlim(0, 1.5 * analytic_k_star())
ax.set_ylim(0, 1.5 * analytic_k_star())

# Demarcate the 45 degree line! 
ax.add_patch(mpl.patches.Arc(xy=(0,0), width=10, height=10, theta1=0, theta2=45))
ax.text(5, 2.5, r'$45^{o}$', color='k')

# Don't forget to label your axes!
ax.set_xlabel(r'$k_{t}$')
ax.set_ylabel(r'$k_{t+1}$')

# Add a title to the plot
ax_title = ax.set_title('Phase Space for the Deterministic Solow Model')
ax_title.set_family('serif')

# add a legend
ax.legend(loc='best', frameon=False)

plt.savefig('Graphics/Solow-Phase-Diagram.png')
plt.show()

# <markdowncell>

# Task 6: Simulating the Solow Model
# ----------------------------------
# 
# Simulating the deterministic version of the Solow model makes use of a Python class called DS (short for dynamical system!).
# 

# <codecell>

"""

Simulating the Solow Model

"""

class solowModel(object):
    
    def __init__(self, params, k=None):
        """

        Represents the model Solow Model 

        Definition of exogenous parameters of the Ramsey Model
            1) A0: initial level of technology
            2) g: growth rate of technology
            3) L0: initial level of population
            4) n: population growth rate
            5) delta: rate of capital depreciation
            6) alpha: capital share of output
            7) s: savings rate
        Attributes: 
            1) params is a dictionary of parameter values for the model, i.e. {'alpha':0.33,...}.  
               Required params are alpha, delta, n, g, s.
            2) k is a number representing the initial condition of the state variable (i.e., capital
               per effective worker.

        """
        # initial value of the state variable
        self.k          = k
        
        # create the dictionary of parameter values
        self.param_dict = params
        
        # dictionary of steady state values        
        self.SS_dict    = {'k_star':self.set_k_star(self.param_dict)} 
        
    def set_k_star(self, params): 
        """
    
        The steady-state level of capital stock per effective worker, k_star, 
        in the Solow model is a function of the exogenous parameters!
    
        """
        # extract params
        s     = params['s']
        n     = params['n']
        g     = params['g']
        alpha = params['alpha']
        delta = params['delta']
    
        return (s / (np.exp(g + n) - np.exp(-delta)))**(1 / (1 - alpha))
    
    def capital(self, k):
        """
    
        Function that takes current periods capital stock per effective worker, 
        k, and returns next period's capital stock per effective worker.
    
        """ 
        # extract params
        n     = self.param_dict['n']
        g     = self.param_dict['g']
        s     = self.param_dict['s']
        alpha = self.param_dict['alpha']
        delta = self.param_dict['delta']
    
        return np.exp(-(g + n)) * (s * k**alpha + (np.exp(-delta)) * k)
    
    def update(self):
        """
        
        Update the state variable, k
        
        """
        self.k = self.capital(self.k)

    def sample_path(self, n=None, tol=None):
        """
        
        Generate path of length n from current value of state;  or
        generante a sample path long enough such that the difference between 
        the steady state value and current value is less than tol.
        
        """
        if n != None:
            path = np.zeros(n)
            for t in range(n):
                path[t] = self.k
                self.update()
            return path
        elif tol != None:
            # compute the steady state value 
            k_star = self.SS_dict['k_star']
            
            # initialize sample path and dist
            path = np.array([self.k])
            dist = np.abs(self.k - k_star)
            while dist > tol:
                self.update()
                path = np.append(path, self.k)
                dist = np.abs(self.k - k_star)   
            return path
        else:
            print "Either length of sample path, n, or desired tolerance must be specified!"
                
        
    def phase_path(self, k0=None, n=None):
        """
        Generate a timepath for k of length n plotted in phase space starting from k = k0.
        """
        self.k = k0
        
        n_even = 2 * int(n / 2)
        path = np.zeros(shape=(n_even, 2))
        
        for t in range(0, n_even, 2):
            path[t, 0] = self.k
            path[t + 1, 0] = self.k
            self.update()
            
        self.k = k0
        path[0, 1] = self.update()
        
        for t in range(1, n_even - 1, 2):
            path[t, 1] = self.k
            path[t + 1, 1] = self.k
            self.update()
        path[n_even - 1, 1] = self.k
        
        return path
        
    def steady(self, k0=None):
        """
        
        Finds the steady state for the Solow economy, using k0 as an initial guess.
        
        """
        
        def solowSS(X):
            out = [self.F(X[0]) - X[0]]
            return out

        return optimize.fsolve(func=solowSS, x0=k0)
       
    def stability(self, k0=None):
        """
        
        Evaluates the partial derivative of capital at steady-state.
        
        """
        # compute the Jacobian
        capital_k = mp.diff(f=self.capital, x=(self.SS_dict['k_star']), n=(1))
        
        return capital_k

# <codecell>

# create a dictionary objects storing the parameter values (Initial calibration corresponds to quarterly data).
params = {'A0':1, 'g':0.005, 'L0':100, 'n':0.0025, 's':0.25, 'delta':0.025, 'alpha':0.33}

# <markdowncell>

# Time path of k plotted in phase space diagram for the deterministic Solow Model

# <codecell>

# Create two instances of the DS class
solow = solowModel(params)

# Initial conditions for the economy's capital per worker, k 
k0 = 1.9 * solow.SS_dict['k_star']
solow.k = k0

# Recreate phase diagram plot
gridmax, gridsize = 2 * solow.SS_dict['k_star'], 500
grid = np.linspace(0, gridmax, gridsize)

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)

ax.plot(grid, capital(grid))
ax.plot(grid, grid, '--', color='k')
ax.plot(solow.SS_dict['k_star'], solow.SS_dict['k_star'], marker='.', markersize=10, color='k')
ax.set_xlim(0, 2 * solow.SS_dict['k_star'])
ax.set_ylim(0, 2 * solow.SS_dict['k_star'])

# Plot path of capital from initial condition
x_vector = solow.phase_path(k0, 350)[:, 0]
y_vector = solow.phase_path(k0, 350)[:, 1]
ax.plot(x_vector, y_vector, '0.5', label='Time path of k')

# Demarcate the 45 degree line! 
ax.add_patch(mpl.patches.Arc(xy=(0,0), width=10, height=10, theta1=0, theta2=45))
ax.text(5, 2.5, r'$45^{o}$', color='k')

# don't forget to add labels!
ax.set_xlabel(r'$k_{t}$')
ax.set_ylabel(r'$k_{t+1}$', rotation='horizontal')

# Add a title to the plot
ax.set_title('Deterministic Solow Model', weight='bold')

# add a legend
ax.legend(loc='best', frameon=False)

plt.savefig('Graphics/Time-path-in-Phase-Space.png')
plt.show()

# <markdowncell>

# In the cells below there we provide code for generating time paths of k for the deterministic Solow model.  First, we simulate 400 periods (or roughly 100 years if using the quarterly parameter values), and then we generate the plot.  Note that it takes over 50 years to get close to steady state given these parameter values.

# <codecell>

## Code for generating time paths of capital per effective worker for the Solow model. ##

# Create an instance of the DS class
solow = solowModel(params)

# Create an array of initial conditions for the economy's capital per worker, k 
initial_conditions = np.linspace(15, 25, num=5)

# Create a new Figure and Axes objects
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)

# for loop plots a sample path for each initial value of k
for k0 in initial_conditions:
    solow.k = k0
    sample_path = solow.sample_path(400)
    ax.plot(sample_path)

# Add the steady state value for capital per effective worker
ax.axhline(solow.SS_dict['k_star'], ls='--', color='k', label=r'$k^{*}$')

# Don't forget to label your axes!
ax.set_xlabel('Time')
ax.set_ylabel('Capital per effective worker, k')

# Add a title to the plot
ax.set_title('Sample paths of capital per effective worker, k')

# Add the legend
ax.legend(loc='best', frameon=False)

plt.savefig('Graphics/Sample-paths-of-capital-in-Solow-Model.png')
plt.show()

# <markdowncell>

# Task 7: (Hint: Useful for the Homework!) Comparative statics  
# ------------------------------------------------------------
# 
# Start the economy in steady-state with current parameter values, and then at t=40, shock the economy by doubling the savings rate. By raising savings, consumption per effective worker must fall initially before eventually rising towards its new steady-state level.  With current parameters, doubling the savings rate results in lower levels of consumption per effective worker for roughly 4 quarters (or one year), followed by permanently higher levels of consumption per effective worker.

# <codecell>

# Recall original parameters (just to make sure!) 
params = {'A0':1, 'g':0.005, 'L0':100, 'n':0.0025, 's':0.25, 'delta':0.025, 'alpha':0.33}

# Create an instance of the DS class
solow = solowModel(params)
solow.k = solow.SS_dict['k_star']

# Create an array of zeros in which to store the sample path...
k_samplePath = np.zeros(400)
y_samplePath = np.zeros(400)
c_samplePath = np.zeros(400)

for i in range(40):
    k_samplePath[i] = solow.k
    y_samplePath[i] = f(solow.k)
    c_samplePath[i] = (1 - params['s']) * f(solow.k)
    solow.update()

# Set a new savings rate and set the new k_star()
params['s'] = 0.5
solow.SS_dict['k_star'] = solow.set_k_star(params)

# Initialize a new     
for i in range(40, 400):
    k_samplePath[i] = solow.k
    y_samplePath[i] = f(solow.k)
    c_samplePath[i] = (1 - params['s']) * f(solow.k)
    solow.update()

# Create a new Figure object
fig = plt.figure(figsize=(12,8))

# create the first subplot
ax1 = fig.add_subplot(311)
ax1.plot(k_samplePath, color='g')    
ax1.set_ylim(20, 65)

# Add the new steady state value for capital per effective worker
ax1.axhline(solow.SS_dict['k_star'], ls='--', color='k', label=r'$k_{1}^{*}$')

# Don't forget to label your axes!
ax1.set_ylabel('Capital, k')

# Add the legend
ax1.legend(loc='best', frameon=False)

# Create a second subplot
ax2 = fig.add_subplot(312, sharex=ax1)
ax2.plot(y_samplePath, color='r')

# Add the new steady state value for output per effective worker
ax2.axhline(f(solow.SS_dict['k_star']), ls='--', color='k', label=r'$y_{1}^{*}$')

# Don't forget to label your axes!
ax2.set_ylabel('Output, y')

# Add the legend
ax2.legend(loc='best', frameon=False)

# Create a third subplot
ax3 = fig.add_subplot(313, sharex=ax1)
ax3.plot(c_samplePath, color='b')
ax3.set_xlabel('Time')

# Add the new steady state value for consumption per effective worker
ax3.axhline((1 - params['s']) * f(solow.SS_dict['k_star']), ls='--', color='k', label=r'$c_{1}^{*}$')

# Don't forget to label your axes!
ax3.set_ylabel('Consumption, c')

# Add the legend
ax3.legend(loc='best', frameon=False)

# Add a title to the plot
plt.figtext(0.5, 0.98,  'Impact of a doubling of the savings rate, s',
            ha='center', color='black', weight='bold', fontsize=20)

# improve the layout
plt.tight_layout()

plt.savefig('Graphics/Shock-to-savings.png')
plt.show()

# <markdowncell>

# Finally, let's examine the impact of a doubling of the savings rate on the growth rate and the level of consumption per worker.  Recall consumption per worker grows at rate g in steady-state.

# <codecell>

# Create a new Figure object
fig = plt.figure(figsize=(12,8))

# Create the first subplot
ax1 = fig.add_subplot(211)

# Plot the grow rate in consumption per worker
ax1.plot(g + np.diff(np.log(c_samplePath)))

# Add the new steady state value for capital per effective worker
ax1.axhline(g, ls='--', color='k', label=r'$c_{bgp}$')

# Label axes
ax1.set_ylabel('Growth rate')

# Add a legend
ax1.legend(loc='best', frameon=False)

# Add a subplot
ax2 = fig.add_subplot(212, sharex=ax1)

# Plot the level of consumption per worker
grid = np.linspace(0, np.size(c_samplePath), np.size(c_samplePath))
ax2.plot(grid, np.log(c_samplePath * A0 * np.exp(g * grid)))

# Add the old balanced growth path for consumption per worker
ax2.plot(np.log(c_samplePath[0] * A0) + g * grid, ls='--', color='k', label=r'$\mathrm{ln}\left(\frac{C}{L}\right)_{bgp}$')

# Label axes
ax2.set_xlabel('Time')
ax2.set_ylabel('ln(C / L)')

# Add a legend
ax2.legend(loc='best', frameon=False)

# Add a title to the plot
plt.figtext(0.5, 0.98,  'Impact of savings shock on consumption per worker',
            ha='center', color='black', weight='bold', fontsize=20)

# improve the layout
plt.tight_layout()

plt.savefig('Graphics/Consumption-per-worker-following-a-shock-to-savings.png')
plt.show()

# <codecell>


