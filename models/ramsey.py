import numpy as np
import sympy as sp
import mpmath as mp
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.integrate import odeint
from scipy.spatial.distance import euclidean

class ramseyModel(object):
    """Class for solving and simulating discrete and continuous-time 
    deterministic, optimal growth models.

    """
    def __init__(self, params, k=None, c=None, timing='discrete'):
        """ Initializes a ramseyModel object.        

        Attributes: 
        
            params: A dictionary of parameters and their values 
                    (i.e., {'theta':2.5, 'alpha':0.33, ...}). 

                     Parameters of the Ramsey Model:
        
                     alpha:      capital share of output
                     beta (rho): discount factor (rate)
                     theta:      coefficient of relative risk aversion 
                     g:          growth rate of technology
                     n:          population growth rate
                     delta:      rate of capital depreciation

            k:       An initial condition for the state variable k, 
                     capital per effective worker.
            c:       An initial condition for the control variable c, 
                     consumption per effective worker.
            timing:  Either 'discrete' or continuous' depending. 
                     Default is 'discrete.'
        
        """
        # current value of state variable, k
        self.k          = k
        # current value of the control variable, c
        self.c          = c
        # are we working in discrete or continuous time?
        if timing not in ['discrete', 'continuous']:
            raise Exception, "Timing must be either 'discrete', " + \
                             "'continuous'!"
        else: 
            self.timing     = timing

        # create the dictionary of parameter values (use copy!)
        self.param_dict = params.copy()
        
        # create the dictionary of steady state values        
        self.SS_dict    = None
        self.get_steadyStateValues()

        ##### Compute jacobian, and eigen* #####
        jac, eigVals, eigVecs, ind = self.checkStability()
        
        self.jacobian = jac    
        self.eigVals  = eigVals
        self.eigVecs  = eigVecs
        self.index    = ind

    ########## Steady state values and stability analysis ##########
    def get_steadyStateValues(self): 
        """Returns a dictionary containing the steady state values of 
        capital, output, consumption, and investment per effective
        worker, the net interest rate, and the real wage per effective
        worker.
                              
        """
        # extract params
        alpha = self.param_dict['alpha']
        delta = self.param_dict['delta']
        theta = self.param_dict['theta']
        n     = self.param_dict['n']
        g     = self.param_dict['g']

        # compute steady state values
        if self.timing == 'discrete':
            if 'beta' not in self.param_dict.keys():
                self.param_dict['beta'] = np.exp(-self.param_dict['rho'])
            beta  = self.param_dict['beta']
            
            kbar = (alpha * beta / ((1 + g)**theta - beta * (1 - delta)))**(1 / (1 - alpha))
            cbar = self.get_output(kbar) - ((1 + g) * (1 + n) - (1 - delta)) * kbar
            ybar = self.get_output(kbar)
            ibar = ybar - cbar
            wbar = self.get_realWage(kbar)
            rbar = self.get_interestRate(kbar)
        
        elif self.timing == 'continuous':
            if 'rho' not in self.param_dict.keys():
                self.param_dict['rho'] = -np.log(self.param_dict['beta'])
            rho   = self.param_dict['rho']
            
            kbar = (alpha / (delta + rho + theta * g))**(1 / (1 - alpha))
            cbar = kbar**alpha - (n + g + delta) * kbar
            ybar = self.get_outputContinuous([kbar, cbar])
            ibar = ybar - cbar
            wbar = self.get_realWageContinuous([kbar, cbar])
            rbar = self.get_interestRateContinuous([kbar, cbar]) 
            
        # compute fraction of output saved in steady state
        self.param_dict['s'] = (ybar - cbar) / ybar

        self.SS_dict = {'k_bar':kbar, 'y_bar':ybar, 'c_bar':cbar, 
                        'i_bar':ibar, 'r_bar':rbar, 'w_bar':wbar}
        
    def get_numericSteadyState(self, k0=None, c0=None):
        """Finds the steady state for the Ramsey economy using fsolve.
        Method is the same for both the continuous and discrete time
        versions of the model.

        Arguments:

            k0: An initial guess for the steady state value of k.
            c0: An initial guess for the steady state value of c.

        Returns: a list containing the steady state values.
            
        """
        # function to be optimized at steady state
        def ramseySS(X, timing=self.timing):
            if self.timing == 'discrete':
                out = [self.get_nextCapital(X[0], X[1]) - X[0]]
                out.append(self.get_nextConsumption(X[0], X[1]) - X[1])

            elif self.timing == 'continuous':
                out = [self.capitalContinuous(X), \
                       self.consumptionContinuous(X)]
                
            return out

        return optimize.fsolve(func=ramseySS, x0=(k0, c0), args=(self.timing))

    def checkStability(self):
        """Computes the Jacobian matrix of partial derivatives and 
        evaluates it at the deterministic steady state, and then 
        calculates the eigenvalues and eigenvectors of the Jacobian.  

        In order for the the steady state to be dynamically stable we 
        need to have:

            1. Same number of stable eigenvalues as pre-determined 
               variables (i.e., state variables).
            2. Same number of unstable eigenvalue as control (i.e., 
               jump variables).

        Returns: A list containing...

            jacobian:     Array of the evaluated partial derivatives.
            eigenvalues:  The eigenvalues of the Jacobian matrix.
            eigenvectors: The eigenvectors of the Jacobian matrix.  
            
        """
        if self.timing == 'discrete':
            # define symbolic variables
            k = sp.var('k')
            c = sp.var('c')
        
            # consumption depends on next period's capital!
            kplus = self.get_nextCapital(k, c)
            ramseySystem = sp.Matrix([self.get_nextCapital(k, c), \
                                      self.get_nextConsumption(kplus, c)])

            # define the Jacobian
            evalDict = {k:self.SS_dict['k_bar'], c:self.SS_dict['c_bar']}
            jac = ramseySystem.jacobian([k, c]).evalf(n=12, subs=evalDict)
            jacobian = np.array(jac).astype('float')

        elif self.timing == 'continuous':
            # define the Jacobian
            SS = (self.SS_dict['k_bar'], self.SS_dict['c_bar'])
            capital_c = mp.diff(self.get_nextCapital, x=SS, n=(0, 1))
            capital_k = mp.diff(self.get_nextCapital, x=SS, n=(1, 0))
            euler_c   = mp.diff(self.get_nextConsumption, x=SS, n=(0, 1))
            euler_k   = mp.diff(self.get_nextConsumption, x=SS, n=(1, 0))

            jacobian = np.array([[capital_k, capital_c], 
                                 [euler_k, euler_c]], dtype='float')
            
        # calculate eigenvalues/vectors
        eigenvalues, eigenvectors = np.linalg.eig(jacobian)

        # which is the eigenvector for the stable eigenvalue
        if eigenvalues[0] < 1:
            index = 0
        elif eigenvalues[1] < 1:
            index = 1
        else:
            raise Exception, 'No stable eigenvalue!'
         
        return [jacobian, eigenvalues, eigenvectors, index]

    ########## Quantity variables ##########
    def get_nextCapital(self, k, c):
        """Evolution of capital stock per effective worker for discrete
        time version of the model.
        
        Arguments:
        
            k: Current capital per effective worker.
            c: Current consumption per effective worker.

        Returns: 

            kplus: Next period's value of capital per effective 
                   worker.
        
        """
        # extract params
        alpha = self.param_dict['alpha']
        delta = self.param_dict['delta']
        g     = self.param_dict['g']
        n     = self.param_dict['n']

        kplus = (1 / ((1 + g) * (1 + n))) * ((1 - delta) * k + \
                                             self.get_output(k) - c)
        return kplus
    
    def capitalContinuous(self, x, t=0):
        """Evolution of capital per effective worker in the continuous
        time version of the model.

        Arguments:

            x: Vector of capital and consumption per effective worker.
               Ordering is x = [k c]
            t: Time.

        """
        # extract params
        delta = self.param_dict['delta']
        g     = self.param_dict['g']
        n     = self.param_dict['n']
        
        return  self.get_outputContinuous(x, t) - x[1] - (n + g + delta) * x[0]
    
    def get_nextConsumption(self, kplus, c):
        """Via the consumption Euler equation, next period's 
        consumption per effective worker can be written as a function 
        of current period consumption and capital stock (both per 
        effective worker).
    
        Arguments:
            kplus: Next period's level of capital per effective worker.
            c:     Current period's level of consumption per effective 
               worker.
        
        Returns: Next period's consumption per effective worker, cplus.
        
        """
        # extract params
        alpha = self.param_dict['alpha']
        beta  = self.param_dict['beta']
        delta = self.param_dict['delta']
        theta = self.param_dict['theta']
        g     = self.param_dict['g']
        n     = self.param_dict['n']
        
        # nest periods return to capital
        rplus = alpha * kplus**(alpha - 1)
        
        # next periods consumption per effective worker
        cplus = (1 / (1 + g)) * (beta * (1 + rplus - delta))**(1 / theta) * c

        return cplus

    def consumptionContinuous(self, x, t=0):
        """Consumption Euler equation describing the evolution of 
        consumption per effective worker for the continuous time 
        version of the model.

        Arguments:

            x: Vector of capital and consumption per effective worker.
               Ordering is x = [k c]
            t: Time

        """
        # extract parameters
        alpha = self.param_dict['alpha']
        rho   = self.param_dict['rho']
        delta = self.param_dict['delta']
        theta = self.param_dict['theta']
        g     = self.param_dict['g']
        n     = self.param_dict['n']
        
        return ((alpha * x[0]**(alpha - 1) - delta - rho - theta * g) / theta) * x[1]

    def get_output(self, k):
        """Output per effective worker in the discrete time version of
        the model.

        Arguments:

            k: Current period's capital per effective worker.

        Returns: 

           y: Current period's output per effective worker.

        """
        # extract params
        alpha = self.param_dict['alpha']

        return k**alpha
    
    def get_outputContinuous(self, x, t=0):
        """Output per effective worker for the continuous time version 
        of the model.

        Arguments:

            x: Vector of capital and consumption per effective worker.
               Ordering is x = [k c]
            t: Time

        """
        # extract parameters
        alpha = self.param_dict['alpha']

        return x[0]**alpha

    def get_investment(self, k, c):
        """Investment per effective worker for the discrete time 
        version of the model.

        Arguments:

            k: Current period's capital per effective worker.
            c: Current period's consumption per effective worker.

        Returns:

            i: Current period's investment per effective worker.

        """
        return self.get_output(k) - self.c
    
    def ramseySystem(self, x, t=0):
        """Continuous-time version of the model can be reduced to a 2D
        system of first-order, non-linear differential equations.

        Method is an input to SciPy's odeint which is used to solve 
        simulate the continuous time version of the model.

        """
        return [self.capitalContinuous(x, t), \
                self.consumptionContinuous(x, t)]
  
    ########## Price variables ##########
    def get_realWage(self, k):
        """Real wage per effective unit of labor for the discrete time 
        version of the model. Factor markets are assumed competitive, 
        so labor is paid its marginal product.

        """
        # extract parameters
        alpha = self.param_dict['alpha']

        return (1 - alpha) * self.get_output(k)

    def get_realWageContinuous(self, x, t=0):
        """Real wage per effective unit of labor for the continuous 
        time version of the model.

        Arguments:

            x: Vector of capital and consumption per effective worker.
               Ordering is x = [k c]
            t: Time

        """
        # extract parameters
        alpha = self.param_dict['alpha']
        
        return (1 - alpha) * self.get_outputContinuous(x)
    
    def get_interestRate(self, k):
        """Real net interest rate for the discrete time version of the
        model. Factor markets are assumed competitive, so capital is 
        paid its marginal product.

        """
        # extract parameters
        alpha = self.param_dict['alpha']

        return alpha * self.get_output(k) / k

    def get_interestRateContinuous(self, x, t=0):
        """Real wage per effective unit of labor for the continuous 
        time version of the model.

        Arguments:

            x: Vector of capital and consumption per effective worker.
               Ordering is x = [k c]
            t: Time

        """
        # extract parameters
        alpha = self.param_dict['alpha']
        delta = self.param_dict['delta']
        
        return alpha * x[0]**(alpha - 1) - delta

    ########## Simulation methods ##########
    def update(self):
        """Iterate forward the discrete time version of the model one
        time step. Remember to update capital per effective worker 
        first as consumption per effective worker in the current period
        depends on next period's capital per effective worker.
        
        """
        # remember to always update k first!
        self.k = self.get_nextCapital(self.k, self.c) 
        self.c = self.get_nextConsumption(self.k, self.c)

    def get_samplePath(self, k0=None, c0=None, T=1000):
        """Generate path of length T from current values of k and c.

        Arguments:
        
            k0:  Initial condition for capital per effective worker.
            c0:  Initial condition for consumption per effective 
                 worker.   
            T:   Length of desired sample path.

        Returns: an array of shape (T,6) representing a sample path of 
        the economy. Column ordering is k, c, y, i, w, r.
        
        """
        # initialize state and control based on user input
        if k0 != None:
            self.k = k0
        if c0 != None:
            self.c = c0
        
        # generate sample path of fixed length T 
        path = np.zeros((int(T), 6))
            
        for t in xrange(int(T)):
            path[t, 0] = self.k # capital
            path[t, 1] = self.c # consumption
            path[t, 2] = self.get_output(self.k) # output
            path[t, 3] = self.get_investment(self.k, self.c) # investment
            path[t, 4] = self.get_realWage(self.k) # real wage
            path[t, 5] = self.get_interestRate(self.k) # interest rate

            self.update()

        return path

    def get_samplePathContinuous(self, k0, c0, T, step_size=0.01):
        """Generate path of length T from current values of k and c.

        Arguments:
        
            k0:        Initial condition for capital per effective 
                       worker.
            c0:        Initial condition for consumption per effective 
                       worker.   
            T:         Length of desired sample path.
            step_size: Length of time step for odeint to use when 
                       generating the time paths of k and c.

        Returns: an array of shape (T, 6) representing a sample path of 
        the economy. Column ordering is k, c, y, i, w, r.
        
        """
        # grid of time points at which to compute the value of k
        t = np.arange(0, T, step_size)

        # use SciPy's built in ode solver to find the time path of k
        path = odeint(self.ramseySystem, [k0, c0], t)

        # output
        output = [self.get_outputContinuous(path[i,:]) for i in \
                  xrange(path.shape[0])]
        path = np.hstack((path, np.array([output]).T))

        # investment
        investment = path[:, 2] - path[:, 1] 
        path = np.hstack((path, np.array([investment]).T))

        # real wage
        real_wage = [self.get_realWageContinuous(path[i,:]) for i in \
                     xrange(path.shape[0])] 
        path = np.hstack((path, np.array([real_wage]).T))

        # interest rate
        interest_rate = [self.get_interestRateContinuous(path[i,:]) \
                         for i in xrange(path.shape[0])] 
        path = np.hstack((path, np.array([interest_rate]).T))
        
        return path
        
    ########## Solving methods ##########
    def get_kLocus(self, k):
        """Values of consumption per effective worker consistent with
        steady state capital per effective worker for discrete time 
        version of the model.

        """
        # extract params
        delta = self.param_dict['delta']
        g     = self.param_dict['g']
        n     = self.param_dict['n']

        # values of c consistent with steady state k
        return self.get_output(k) - (((1 + g) * (1 + n)) - (1 - delta)) * k

    def get_kLocusContinuous(self, x, t=0):
        """Values of consumption per effective worker consistent with
        steady state capital per effective worker for continuous time 
        version of the model.

        """
        # extract params
        delta = self.param_dict['delta']
        g     = self.param_dict['g']
        n     = self.param_dict['n']

        # values of c consistent with steady state k
        return self.get_outputContinuous(x) - (n + g + delta) * x[0]

    def solve_localPerturbation(self):
        """1st and 2nd order approximations around steady state.

        TODO:

            Implement for both discrete and continuous time!

        """
        pass
    
    def solve_forwardShoot(self, k0=None, c0=None, tol=1.5e-08):
        """Computes the full, non-linear saddle path for the discrete 
        time version of the Ramsey model using the 'forward shooting' 
        algorithm (see Judd (1992) p. 357 for details).

        Arguments:

            k0: initial value for capital per effective worker, k.
            c0: initial guess of the optimal choice for consumption
                per effective worker, c.
            tol: how close do you want to be to steady state before 
                you are "close enough." 
               
        """
        # catch user error
        if self.timing == 'continuous':
            raise Exception, "Forward shooting is unstable method " + \
                             " for continuous time model. Use " + \
                             "'solve_reverseShoot' instead!" 
                             
        # compute steady state values
        k_bar, c_bar = self.SS_dict['k_bar'], self.SS_dict['c_bar']
        
        if k0 <= k_bar:
            c_l = 0
            c_h = self.get_kLocus(k0)
        else:
            c_l = self.get_kLocus(k0)
            c_h = (1 - self.param_dict['delta']) * k0 + self.get_output(k0)

        c0 = (c_h + c_l) / 2
        self.k, self.c = k0, c0
    
        # Initialize a counter
        count  = 0
        n_iter = 0
        
        # Forward Shooting Algorithm
        while True:
            self.update()
            dist = euclidean([self.k, self.c], [k_bar, c_bar])
            count = count + 1
            if k0 <= k_bar:
                if self.k > k_bar:
                    if dist < tol:
                        break
                    else: # initial c_l too low!
                        c_l = c0
                        c0 = (c_h + c_l) / 2
                        self.k, self.c = k0, c0
                        count = 0
                elif self.c > self.get_kLocus(self.k):
                    if dist < tol:
                        break
                    else: # initial c_h too high!
                        c_h = c0 
                        c0 = (c_h + c_l) / 2
                        self.k, self.c = k0, c0
                        count = 0
            else:
                if self.k < k_bar:
                    if dist < tol:
                        break
                    else: # initial c_l too high!
                        c_h = c0 
                        c0 = (c_h + c_l) / 2
                        self.k, self.c = k0, c0
                        count = 0
                elif self.c < self.get_kLocus(self.k):
                    if dist < tol:
                        break
                    else: # initial c_l too low!
                        c_l = c0
                        c0 = (c_h + c_l) / 2
                        self.k, self.c = k0, c0
                        count = 0
                
        self.k, self.c = k0, c0
        saddlePath = self.get_samplePath(T=count)

        return [self.c, saddlePath, count, dist]

    def solve_reverseShoot(self, k0, h=1e-5, T=1000):
        """Computes the full, non-linear saddle path for a continuous
        time version of the Ramsey model using the 'reverse shooting' 
        algorithm (see Judd (1992) section 10.7 Infinite-Horizon 
        Optimal Control and Reverse Shooting, p. 355-361 for details).

        Inputs:

            k0:        Initial value for capital per effective worker.
            h:         Step size.
            T:         Length of time path to plot

        """
        # compute steady state values
        k_bar, c_bar = self.SS_dict['k_bar'], self.SS_dict['c_bar']

        # local slope of optimal policy evaluated at steady state
        Ck_prime = self.eigVecs[1, self.index] / self.eigVecs[0, self.index]

        # initial conditions
        if k0 > k_bar:
                init = [k_bar + h, c_bar + h * Ck_prime]
        elif k0 < k_bar:
                init = [k_bar - h, c_bar - h * Ck_prime]

        # reverse time...stable manifold is now unstable!
        def inverseRamseySystem(x, t=0):
            """Reversed Ramsey system of differential eqns."""
            return [-self.capitalContinuous(x, t), \
                    -self.consumptionContinuous(x, t)]

        # grid of time points at which to compute the value of k
        t = np.arange(0, T, h)
        
        # integrate backwards to k0
        saddlePath = odeint(inverseRamseySystem, init, t)

        return saddlePath

    ########## Plotting methods ##########
    def plot_phaseSpace(self, gridmax, N=500):
        """Phase space diagram for the Ramsey economy.

        Arguments:

            gridmax: maximum value of k to use in creating the plot.
            N: (default=500) number of grid points to plot. 

        Returns a list containing the matplotlib objects for the plot.

        """
        # Create a grid of points for plotting
        grid = np.linspace(0, gridmax, N)

        # Create the plots
        ax        = plt.subplot(111)

        if self.timing == 'discrete':
            k_locus = ax.plot(grid, self.get_kLocus(grid), '-', \
                              color='orange', label=r'$\Delta k=0$')[0]
        elif self.timing == 'continuous':
            k_locus = ax.plot(grid, self.get_kLocusContinuous([grid, np.ones(N)]), \
                              '-', color='orange', label=r'$\Delta k=0$')[0]
        
        c_locus   = ax.axvline(self.SS_dict['k_bar'], color='black', 
                               label=r'$\Delta c=0$')
        
        ss_marker = ax.plot(self.SS_dict['k_bar'], self.SS_dict['c_bar'], 
                            marker='.', markersize=10, color='k')[0]

        # Add arrows to indicate out of steady-state dynamics
        x_len = 0.25 * self.SS_dict['k_bar'] 
        y_len = 0.25 * self.SS_dict['c_bar']   

        ax.arrow(x=0.5 * self.SS_dict['k_bar'], 
                 y=0.5 * self.SS_dict['c_bar'], dx=0, dy=y_len)
        ax.arrow(x=0.5 * self.SS_dict['k_bar'], 
                 y=0.5 * self.SS_dict['c_bar'], dx=x_len, dy=0)

        ax.arrow(x=0.5 * self.SS_dict['k_bar'] + x_len, 
                 y=1.5 * self.SS_dict['c_bar'], dx=0, dy=y_len)
        ax.arrow(x=0.5 * self.SS_dict['k_bar'] + x_len, 
                 y=1.5 * self.SS_dict['c_bar'], dx=-x_len, dy=0)

        ax.arrow(x=1.5 * self.SS_dict['k_bar'], 
                 y=0.5 * self.SS_dict['c_bar'] + y_len, dx=0, dy=-y_len)
        ax.arrow(x=1.5 * self.SS_dict['k_bar'], 
                 y=0.5 * self.SS_dict['c_bar'] + y_len, dx=x_len, dy=0)

        ax.arrow(x=1.5 * self.SS_dict['k_bar'] + x_len, 
                 y=1.5 * self.SS_dict['c_bar'] + y_len, dx=0, dy=-y_len)
        ax.arrow(x=1.5 * self.SS_dict['k_bar'] + x_len, 
                 y=1.5 * self.SS_dict['c_bar'] + y_len, dx=-x_len, dy=0)

        # axes, labels, title, legend, etc
        ax.set_xlabel('$k_t$', fontsize=15)
        ax.set_ylim(0, 2 * self.SS_dict['c_bar'])
        ax.set_ylabel('$c_t$', rotation='horizontal', fontsize=15)
        ax.set_title('Phase Space for the Ramsey Model', fontsize=20, 
                     weight='bold')
        ax.legend(loc='best', frameon=False)
        
        return [ax, k_locus, c_locus, ss_marker]
