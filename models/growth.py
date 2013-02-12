import numpy as np
import mpmath as mp
import sympy as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
from scipy.optimize import fminbound
from scipy.interpolate import interp1d
from pyeconomics.utility.functions import CRRA, CARA

def checkParametersCRRA(params):
    """Checks to make sure that lifetime utility will converge for 
    given set of model parameters when household preferences are CRRA.

    """
    # extract parameters
    alpha = params['alpha']
    beta  = params['beta']
    delta = params['delta']
    theta = params['theta']
    n     = params['n']
    g     = params['g']

    if theta > 1:
        return beta * (1 + n) < 1
    else:
        return beta * (1 + g)**(1 - theta) * (1 + n) < 1

class Model(object):
    """Class for solving and simulating discrete time deterministic or
    stochastic growth models. Deterministic versions are based on the
    classic work of Solow (1956) and Ramsey (1928).

    """
    
    def __init__(self, params, k=None, c=None, **kwargs):
        """ Summary of exogenous parameters:
        
            g:     growth rate of technology
            n:     population growth rate
            delta: rate of capital depreciation
            alpha: capital share of output
            s:     savings rate
            theta: (optional) coefficient of relative (or absolute) 
                   risk aversion, depending on utility specification. 
                   The inverse of this parameter is the intertemporal
                   elasticity of substitution. 
            beta:  (optional) discount factor 
            rho:   (optional) persistence of the technology shock
            sigma: (optional) standard deviation of the technology 
                   disturbance
            
        Required arguments: 
        
            params: a dictionary of parameter values for the model, 
                    i.e. {'alpha':0.33,...}.

        Optional arguments: 
            k: A number representing the initial condition of the 
               state variable (i.e., capital per effective worker).
            c: A number representing an initial choice for the 
               control (i.e., consumption per effective worker). 

        Optional keyword arguments:
            'utility':    Household preferences. One of either CRRA or
                          CARA.
            'stochastic': If you wish to simulate a stochastic 
                          version of the Solow model, then you will
                          need set this flag to True.
            'seed':       Seed value for the random number 
                          generator (required if 'stochastic'=True).

        """
        # initial value of the state variable and the control
        self.k          = k
        self.c          = c
        
        # create the dictionary of parameter values (use copy!)
        self.param_dict = params.copy()

        # create an empty dictionary to hold SS values
        self.SS_dict    = None
        
        # household utility function 
        self.utility    = kwargs.get('utility', None)

        # validity check on utility
        if self.utility != None and self.utility not in [CRRA, CARA]:
            raise Exception, 'Invalid utility specification' 

        # if no savings rate, then default to Ramsey savings rate
        if 's' not in self.param_dict.keys():
            self.param_dict['s'] = self.get_ramseySavingsRate()
        else:
            pass
        
        ##### optional keyword arguments for stochastic model #####

        # by default you are simulating a deterministic model...
        self.stochastic = kwargs.get('stochastic', False)
        self.seed       = kwargs.get('seed', None)

        # ...but want stochastic model to nest deterministic model.
        self.z          = kwargs.get('z', 1.0)
        self.e          = kwargs.get('e', 0.0)

        # if simulating a stochastic model, the seed must be set!
        if self.stochastic == True and self.seed == None:
            raise Exception, "You must specify a seed for the RNG!"
        elif self.stochastic == True and self.seed != None:
            # sets the seed
            np.random.seed(self.seed)
            # create a disturbance term for the technology shock
            self.eps = stats.norm(0, self.param_dict['sigma'])
        else:
            pass
        
    def get_steadyStateValues(self, which='solow'): 
        """Returns a dictionary containing the steady state values of 
        capital, output, consumption, and investment per effective
        worker.
                              
        """
        # extract params
        n     = self.param_dict['n']
        g     = self.param_dict['g']
        alpha = self.param_dict['alpha']
        delta = self.param_dict['delta']

        # compute steady state of Solow model
        if which=='solow':
            s     = self.param_dict['s']

            kstar = (s / (((1 + g) * (1 + n))**alpha - \
                          ((1 + g) * (1 + n))**(alpha - 1) * \
                          (1 - delta)))**(1 / (1 - alpha))
            ystar = (1 / ((1 + g) * (1 + n)))**alpha * kstar**alpha
            cstar = (1 - s) * ystar
            istar = s * ystar
            rstar = (alpha / ((1 + g) * (1 + n))) * (ystar / kstar) - \
                    delta
            wstar = (1 - alpha) * ystar

        # compute the steady state of the Ramsey model
        elif which=='ramsey':
            beta  = self.param_dict['beta']
            theta = self.param_dict['theta']
            
            kstar = (1 + g) * (1 + n) * (alpha * beta / \
                    ((1 + g)**theta - beta * (1 - delta)))**(1 / (1 - alpha))
            ystar = (1 / ((1 + g) * (1 + n)))**alpha * kstar**alpha
            cstar = ystar + (((1 - delta) / ((1 + g) * (1 + n))) - 1) * kstar
            istar = ystar - cstar
            rstar = (alpha / ((1 + g) * (1 + n))) * (ystar / kstar) - \
                    delta
            wstar = (1 - alpha) * ystar
        else:
            raise Exception, "Which steady state values fo you want? " +\
                             "Must be one of either 'solow' or 'ramsey'!"
                             
        self.SS_dict = {'k_star':kstar, 'y_star':ystar, 'c_star':cstar, 
                        'i_star':istar, 'r_star':rstar, 'w_star':wstar}

    def get_ramseySavingsRate(self):
        """Returns: the steady state savings rate of the Ramsey 
        economy.        

        """
        # extract params
        n     = self.param_dict['n']
        g     = self.param_dict['g']
        alpha = self.param_dict['alpha']
        delta = self.param_dict['delta']
        beta  = self.param_dict['beta']
        theta = self.param_dict['theta']

        if self.utility == CRRA:
            s = alpha * beta * ((1 + g) * (1 + n) - (1 - delta)) / \
                ((1 + g)**theta - beta * (1 - delta))
        else:
            pass
        
        return s
    
    def checkStability(self, policy):
        """Computes the Jacobian matrix of partial derivatives and 
        evaluates it at the deterministic steady state, and then 
        calculates the eigenvalues and eigenvectors of the Jacobian.  

        In order for the the steady state to be dynamically stable we 
        need to have:

            1. same number of stable eigenvalues as pre-determined 
               variables (i.e., state variables)
            2. same number of unstable eigenvalue as control (i.e., 
               jump variables).

        Returns: A list containing...

            jacobian:     Array of the evaluated partial derivatives.
            eigenvalues:  The eigenvalues of the Jacobian matrix.
            eigenvectors: The eigenvectors of the Jacobian matrix.  
            
        """
        # define symbolic variables
        k = sp.var('k')
        c = sp.var('c')
        
        if policy == 'solow':
            # compute the derivative
            k_k = sp.diff(self.get_newCapital(k, self.solowPolicy(k)), k)
            # evaluate at the steady state
            k_kEval = k_k.evalf(subs={k:self.SS_dict['k_star']})

            # for 'solow' policy, Jacobian is a scalar!
            jacobian = k_kEval 

            # eigenvalues and eigenvectosr are trivial
            eigenvalues, eigenvectors = k_kEval, 1

        elif 'ramsey' in policy:
            if self.utility == CRRA:
                # compute the derivatives
                k_k = sp.diff(self.get_newCapital(k, c), k)
                k_c = sp.diff(self.get_newCapital(k, c), c)
                c_k = sp.diff(self.ramseyEulerPolicy(self.get_newCapital(k, c), c), k)
                c_c = sp.diff(self.ramseyEulerPolicy(self.get_newCapital(k, c), c), c)

            elif self.utility == CARA:
                pass
            
            # evaluate at the steady state
            evalDict = {k:self.SS_dict['k_star'], c:self.SS_dict['c_star']}
            k_kEval = k_k.evalf(subs=evalDict)
            k_cEval = k_c.evalf(subs=evalDict)
            c_kEval = c_k.evalf(subs=evalDict)
            c_cEval = c_c.evalf(subs=evalDict)

            # define the Jacobian
            jacobian = np.array([[k_kEval, k_cEval], 
                                 [c_kEval, c_cEval]], dtype='float')
        
            # calculate eigenvalues/vectors
            eigenvalues, eigenvectors = np.linalg.eig(jacobian)
            
        else:
            raise Exception, 'Invalid policy!'
        
        return [jacobian, eigenvalues, eigenvectors]

    def get_samplePath(self, policy, k0=None, c0=None, T=None, \
                       tol=None, **kwargs):
        """Generate path of length T from current value of state;  or
        generate a sample path long enough such that the difference 
        between the steady state value and current value is less than 
        tol (useful for welfare comparisons).

        Policy must be one of the following:

            'solow':           Household chooses consumption using the 
                               "rule-of-thumb", 'Solow' policy.
            'ramseyLinear':    Housdhold chooses consumption using
                               a linear approximation of the optimal 
                               decison rule.
            'ramseyEuler':     Household chooses consumption using the 
                               consumption Euler equation.
            'ramseyBellmanV':  Household chooses consumption using the 
                               policy function found using value 
                               iteration.
            'ramseyBellmanP':  Household chooses consumption using the 
                               policy function found using value 
                               iteration. 
        Arguments:
        
            k0:  Initial condition for capital per effective worker.
            c0:  Initial condition for consumption per effective 
                 worker (optional).   
            T:   Length of desired sample path.
            tol: Desired tolerance for convergence to steady state.

        Optional keyword arguments:

            'approx': (default: 'None') If you want to use simulate 
                      using the 'ramseyLinear' or 'ramseyQuadratic' 
                      consumption policies, then you may also want to 
                      approximate the equation of motion for capital 
                      as well. Approximation scheme must be one of 
                      'linear' or 'quadratic.'

        Returns: an array of shape (T,4) representing a sample path of 
        the economy. Column ordering is k, c, z, e.
        
        """
        # unpack optional keyword args
        approx = kwargs.get('approx', None)
        
        # catch possible user errors
        if self.stochastic == True and tol != None:
            raise Exception, "tol must be 'None' for stochastic model."
        if T != None and tol != None:
            raise Exception, "Either N or tol must be specified (not both!)."
        if policy not in ['solow', 'ramseyEuler', 'ramseyLinear', 
                          'ramseyBellmanV', 'ramseyBellmanP']:
            raise Exception, "Invalid policy!"
          
        # confirm that SS_dict has been populated
        if self.SS_dict == None:
            self.get_steadyStateValues()
        else:
            pass

        # initialize state and control based on user input
        if k0 != None:
            self.k = k0
        else:
            pass
            
        # generate sample path of fixed length T 
        if T != None:
            path = np.zeros((int(T), 4))
            
            for t in xrange(int(T)):
                path[t, 0] = self.k # capital
                path[t, 1] = self.c # consumption
                path[t, 2] = self.z # technology shock
                path[t, 3] = self.e # disturbance
                self.update(policy, approx)

        # generates sample path to achieve specified convergence tol
        elif tol != None:
            # compute the steady state value 
            k_star = self.SS_dict['k_star']
            
            # initialize a counter and a distance
            n_iter = 0
            dist = np.abs(self.k - k_star)
            
            while dist > tol:
                if n_iter == 0:
                    # Note that column ordering is k, c, z, e!
                    path = np.array([[self.k, self.c, self.z, self.e]])
                else:
                    step = np.array([[self.k, self.c, self.z, self.e]])
                    path = np.vstack((path, step))
                self.update(policy, approx)
                dist = np.abs(self.k - k_star)   
                n_iter += 1
        return path

    def update(self, policy, approx=None):
        """Update the state variable, k, the control, c, as well as 
        the technology shocks and their innovations (if stochastics 
        are turned on!)

        """
        # Unpack optional keyword args
        approx = kwargs.get('approx', None)
                
        # deterministic cases first!
        if self.stochastic == False:
            if policy == 'solow':
                self.c = self.solowPolicy(self.k)
                self.k = self.get_newCapital(self.k, self.c)

            elif policy == 'ramseyLinear':
                self.c = self.ramseyLinearPolicy(self.k)
                self.k = self.get_newCapital(self.k, self.c, approx)
                
            elif policy == 'ramseyEuler':
                # need to update k first!
                self.k = self.get_newCapital(self.k, self.c)
                self.c = self.ramseyEulerPolicy(self.k, self.c)
            
            elif policy == 'ramseyBellmanV':
                self.c = self.ramseyBellmanPolicy(self.k, self.c)
                self.k = self.get_newCapital(self.k, self.c)
                
            elif policy == 'ramseyBellmanP':
                self.c = self.ramseyEulerPolicy(self.k, self.c)
                self.k = self.get_newCapital(self.k, self.c)

        elif self.stochastic == True:
            if policy == 'solow':
                # draw a new value of the disturbance
                self.e = self.eps.rvs()
                # update the technology shock
                self.z = self.get_newTechShock(self.z, self.e)
                # update consumption per effective worker 
                self.c = self.solowPolicy(self.k, self.z)
                # finally, update capital per effective worker
                self.k = self.get_newCapital(self.k, self.c, self.z)
            else:
                raise Exception, "'solow' is only valid policy for ", + \
                                 "stochastic models!"
        else:
             pass   

    ##### methods implementing the various consumption policies #####
    def solowPolicy(self, k, z=1.0):
        """This rule implements the "Solow" policy. The "Solow" policy 
        is a simple "rule of thumb" policy where the household chooses 
        to consume a fixed fraction of its output in every period.

        Inputs:

            k: current period's stock of capital per effective worker.
            z: (optional) current value of the technology shock.

        Returns: current period's choice of consumption per effective
        worker.

        """
        # extract parameters
        s     = self.param_dict['s']
        alpha = self.param_dict['alpha']

        if self.stochastic == False:
            c = (1 - s) * self.get_output(k)

        elif self.stochastic == True:
            c = (1 - s) * self.get_output(k, z)
                
        return c

    def ramseyLinearPolicy(self, k):
        """Computes a first order (i.e., linear approximation to the 
        non-linear optimal consumption policy in the neighborhood of 
        the deterministic steady state.
    
        Inputs:

            k: current period's level of capital per effective worker
        
        Returns: current period's consumption per effective worker.
    
        """

        # check stability of the model
        eigVals, eigVecs = self.checkStability('ramseyEuler')[1:]

        # find the stable eigenvalue
        if np.abs(eigVals[0]) < 1:
            i = 0
        elif np.abs(eigVals[1]) < 1:
            i = 1
        else:
            raise Exception, 'Ahhh! Jacobian has no stable eigenvalue!'

        # steady state values
        cstar = self.SS_dict['c_star']
        kstar = self.SS_dict['k_star']
        
        return cstar + (eigVecs[1,i] / eigVecs[0,i]) * (k - kstar)

    def ramseyQuadraticPolicy(self):
        pass
    
    def ramseyEulerPolicy(self, kNext, c, z=1.0, zplus=1.0):
        """Via the consumption Euler equation, next period's consumption 
        per effective worker can be written as a function of current 
        period consumption and capital stock (both per effective worker).
    
        Inputs:

            k: future value of capital per effective worker
            c: current value of consumption per effective worker
        
        Returns: next period's consumption per effective worker, cplus
        
        """
        # extract params
        n     = self.param_dict['n']
        g     = self.param_dict['g']
        alpha = self.param_dict['alpha']
        delta = self.param_dict['delta']
        beta  = self.param_dict['beta']
        theta = self.param_dict['theta']
            
        # next period's return to capital
        try:
            rplus = alpha * (1 / ((1 + g) * (1 + n)))**(alpha - 1) * \
                    kNext**(alpha - 1)
        except ValueError:
            rplus = np.nan

        # euler equation depends on utility function!
        if self.utility == CRRA:
            cplus = (1 / (1 + g)) * (beta * (1 + rplus - delta))**(1 / theta) * c
        elif self.utility == CARA:
            pass
        else:
            pass
        
        return cplus

    def ramseyBellmanPolicy(self):
        pass

    ##### method defining the evolution of capital #####
    def get_output(self, k, z=1.0):
        """Returns output per effective worker. For deterministic 
        models, output is simply a function of capital per effective
        worker:

        f(k_{t-1}) = \left(\frac{1}{(1+g)(1+n)}\right)^{\alpha}k_{t-1}^{\alpha}

        For stochastic growth models, output is a function of the 
        technology shock as well as the current stock of capital per 
        effective worker.

        f(k_{t-1}, z_t) = \left(\frac{1}{(1+n)z_t}\right)^{\alpha}k_{t-1}^{\alpha}
        
        """
        # extract params
        n     = self.param_dict['n']
        g     = self.param_dict['g']
        alpha = self.param_dict['alpha']

        if self.stochastic == False:
            y = (1 / ((1 + g) * (1 + n)))**alpha * k**alpha

        elif self.stochastic == True:
            y = (1 / ((1 + n) * z))**alpha * k**alpha
        
        return y
    
    def get_newCapital(self, k, c, z=1.0, **kwargs):
        """Function that maps the capital stock per effective worker 
        at end of period t-1, k_{t-1}, to the capital stock per 
        effective worker at the end of period t, k_t. In the 
        deterministic case (i.e., when there are no technology shocks), 
        the equation of motion is:
        
        k_t = \\frac{(1-delta)}{(1+g)(1+n)}k_{t-1}+f(k_{t-1})-c_t

        In the stochastic case, the equation of motion also needs to 
        take into account the values of technology shocks.
        
        
        The houshold's stock of captial per effective worker in period
        t+1 depends on its choice of consumption in period t. Note 
        that the determinisitc model is a special case of the 
        stochastic model where z_t = 1.

        Inputs: 

            k:     capital per effective worker
            c:     consumption per effective worker
            z:     (default = 1.0) technology shock

        Returns: next period's value of capital per effective worker,
        kNext.

        BUGS: 

            * Linear appoximation scheme not working properly! 

        TODO:

            * implement quadratic approximation scheme.
            
        """ 
        # extract params
        n     = self.param_dict['n']
        g     = self.param_dict['g']
        s     = self.param_dict['s']
        alpha = self.param_dict['alpha']
        delta = self.param_dict['delta']

        # how are we approximating the equation of motion?
        approx = kwargs.get('approx')

        if self.stochastic == False:
            if approx == 'linear':
                # check stability of the model
                jacobian = self.checkStability('ramseyEuler')[0]

                # steady state values
                cstar = self.SS_dict['c_star']
                kstar = self.SS_dict['k_star']

                # evolution of linearized capital per effective worker
                kNext = kstar + jacobian[0,0] * (k - kstar) + \
                                jacobian[0,1] * (c - cstar) 

            elif approx == 'quadratic':
                pass

            elif approx == None:
                kNext = ((1 - delta) / ((1 + g) * (1 + n))) * k + \
                        self.get_output(k) - c

            else:
                raise Exception, "'approx' must be either 'linear'" + \
                                 " or 'quadratic'!"
            
        elif self.stochastic == True:
            kNext = ((1 - delta) / ((1 + n) * z)) * k + \
                    self.get_output(k, z) - c 

        return kNext

    ##### methods for generating technology shocks #####
    def get_newTechShock(self, z, eplus):
        """Defines a stochastic process representing technology shocks.

        Inputs:

            z: previous value of the technology shock.
            e: a draw from a N(0, sigma) RV.

        Returns: current value for the technology shock.
        
        """
        # extract params
        g   = self.param_dict['g']
        rho = self.param_dict['rho']

        # technology shock is log-normal
        zplus =  (1 + g)**(1 - rho) * z**rho * np.exp(eplus) 
        return zplus
        
    ########## methods for finding the Ramsey policy ##########
    def solve_forwardShoot(self, k0, c0=None, tol=1.0e-04, max_iter=1e6):
        """Computes the full, non-linear saddle path for the Ramsey 
        model using the 'forward shooting' algorithm (see Judd (1992) 
        p. 357 for details).

        Inputs:

            k0:       initial value for capital per effective worker.
            c0:       initial guess of the optimal choice for 
                      consumption per effective worker.
            tol:      how close do you want to be to steady state 
                      before you are "close enough." 
            max_iter: maximum number of iterations until you give up.

        Algorithm may fail to converge when k0 is very small because 
        initial consumption guess leads to k going negative!

        """
        # catch possible user errors
        if self.stochastic == True:
            raise Exception, "forward shooting is not valid for " + \
                             "stochastic model!"

        # extract params
        alpha = self.param_dict['alpha']
        delta = self.param_dict['delta']
        n     = self.param_dict['n']
        g     = self.param_dict['g']
        
        # values of c consistent with steady state k
        kLocus = lambda k: self.get_output(k) + \
                 (((1 - delta) / ((1 + g) * (1 + n))) - 1) * k 
                          
        # compute steady state values
        k_star = self.SS_dict['k_star']
        c_star = self.SS_dict['c_star']

        # set bounds on feasible consumption values
        if k0 <= k_star:
            c_l = 0
            c_h = kLocus(k0)
        else:
            c_l = c_star
            c_h = self.get_output(k0)

        # use 'Solow' policy as initial guess for c0 (if none provided)
        if c0 == None:
            c0 = self.solowPolicy(k0)

        # set current values of state and control 
        self.k, self.c = k0, c0
    
        # Initialize a counter
        count    = 0
        n_iter   = 0

        if k0 <= k_star:
            while n_iter < max_iter:
                self.update('ramseyEuler')
                dist = np.sqrt((self.k - k_star)**2 + (self.c - c_star)**2)
                count = count + 1

                if self.k > k_star:
                    if dist < tol:
                        break
                    else: # initial c_l too low!
                        c_l = c0
                        c0 = (c_h + c_l) / 2.
                        self.k, self.c = k0, c0
                        count = 0
                        
                elif self.c > c_star:
                    if dist < tol:
                        break
                    else: # initial c_h too high!
                        c_h = c0 
                        c0 = (c_h + c_l) / 2.
                        self.k, self.c = k0, c0
                        count = 0
                n_iter += 1
                
        else:
            while n_iter < max_iter:
                self.update('ramseyEuler')
                dist = np.sqrt((self.k - k_star)**2 + (self.c - c_star)**2)
                count = count + 1 

                if self.k < k_star:
                    if dist < tol:
                        break
                    else: # initial c_l too high!
                        c_h = c0 
                        c0 = (c_h + c_l) / 2
                        self.k, self.c = k0, c0
                        count = 0
                elif self.c < c_star:
                    if dist < tol:
                        break
                    else: # initial c_l too low!
                        c_l = c0
                        c0 = (c_h + c_l) / 2
                        self.k, self.c = k0, c0
                        count = 0
                n_iter += 1

        # Algorithm failed to converge! 
        if n_iter == max_iter:
            print "Convergence not achieved after", max_iter, "iterations!"        
            print "Lower threshold for c:", c_l
            print "Upper threshold for c:", c_h
        # Success!
        else:
            self.k, self.c = k0, c0
            solutionPath = self.get_samplePath('ramseyEuler', T=count)
            return [c0, solutionPath, count, dist, n_iter]

    def solve_valueIteration(self, k_upper=None, N=1e3, tol=1e-4, **kwargs):
        """Solves for the optimal 'Ramsey' consumption policy using 
        value iteration.

        """
        # extract parameters
        beta  = self.param_dict['beta']
        theta = self.param_dict['theta']
        g     = self.param_dict['g']
        n     = self.param_dict['n']
    
        # what kind of spacing to use for grid?
        spacing = kwargs.get('spacing', 'linear')
         
        # specify interpolation scheme
        kind    = kwargs.get('kind', 'linear')

        ##### Create the grid of points #####

        # upper bound on the state space
        kBar = float(k_upper)

        # grid of value of the state variable
        if spacing == 'linear':
            kGrid = np.linspace(1.0 / N, kBar, N)
        elif spacing == 'log':
            kGrid = np.logspace(1.0 / N, np.log(kBar + 1), N, \
                                base=np.exp(1)) - 1
        else:
            raise Exception, "Grid spacing must be 'linear' or 'log'."
        
        ##### Value iteration algorithm #####
        if self.utility == CRRA:
            currentValue = lambda c: self.utility(c, self.param_dict)
        else:
            pass

        # initialize the error and counter
        error  = 1
        n_iter = 0

        # parameter restriction requires that this is < 1
        discountFactor = beta * (1 + g)**(1 - theta) * (1 + n)

        while True:
            nextValue, nextPol = self.deterministicBellman(currentValue, \
                                                           kGrid, kind)
            error = np.max(np.abs(currentValue(kGrid) - nextValue(kGrid)))
            n_iter += 1
            if error < tol * (1 - discountFactor):
                finalValue, finalPol = nextValue, nextPol
                print "After", n_iter, "iterations, the final error " + \
                      "is", error
                break
            else:
                currentValue = nextValue 
                if n_iter % 50 == 0:
                    print "After", n_iter, "iterations, the error " + \
                          "is", error
        
        return [finalValue, finalPol, n_iter, error]

    def deterministicBellman(self, w, grid, kind='linear'):
        """The Bellman operator for the deterministic Ramsey model.
        Takes as input an initial guess of the functional form of the 
        true value function. Default initial guess is to use the 
        household's utility.
                
        Returns a list containing two callable functions. The first is
        the value function, the second is the policy function.
            
        """
        # extract the parameters
        beta  = self.param_dict['beta']
        theta = self.param_dict['theta']
        g     = self.param_dict['g']
        n     = self.param_dict['n']
    
        # parameter restriction requires that this is < 1
        discountFactor = beta * (1 + g)**(1 - theta) * (1 + n)

        def maximum(h, a, b):
            """Finds the maximum value of the function h on the 
            interval [a,b]. Returns a list containing the maximum 
            value of h and its maximizer.
    
            """
            pol, val = fminbound(lambda x: -h(x), a, b, \
                                 full_output=True)[0:2]
            return [-val, pol]

        def gamma(k):
            """Bounds the feasibility set for choice of consumption 
            per effective worker. In this model investment is 
            reversible and therefore the household can consume its
            undepreciated capital stock.

            """
            # extract parameters
            alpha = self.param_dict['alpha']
            delta = self.param_dict['delta']

            # output for given value of k
            y = self.get_output(k)
            
            return (0, y + ((1 - delta) / ((1 + g) * (1 + n))) * k)

        # containers for values and policy
        vals = []
        pols = []

        for k in grid:
            # current value function
            h = lambda c: self.utility(c, self.param_dict) + \
                discountFactor * w(self.get_newCapital(k, c))
            # value of the c that maximizes the current value function
            lower, upper = gamma(k)
            updatevector = maximum(h, lower, upper)
            # update the value function with maximum value of h
            vals.append(updatevector[0])
            # update the policy function with the maximizer of h
            pols.append(updatevector[1])

        # use interpolation to estimate value and policy funcs
        tmp_valueFunc  = interp1d(grid, vals, kind, bounds_error=False)
        tmp_policyFunc = interp1d(grid, pols, kind, bounds_error=False)

        return [tmp_valueFunc, tmp_policyFunc]

    def solve_policyIteration(self):
        pass

    ########## methods for computing distributions ##########
    def get_marginalDist(self, k0, c0, z0, e0, T=None, N=None):
        """Returns N draws of k_T, c_T, z_T, e_T starting from initial 
        values k0, c0, z0, and e0.
        
        """
        samples = np.zeros(shape=(n, 4))
        for i in xrange(N):
            self.k = k0
            self.c = c0
            self.z = z0
            self.e = e0
            for t in xrange(T):
                self.update(method)
            samples[i, 0] = self.k
            samples[i, 1] = self.c
            samples[i, 2] = self.z
            samples[i, 3] = self.e
            
        return samples

    ########## methods for generating plots ##########   
    def plot_solowDiagram(self, gridmax, N=500):
        """Generates the classic Solow diagram. 

        Inputs:

            1. gridmax: maximum value of k to use in creating the plot.
            2. N: (default=500) number of grid points to plot. 

        Returns a list containing the matplotlib objects for the plot.

        """
        # extract params
        n     = self.param_dict['n']
        g     = self.param_dict['g']
        s     = self.param_dict['s']
        alpha = self.param_dict['alpha']
        delta = self.param_dict['delta']
        
        # create the grid of points to plot
        grid = np.linspace(0, gridmax, N)

        # Create the plots
        ax           = plt.subplot(111)
        output       = ax.plot(grid, grid**alpha, 'r-')[0]
        actualInv    = ax.plot(grid, s * grid**alpha, 'g-')[0]
        breakEvenInv = ax.plot(grid, ((1 + g) * (1 + n) - (1 - delta)) * grid, 
                               'b-')[0]

        # axes, labels, title, legend, etc
        ax.set_xlabel('Capital per effective worker, $k$', fontsize=15)
        ax.set_xlim(0, gridmax)
        ax.set_ylabel('$y$, $i$, and $c$', rotation='horizontal', 
                      fontsize=15)
        ax.set_title('Classic Solow Diagram\n' + \
                     '$s=%.2f, n=%.4f, g=%.3f, \delta=%.3f$' \
                     %(s,n,g,delta), weight='bold', fontsize=20)
        
        return [ax, output, actualInv, breakEvenInv]

    def plot_solowPhaseSpace(self, gridmax, N=500, k0=None, T=None):
        """Phase diagram for the Solow economy.  

        Inputs:

            1. gridmax: maximum value of k to use in creating the plot.
            2. N: (default=500) number of grid points to plot. 
            3. k0: initial value of capital per effective worker, k.
            4. T: desired Length of time path for k.

        If optional inputs k0 and T are specified, then a timepath for
        capital per effective worker, k, of length N starting from k0 
        will be added to the plot.
        
        Returns a list containing the matplotlib objects for the plot.

        """
        # extract params
        n     = self.param_dict['n']
        g     = self.param_dict['g']
        s     = self.param_dict['s']
        alpha = self.param_dict['alpha']
        delta = self.param_dict['delta']
        
        # create the grid of points to plot
        grid = np.linspace(0, gridmax, N)

        # Create the plot
        ax = plt.subplot(111)

        kplus = self.get_newCapital(grid, self.solowPolicy(grid))
        kLocus = ax.plot(grid, kplus, 'g-')[0]
        ax.plot(grid, grid, '--', color='k')
        
        # Demarcate the 45 degree line!
        patch = mpl.patches.Arc((0,0), width=self.SS_dict['k_star'] / 2,
                                height=self.SS_dict['k_star'] / 2, 
                                theta1=0, theta2=45)
        ax.add_patch(patch)
        ax.text(self.SS_dict['k_star'] / 4, 
                self.SS_dict['k_star'] / 8, r'$45^{o}$', color='k')

        # axes, labels, title, legend, etc
        ax.set_xlabel(r'$k_{t}$')
        ax.set_xlim(0, 1.5 * self.SS_dict['k_star'])
        ax.set_ylabel(r'$k_{t+1}$', rotation='horizontal')
        ax.set_ylim(0, 1.5 * self.SS_dict['k_star'])
        ax.set_title('Phase Space for the Solow Model\n' + \
                     '$s=%.2f, n=%.4f, g=%.3f, \delta=%.3f$' \
                     %(s,n,g,delta), weight='bold', fontsize=20)
        ax.axis('tight')
        
        # if k0 and N are specified, then add a phase path
        if k0 != None and T != None:
            phase_path = self.get_solowPhasePath(k0, T)
            x_vec = phase_path[:, 0]
            y_vec = phase_path[:, 1]
            timePath = ax.plot(x_vec, y_vec)[0]

            return [ax, kLocus, timePath]
        else:
            return [ax, kLocus]
    
    def get_solowPhasePath(self, k0=None, T=None):
        """Generate and plot a timepath for k of length T plotted in 
        phase space starting from k = k0.

        """
        # initial conditions for the phase path
        self.k = k0
        self.c = self.solowPolicy(k0)

        # define a container for the path
        n_even = 2 * int(T / 2)
        path = np.zeros(shape=(n_even, 2))

        # compute the path
        for t in range(0, n_even, 2):
            path[t, 0] = self.k
            path[t + 1, 0] = self.k
            self.update('solow')
            
        self.k = k0
        path[0, 1] = 0
        
        for t in range(1, n_even - 1, 2):
            self.update('solow')
            path[t, 1] = self.k
            path[t + 1, 1] = self.k
        path[n_even - 1, 1] = self.k
        
        return path

    def plot_ramseyPhaseSpace(self, gridmax, N=500):
        """Phase diagram for the Ramsey economy.

        Inputs:

            1. gridmax: maximum value of k to use in creating the plot.
            2. N: (default=500) number of grid points to plot. 

        Returns a list containing the matplotlib objects for the plot.

        """
        # extract params
        n     = self.param_dict['n']
        g     = self.param_dict['g']
        alpha = self.param_dict['alpha']
        delta = self.param_dict['delta']

        # values of c consistent with steady state k
        kLocus = lambda k: self.get_output(k) + \
                 (((1 - delta) / ((1 + g) * (1 + n))) - 1) * k

        # Create a grid of points for plotting
        grid = np.linspace(0, gridmax, N)

        # Create the plots
        ax        = plt.subplot(111)
        k_locus   = ax.plot(grid, kLocus(grid), '-', color='orange', 
                            label=r'$\Delta k=0$')[0]
        
        c_locus   = ax.axvline(self.SS_dict['k_star'], color='black', 
                               label=r'$\Delta c=0$')
        
        ss_marker = ax.plot(self.SS_dict['k_star'], self.SS_dict['c_star'], 
                            marker='.', markersize=10, color='k')[0]

        # Add arrows to indicate out of steady-state dynamics
        x_len = 0.25 * self.SS_dict['k_star'] 
        y_len = 0.25 * self.SS_dict['c_star']   

        ax.arrow(x=0.5 * self.SS_dict['k_star'], 
                 y=0.5 * self.SS_dict['c_star'], dx=0, dy=y_len)
        ax.arrow(x=0.5 * self.SS_dict['k_star'], 
                 y=0.5 * self.SS_dict['c_star'], dx=x_len, dy=0)

        ax.arrow(x=0.5 * self.SS_dict['k_star'] + x_len, 
                 y=1.5 * self.SS_dict['c_star'], dx=0, dy=y_len)
        ax.arrow(x=0.5 * self.SS_dict['k_star'] + x_len, 
                 y=1.5 * self.SS_dict['c_star'], dx=-x_len, dy=0)

        ax.arrow(x=1.5 * self.SS_dict['k_star'], 
                 y=0.5 * self.SS_dict['c_star'] + y_len, dx=0, dy=-y_len)
        ax.arrow(x=1.5 * self.SS_dict['k_star'], 
                 y=0.5 * self.SS_dict['c_star'] + y_len, dx=x_len, dy=0)

        ax.arrow(x=1.5 * self.SS_dict['k_star'] + x_len, 
                 y=1.5 * self.SS_dict['c_star'] + y_len, dx=0, dy=-y_len)
        ax.arrow(x=1.5 * self.SS_dict['k_star'] + x_len, 
                 y=1.5 * self.SS_dict['c_star'] + y_len, dx=-x_len, dy=0)

        # axes, labels, title, legend, etc
        ax.set_xlabel('$k_t$', fontsize=15)
        ax.set_ylim(0, 2 * self.SS_dict['c_star'])
        ax.set_ylabel('$c_t$', rotation='horizontal', fontsize=15)
        ax.set_title('Phase Space for the Ramsey Model', fontsize=20, 
                     weight='bold')
        ax.legend(loc='best', frameon=False)
        
        return [ax, k_locus, c_locus, ss_marker]

    ########## Methods for welfare analysis ##########
    def get_lifetimeUtility(self, c, A0=1.0, L0=1.0, H=1.0):
        """Computes lifetime utility associated with a consumption 
        stream. Note that the function automatically accounts for 
        growth in technology and population. 
    
        Inputs: 
            c:  an array representing a consumption stream in per 
                effective worker units.
            A0: (default = 1) some value for the initial level of 
                technology.
            L0: (default = 1) some value for the initial size of 
                labor force.
            H:  (default = 1) size of the representative household.
        
        Returns: A list containing...
        
            1. Present discount value, in terms of utility, of the 
               consumption stream.
            2. The path of discounted flow utility (for plotting)
    
        """
        # extract the params
        g    = self.param_dict['g']
        n    = self.param_dict['n']
        beta = self.param_dict['beta']

        # time path
        t = np.arange(0, np.size(c), 1)
        
        # need to inflate consumption per effective worker by technology
        tech_path = A0 * (1 + g)**t
    
        # need to discount future utility from consumption 
        discount_factor = (L0 / H) * (beta * (1 + n))**t
    
        # compute the path of utility
        utility_path = discount_factor * self.utility(c * tech_path, self.param_dict)

        return [np.sum(utility_path), utility_path]

