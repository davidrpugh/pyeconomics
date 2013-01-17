import numpy as np
import mpmath as mp
from scipy import optimize

class solowModel(object):
    
    def __init__(self, params, k=None):
        """A discrete time version of the Solow Model. 

        Firms:
        
        In the Solow model there are a large number of identical firms, 
        each having access to a common constant returns to scale (CRTS) 
        Cobb-Douglas production technology:

        Y_{t} = F(K_{t}, A_{t}L_{t}) 
              = K_{t}^{alpha}(A_{t}L_{t})^{1 - alpha}

        or in intensive form:

        y_{t} = \frac{Y_{t}}{A_{t}L_{t}} 
              = \left(\frac{K_{t}}{A_{t}L_{t}}\right)^{\alpha} 
              = k_{t}^{\alpha} 
              = f(k_{t})

        where 0 < \alpha < 1 and thus f'(k) > 0 and f''(k) < 0.

        Firms choose to hire workers and rent capital in order to maximize 
        profits.  All markets (i.e., both for factors of production and  
        final output goods) are competitive. Firms take technology, A, as given
        when solving their optimization problems. As in the Solow model, A is
        assumed to grow exogenously as rate g.

        A_{t} = A_{0}e^{gt}

        Summary of exogenous parameters of the Ramsey Model
        
            1. A0: initial level of technology
            2. g: growth rate of technology
            3. L0: initial level of population
            4. n: population growth rate
            5. delta: rate of capital depreciation
            6. alpha: capital share of output
            7. s: savings rate
            
        Attributes: 
        
            1. params: a dictionary of parameter values for the model, 
               i.e. {'alpha':0.33,...}.
            2. k: a number representing the initial condition of the state 
               variable (i.e., capital per effective worker.

        """
        # initial value of the state variable
        self.k          = k
        
        # create the dictionary of parameter values
        self.param_dict = params
        
        # dictionary of steady state values        
        self.SS_dict    = {'k_star':self.set_k_star()} 

    def u(self, C):
        """Representative agent has constant relative risk aversion (CRRA) 
        preferences over consumption per worker, C.

        For \theta \neq 1, CRRA preferences are:
        
        u(C_{t}) = \frac{C_{t}^{1 - \theta} - 1}{1 - \theta}

        For \theta = 1, CRRA preferences are equivalent to:

        u(C_{t} = ln C_{t}

        Inputs: consumption per worker, C.

        Returns: utility from consuming C, u(C).
        
        """
        # extract the params
        theta = self.param_dict['theta']
    
        if theta != 1:
            return (C**(1 - theta) - 1) / (1 - theta)
        else:
            return np.log(C)

    def lifetimeUtility(self, c):
        """Computes lifetime utility associated with a consumption stream. Note 
        that the function automatically accounts for growth in technology and 
        population. 
    
        Inputs: an array, c representing a consumption stream in per effective 
                worker units.
        
        Returns: A list containing...
        
            1. Present discount value, in terms of utility, of the consumption 
               stream.
            2. The path of discounted flow utility (for plotting)
    
        """
        # extract the params
        A0  = self.param_dict['A0']
        g   = self.param_dict['g']
        L0  = self.param_dict['L0']
        H   = self.param_dict['H']
        n   = self.param_dict['n']
        rho = self.param_dict['rho']
    
        # need to inflate consumption per effective worker by technology
        tech_path = A0 * np.exp(g * np.arange(0, np.size(c), 1))
    
        # need to discount future utility from consumption appropriately 
        discounted_u = (L0 / H) * np.exp((n - rho) * np.arange(0, np.size(c), 1))
    
        # compute the path of utility
        utility_path = discounted_u * self.u(c * tech_path)

        return [np.sum(utility_path), utility_path]
    
    def set_k_star(self): 
        """
    
        The steady-state level of capital stock per effective worker, k_star, 
        in the Solow model is a function of the exogenous parameters!
    
        """
        # extract params
        s     = self.param_dict['s']
        n     = self.param_dict['n']
        g     = self.param_dict['g']
        alpha = self.param_dict['alpha']
        delta = self.param_dict['delta']
    
        return (s / (np.exp(g + n) - np.exp(-delta)))**(1 / (1 - alpha))
    
    def capital(self, k):
        """Function that takes current periods capital stock per effective 
        worker, k, and returns next period's capital stock per effective 
        worker.

        k_{t+1} = e^{-(n + g)}(sk_{t}^{\alpha} + e^{-\delta}k)
        
        Inputs: capital per effective worker, k.

        Returns: next period's value of capital per effective worker,
                 kplus.
    
        """ 
        # extract params
        n     = self.param_dict['n']
        g     = self.param_dict['g']
        s     = self.param_dict['s']
        alpha = self.param_dict['alpha']
        delta = self.param_dict['delta']

        # next period's capital per effective worker
        kplus = np.exp(-(g + n)) * (s * k**alpha + (np.exp(-delta)) * k) 

        return kplus
    
    def update(self):
        """Update the state variable, k according to

        kplus = capital(k)

        """
        
        self.k = self.capital(self.k)

    def sample_path(self, N=None, tol=None):
        """Generate path of length n from current value of state;  or
        generate a sample path long enough such that the difference between 
        the steady state value and current value is less than tol.

        Inputs:
         
            1. N: length of sample path.
            2. tol: desired tolerance for conversion to steady state.

        Returns: an array of shape (N,) representing a sample path of the 
                 Solow economy.
        
        """
        if N != None:
            path = np.zeros(N)
            
            for t in xrange(N):
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
            print "Either N or tol must be specified (but not both!)."
                
    def phase_path(self, k0=None, N=None):
        """Generate a timepath for k of length n plotted in phase space 
        starting from k = k0.

        """
        self.k = k0
        
        n_even = 2 * int(N / 2)
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
        
    def get_numericSteadyState(self, k0=None):
        """Finds the steady state for the Solow economy using fsolve.

        Inputs: an initial guess for the steady state value of k, k0.

        Returns: a list containing the steady state value of k.

        """
        # function to be optimized
        def solowSS(X):
            out = [self.F(X[0]) - X[0]]
            return out

        return optimize.fsolve(func=solowSS, x0=k0)
       
    def checkStability(self):
        """Computes the Jacobian matrix of partial derivatives and evaluates
        it at steady state, and then calculates the eigenvalues and 
        eigenvectors of the Jacobian.  

        In order for the the steady state to be dynamically stable, we need 
        to have one stable eigenvalue (i.e., one eigenvalue less than unity).

        Returns: A list containing the Jacobian evaluated at steady state.
            
        """
        # compute the Jacobian
        capital_k = mp.diff(f=self.capital, x=(self.SS_dict['k_star']), n=(1))
        jacobian = np.array([capital_k])
        
        return [jacobian]
