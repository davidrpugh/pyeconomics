import numpy as np
import mpmath as mp
from scipy import optimize, stats

class Model(object):
    
    def __init__(self, params, k, **kwargs):
        """A discrete time version of the classic Solow (1956) model of
        economic growth. 

        Summary of exogenous parameters of the Solow Model:
        
            1. g:     growth rate of technology
            2. n:     population growth rate
            3. delta: rate of capital depreciation
            4. alpha: capital share of output
            5. s:     savings rate
            6. theta: (optional) coefficient of relative risk aversion
            7. beta:  (optional) discount factor 
            8. rho:   (optional) persistence of the technology shock
            9. sigma: (optional) standard deviation of the technology 
                      disturbance
            
        Required attributes: 
        
            1. params: a dictionary of parameter values for the model, 
                       i.e. {'alpha':0.33,...}.
            2. k:      a number representing the initial condition of 
                       the state variable (i.e., capital per effective 
                       worker.

        Optional kwargs: 

            1. 'stochastic': If you wish to simulate a stochastic 
                             version of the Solow model, then you will
                             need set this flag to True.
            2. 'seed':       a seed value for the random number 
                             generator (required if 'stochastic'=True).
         
        """
        # initial value of the state variable
        self.k          = k
        
        # create the dictionary of parameter values
        self.param_dict = params
        
        # dictionary of steady state values        
        self.SS_dict    = {'k_star':self.set_k_star()} 

        # optional keyword arguments (for stochastic model only!)
        self.stochastic = kwargs.get('stochastic', False)
        self.seed       = kwargs.get('seed', None)
        self.z0         = kwargs.get('z0', 1.0)
         
        # if simulating a stochastic model, the seed must be set!
        if self.stochastic == True and self.seed == None:
            assert "You forgot to set a seed for the RNG!"
        elif self.stochastic == True and self.seed != None:
            # set the seed
            np.random.seed(self.seed)
            # create a disturbance term for the technology shock
            self.eps = stats.norm(0, self.param_dict['sigma'])
        else:
            pass
        
    def get_utility(self, C):
        """Representative agent has constant relative risk aversion 
        (CRRA) preferences over consumption per worker, C.

        For \theta \neq 1, CRRA preferences are:
        
        u(C_t) = \frac{C_t^{1 - \theta} - 1}{1 - \theta}

        For \theta = 1, CRRA preferences are equivalent to:

        u(C_t) = ln C_t

        The parameter \theta has a dual role.  It is the coefficient 
        of relative risk aversion (the higher is theta, the more risk 
        averse is the representative agent); and it is also the inverse
        of the intertemporal elasticity of substitution (the higher is 
        theta, the less willing is the representative agent to shift
        consumption back and forth across time).
        
        Inputs: consumption per worker, C.

        Returns: utility from consuming C, u(C).
        
        """
        # extract the params
        theta = self.param_dict['theta']
    
        if theta != 1:
            return (C**(1 - theta) - 1) / (1 - theta)
        else:
            return np.log(C)

    def get_lifetimeUtility(self, c, A0=1, L0=1, H=1):
        """Computes lifetime utility associated with a consumption 
        stream. Note that the function automatically accounts for 
        growth in technology and population. 
    
        Inputs: 
            1. an array, c representing a consumption stream in per 
               effective worker units.
            2. A0: (default = 1) some value for the initial level of 
               technology.
            3. L0: (default = 1) some value for the initial size of 
               labor force.
            4. H: (default = 1) size of the representative household.
        
        Returns: A list containing...
        
            1. Present discount value, in terms of utility, of the 
               consumption stream.
            2. The path of discounted flow utility (for plotting)
    
        """
        # extract the params
        g   = self.param_dict['g']
        n   = self.param_dict['n']
        beta = self.param_dict['beta']

        # time path
        time = np.arange(0, np.size(c), 1)
        
        # need to inflate consumption per effective worker by technology
        tech_path = A0 * (1 + g)**time
    
        # need to discount future utility from consumption 
        discount_factor = (L0 / float(H)) * (beta * (1 + n))**time
    
        # compute the path of utility
        utility_path = discount_factor * self.get_utility(c * tech_path)

        return [np.sum(utility_path), utility_path]
    
    def set_k_star(self): 
        """The steady-state level of capital stock per effective 
        worker, k_star, in the Solow model is a function of the 
        exogenous parameters!
    
        """
        # extract params
        s     = self.param_dict['s']
        n     = self.param_dict['n']
        g     = self.param_dict['g']
        alpha = self.param_dict['alpha']
        delta = self.param_dict['delta']
    
        return (s / ((1 + g) * (1 + n) - (1 -delta)))**(1 / (1 - alpha))
    
    def get_newCapital(self, k, zplus):
        """Function that takes end of period t's capital stock per 
        effective worker, k_t, and and a contemporaneous technology 
        shock, z_{t+1}, and returns the end of period t+1's capital 
        stock per effective worker, kplus.

        k_{t+1} = \frac{1}{(1 + g)(1 + n)z_{t+1}}[(1 - \delta)k_t + sk_t^{\alpha}]
        
        Inputs: 

            1. capital per effective worker, k.
            2. technology shock, z (defaults to 1 in every period if
               self.deterministic == True).

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
        kplus = (1 / ((1 + g) * (1 + n) * z)) * ((1 - delta) * k + s * k**alpha) 

        return kplus
        
    def get_newTechShock(self, z, eplus):
        """Defines a stochastic process representing technology shocks.

        Inputs:

            1. z: previous value of the technology shock, z.
            2. e: a draw from a N(0, sigma) RV.

        Returns: current value for the technology shock.
        
        """
        rho = self.param_dict['rho']

        # technology shock is log-normal
        zplus = self.z**rho * np.exp(eplus)
        return zplus
    
    def update(self):
        """Update the state variable, k, as well as the technology 
        shocks and their innovations (if stochastics are turned on!)

        """
        # always update in this order!
        if self.stochastic == False:
            self.e = 0
            self.z = 1
        else:
            self.e = self.eps.rvs()
            self.z = self.get_newTechShock(self.z, self.e)
            self.k = self.get_newCapital(self.k, self.z)

    def get_samplePath(self, N=None, tol=None):
        """Generate path of length n from current value of state;  or
        generate a sample path long enough such that the difference 
        between the steady state value and current value is less than 
        tol (useful for welfare comparisons).

        Inputs:
         
            1. N: length of sample path.
            2. tol: desired tolerance for convergence to steady state.

        Returns: an array of shape (N,3) representing a sample path of 
                 the Solow economy.
        
        """
        if N !=None and tol != None:
            assert "Either N or tol must be specified (but not both!)."
        if N != None:
            path = np.array([]).reshape((t, 3))
            
            for t in xrange(N):
                path[t, 0] = np.append(path, self.k)
                path[t, 1] = np.append(path, self.z)
                path[t, 2] = np.append(path, self.e)
                self.update()
            
        elif tol != None:
            # compute the steady state value 
            k_star = self.SS_dict['k_star']
            
            # initialize sample path and dist
            path = np.array([[self.k, self.z, self.e]])
            dist = np.abs(self.k - k_star)
            
            while dist > tol:
                self.update()
                path[t, 0] = np.append(path, self.k)
                path[t, 1] = np.append(path, self.z)
                path[t, 2] = np.append(path, self.e)
                dist = np.abs(self.k - k_star)   

        return path

    def get_marginalDist(self, k0, z0, e0, T=None, N=None):
        """Returns n draws of k_T, z_T, e_T starting from initial 
        values k0, z0, and e0.
        
        """
        samples = np.zeros(shape=(n, 3))
        for i in range(n):
            self.k = k0
            self.z = z0
            self.e = e0
            for t in range(T):
                self.update()
            samples[i, 0] = self.k
            samples[i, 1] = self.z
            samples[i, 2] = self.e
        return samples
        
    def plot_phaseSpace(self, k0=None, N=None):
        """Generate and plot a timepath for k of length N plotted in 
        phase space starting from k = k0.

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
        Requires an initial guess for the steady state value of k, k0.
        Returns a list containing the steady state value of k.

        """
        # function to be optimized
        def solowSS(X):
            out = [self.F(X[0]) - X[0]]
            return out

        return optimize.fsolve(func=solowSS, x0=k0)
       
    def checkStability(self):
        """Computes the Jacobian matrix of partial derivatives and 
        evaluates it at steady state, and then calculates the 
        eigenvalues and eigenvectors of the Jacobian.  

        In order for the the steady state to be dynamically stable, 
        we need to have one stable eigenvalue (i.e., one eigenvalue 
        less than unity).

        Returns: A list containing the Jacobian evaluated at steady 
        state.
            
        """
        # compute the Jacobian
        capital_k = mp.diff(f=self.capital, x=(self.SS_dict['k_star']), 
                            n=(1))
        jacobian = np.array([capital_k])
        
        return [jacobian]
