import numpy as np
import mpmath as mp
from scipy import optimize

class ramseyModel(object):
    """A discrete time version of the Ramsey-Cass-Koopmans model. 

    Firms:
        
    In the Ramsey-Cass-Koopmans model there are a large number of identical 
    firms, each having access to a common constant returns to scale (CRTS) 
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

    Households:

    There are a large number of identical households. The size of each 
    household grows exogenously at rate n.

    \frac{L_{t}}{H} = \frac{L_{0}}{H}e^{nt}

    Summary of the exogenous parameters of the Ramsey Model:
        
            1. theta: coefficient of relative risk aversion 
            2. beta (rho): discount factor (rate)
            3. A0: initial level of technology
            4. g: growth rate of technology
            5. L0: initial level of population
            6. n: population growth rate
            7. H: household size
            8. delta: rate of capital depreciation
            9. alpha: capital share of output
            
    """
    def __init__(self, params, k=None, c=None):
        """ Initializes a ramseyModel object.        

        Attributes: 
        
            1. params: a dictionary of parameters and their values 
               (i.e., {'theta':2.5, 'alpha':0.33, ...}). 
            2. k: an initial condition for the state variable k, capital per 
               effective worker.
            3. c: an initial condition for the control variable c, consumption 
               per effective worker.
        
        """
        # current value of state variable, k
        self.k            = k
        # current value of the control variable, c
        self.c            = c
        # dictionary of parameter values
        self.param_dict   = params
        # dictionary of steady state values        
        self.SS_dict      = {'k_star':self.set_k_star(), 
                             'c_star':self.set_c_star(),
                             's_star':self.set_s_star()}
        
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
        """The steady-state level of capital stock per effective worker, k^*, 
        in the Ramsey model is a function of the exogenous parameters.

        k^* = \left(\frac{\alpha e^{-\rho}}{e^{\theta g} - e^{-(\rho + \delta)}}\right)^{\frac{1}{1-\alpha}}

        Returns: the steady state level of capital per effective worker.
        
        """
        # extract params
        n     = self.param_dict['n']
        g     = self.param_dict['g']
        alpha = self.param_dict['alpha']
        delta = self.param_dict['delta']
        rho   = self.param_dict['rho']
        theta = self.param_dict['theta']

        num = alpha * np.exp(-rho) 
        denom = (np.exp(theta * g) - np.exp(-(rho + delta)))
            
        return (num / denom)**(1 / (1 - alpha))
    
    def set_c_star(self): 
        """The steady-state level of consumption per effective worker, c_star, 
        in the Ramsey model is a direct function of the exogenous parameters 
        and the steady-state level of capital stock per effective worker.

        c^* = k^{*\alpha} + (e^{-\delta} - e^{n+g})k^*

        Returns: the steady state level of consumption per effective worker.
        
        """
        # extract params
        n      = self.param_dict['n']
        g      = self.param_dict['g']
        alpha  = self.param_dict['alpha']
        delta  = self.param_dict['delta']
        k_star = self.set_k_star()
        
        return k_star**alpha + (np.exp(-delta) - np.exp(n + g)) * k_star
    
    def set_s_star(self):
        """Steady state savings rate of the Ramsey economy is a direct 
        function of the exogenous parameters of the model.

        s^* = \frac{\alpha e^{-\rho}(e^{-\delta} - e^{n+g})}{e^{\theta g} - e^{-(\rho + \delta)}} 

        Returns: the steady state savings rate.
        
        """
        # extract params
        n     = self.param_dict['n']
        g     = self.param_dict['g']
        alpha = self.param_dict['alpha']
        delta = self.param_dict['delta']
        rho   = self.param_dict['rho']
        theta = self.param_dict['theta']
        
        s     = alpha * np.exp(-rho) * (np.exp(n + g) - np.exp(-delta)) / \
                (np.exp(theta * g) - np.exp(-(rho + delta)))
        return s
    
    def capital(self, k, c):
        """Next period's capital stock per effective worker can be written as a 
        function of current period consumption and capital stock (both per 
        effective worker):

        k_{t+1} = e^{-(n + g)}(k_{t}^{\alpha} + e^{-\delta}k_{t} - c_{t}
        
        Inputs:
        
            1. k: current period's level of capital per effective worker
            2. c: current period's level of consumption per effective worker
        
        Returns: next period's capital stock per effective worker, kplus.
        
        """
        # extract params
        n     = self.param_dict['n']
        g     = self.param_dict['g']
        alpha = self.param_dict['alpha']
        delta = self.param_dict['delta']

        # next period's capital per effective worker
        kplus = np.exp(-(n + g)) * (k**alpha + np.exp(-delta) * k - c)

        return kplus
    
    def euler(self, k, c):
        """
    
        Via the consumption Euler equation, next period's consumption per
        effective worker can be written as a function of current period 
        consumption and capital stock (both per effective worker).
    
        Inputs:
            1) k: next period's level of capital per effective worker
            2) c: current period's level of consumption per effective worker
        
        Returns: next period's consumption per effective worker, cplus
        
        """
        # extract params
        n     = self.param_dict['n']
        g     = self.param_dict['g']
        alpha = self.param_dict['alpha']
        delta = self.param_dict['delta']
        rho   = self.param_dict['rho']
        theta = self.param_dict['theta']

        # nest periods return to capital
        rplus = alpha * self.capital(k, c)**(alpha - 1)
        
        # next periods consumption per effective worker
        cplus = np.exp(-g) * (np.exp(-rho) * (np.exp(-delta) + rplus))**(1 / theta) * c

        return cplus
    
    def update(self):
        """Update the state variables according to: 

        kplus = capital(k, c)
        cplus = euler(k, c)
        
        It is CRUCIAL that k is updated PRIOR to updating c (otherwise the 
        timing of the model will be wrong!).

        """
        # remember to always update k first!
        self.k = self.capital(self.k, self.c) 
        self.c = self.euler(self.k, self.c)

    def sample_path(self, N=None):
        """Generates a sample path of the Ramsey economy of length N starting 
        from the current state.

        Inputs: length of sample path, N.

        Returns: an array of shape (N, 2) representing a sample path of the 
                 Ramsey economy.
                 
        """
        path = np.zeros(shape=(N, 2))
        
        for t in xrange(N):
            path[t, 0] = self.k
            path[t, 1] = self.c
            self.update()
        
        return path

    def get_numericSteadyState(self, k0=None, c0=None):
        """Finds the steady state for the Ramsey economy using fsolve.

        Inputs:

            1. k0: an initial guess for the steady state value of k.
            2. c0: an initial guess for the steady state value of c.

        Returns: a list containing the steady state values.

        """
        # function to be optimized at steady state
        def ramseySS(X):
            out = [self.capital(X[0], X[1]) - X[0]]
            out.append(self.euler(X[0], X[1]) - X[1])
            return out

        return optimize.fsolve(func=ramseySS, x0=(k0, c0))
        
    def checkStability(self):
        """Computes the Jacobian matrix of partial derivatives and evaluates
        it at steady state, and then calculates the eigenvalues and 
        eigenvectors of the Jacobian.  

        In order for the the steady state to be dynamically stable, we need 
        to have one stable eigenvalue (i.e., one eigenvalue less than unity) 
        and one unstable eigenvalue (i.e., one eigenvalue greater than unity).

        Returns: A list containing...

            1. jacobian: a (2, 2) array of the evaluated partial derivatives.
            2. eigenvalues: the eigenvalues of the Jacobian matrix.
            3. eigenvectors: the eigenvectors of the Jacobian matrix.
            
        """
        # want to evaluate partial derivatives at steady state
        SS = (self.SS_dict['k_star'], self.SS_dict['c_star'])
            
        # calculate partial derivatives
        capital_c = mp.diff(f=self.capital, x=SS, n=(0, 1))
        capital_k = mp.diff(f=self.capital, x=SS, n=(1, 0))
        euler_c   = mp.diff(f=self.euler, x=SS, n=(0, 1))
        euler_k   = mp.diff(f=self.euler, x=SS, n=(1, 0))
        
        # define the Jacobian
        jacobian = np.array([[capital_k, capital_c], 
                             [euler_k, euler_c]], dtype='float')
        
        # calculate eigenvalues/vectors
        eigenvalues, eigenvectors = np.linalg.eig(jacobian)
        
        return [jacobian, eigenvalues, eigenvectors]

    def forward_shoot(self, k0=None, c0=None, tol=1.5e-08):
        """Computes the full, non-linear saddle path for the Ramsey model 
        using the 'forward shooting' algorithm (see Judd (1992) p. 357 for 
        details).

        Inputs:

            1. k0: initial value for capital per effective worker, k.
            2. c0: initial guess of the optimal choice for consumption
               per effective worker, c.
            3. tol: how close do you want to be to steady state before 
               you are "close enough." 
               
        """
        # compute steady state values
        k_star, c_star = self.SS_dict['k_star'], self.SS_dict['c_star']
        
        if k0 <= k_star:
            c_l = 0
            c_h = c_star
        else:
            c_l = c_star
            c_h = k0**alpha
        c0 = (c_h + c_l) / 2
        self.k, self.c = k0, c0
    
        # Initialize a counter
        count  = 0
        n_iter = 0
        
        # Forward Shooting Algorithm
        while 1:
            self.update()
            dist = np.abs(((self.k - k_star)**2 + (self.c - c_star)**2)**0.5)
            count = count + 1
            if k0 <= k_star:
                if self.k > k_star:
                    if dist < tol:
                        break
                    else: # initial c_l too low!
                        c_l = c0
                        c0 = (c_h + c_l) / 2
                        self.k, self.c = k0, c0
                        count = 0
                if self.c > c_star:
                    if dist < tol:
                        break
                    else: # initial c_h too high!
                        c_h = c0 
                        c0 = (c_h + c_l) / 2
                        self.k, self.c = k0, c0
                        count = 0
            else:
                if self.k < k_star:
                    if dist < tol:
                        break
                    else: # initial c_l too high!
                        c_h = c0 
                        c0 = (c_h + c_l) / 2
                        self.k, self.c = k0, c0
                        count = 0
                if self.c < c_star:
                    if dist < tol:
                        break
                    else: # initial c_l too low!
                        c_l = c0
                        c0 = (c_h + c_l) / 2
                        self.k, self.c = k0, c0
                        count = 0
                
        self.k, self.c = k0, c0
        solutionPath = self.sample_path(count)

        return [self.c, solutionPath, count, dist]
