import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import optimize, stats

class Model(object):
    
    def __init__(self, params, k=None, **kwargs):
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
        
        # create the dictionary of parameter values (use copy!)
        self.param_dict = params.copy()
        
        # dictionary of steady state values        
        self.SS_dict    = self.get_steadyStateValues() 

        # optional keyword arguments for stochastic model!
        self.stochastic = kwargs.get('stochastic', False)
        self.seed       = kwargs.get('seed', None)
        self.z          = kwargs.get('z', 1.0)
        self.e          = kwargs.get('e', 0.0)
         
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

        # optional keyword arguments for welfare analysis
        self.utility    = kwargs.get('utility', None)
    
    def get_newCapital(self, k, zplus=1):
        """Function that takes end of period t's capital stock per 
        effective worker, k_t, and a contemporaneous technology 
        shock, z_{t+1}, and returns the end of period t+1's capital 
        stock per effective worker, kplus.

        k_{t+1} = \\frac{1}{(1 + g)(1 + n)z_{t+1}}[(1 - \delta)k_t + sk_t^{\alpha}]
        
        Inputs: 

            1. capital per effective worker, k.
            2. technology shock, z (default = 1).

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
        kplus = (1 / ((1 + g) * (1 + n) * zplus)) * ((1 - delta) * k + s * k**alpha) 

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

    def get_samplePath(self, T=None, tol=None):
        """Generate path of length T from current value of state;  or
        generate a sample path long enough such that the difference 
        between the steady state value and current value is less than 
        tol (useful for welfare comparisons).

        Inputs:
         
            1. T: length of sample path.
            2. tol: desired tolerance for convergence to steady state.

        Returns: an array of shape (T,3) representing a sample path of 
                 the Solow economy.
        
        """
        if T !=None and tol != None:
            assert "Either N or tol must be specified (but not both!)."
        elif T != None:
            path = np.zeros((T, 3))
            
            for t in xrange(T):
                path[t, 0] = self.k
                path[t, 1] = self.z
                path[t, 2] = self.e
                self.update()
            
        elif tol != None:
            # compute the steady state value 
            k_star = self.SS_dict['k_star']
            
            # initialize sample path and dist
            path = np.array([[self.k, self.z, self.e]])
            dist = np.abs(self.k - k_star)
            
            while dist > tol:
                self.update()
                step = np.array([[self.k, self.z, self.e]])
                path = np.vstack((path, step))
                dist = np.abs(self.k - k_star)   

        return path

    def get_marginalDist(self, k0, z0, e0, T=None, N=None):
        """Returns N draws of k_T, z_T, e_T starting from initial 
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

    def get_steadyStateValues(self): 
        """Returns a dictionary containing the steady state values of 
        capital, output, consumption, and investment per effective
        worker.
    
        """
        # extract params
        s     = self.param_dict['s']
        n     = self.param_dict['n']
        g     = self.param_dict['g']
        alpha = self.param_dict['alpha']
        delta = self.param_dict['delta']
    
        kstar = (s / ((1 + g) * (1 + n) - (1 -delta)))**(1 / (1 - alpha))
        ystar = kstar**alpha
        cstar = (1 - s) * ystar
        istar = s * ystar
        
        SSdict = {'k_star':kstar, 'y_star':ystar, 
                  'c_star':cstar, 'i_star':istar}
        return SSdict
    
    def get_numericSteadyState(self, k0=None):
        """Finds the steady state for the Solow economy using fsolve.
        Requires an initial guess for the steady state value of k, k0.
        Returns a list containing the steady state value of k.

        """
        # function to be optimized
        def solowSS(X):
            out = [self.get_newCapital(X[0]) - X[0]]
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
        capital_k = mp.diff(f=self.get_newCapital, 
                            x=(self.SS_dict['k_star']), n=(1))
        jacobian = np.array([capital_k], 'float')
        
        return [jacobian]

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

        # Create the plot
        ax = plt.subplot(111)
        output = ax.plot(grid, grid**alpha, 'r-')[0]
        actualInv = ax.plot(grid, s * grid**alpha, 'g-')[0]
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

    def plot_phaseSpace(self, gridmax, N=500, k0=None, T=None):
        """Phase diagram for the Solow model. 

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
        kLocus = ax.plot(grid, self.get_newCapital(grid), 'g-')[0]
        ax.plot(grid, grid, '--', color='k')
        
        # Demarcate the 45 degree line! 
        ax.add_patch(mpl.patches.Arc(xy=(0,0), 
                                     width=self.SS_dict['k_star'] / 2,
                                     height=self.SS_dict['k_star'] / 2, 
                                     theta1=0, theta2=45))
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
            phase_path = self.get_phasePath(k0, T)
            x_vec = phase_path[:, 0]
            y_vec = phase_path[:, 1]
            timePath = ax.plot(x_vec, y_vec)[0]

            return [ax, kLocus, timePath]
        else:
            return [ax, kLocus]
    
    def get_phasePath(self, k0=None, T=None):
        """Generate and plot a timepath for k of length T plotted in 
        phase space starting from k = k0.

        """
        self.k = k0
        n_even = 2 * int(T / 2)
        path = np.zeros(shape=(n_even, 2))
        
        for t in range(0, n_even, 2):
            path[t, 0] = self.k
            path[t + 1, 0] = self.k
            self.update()
            
        self.k = k0
        path[0, 1] = 0
        
        for t in range(1, n_even - 1, 2):
            self.update()
            path[t, 1] = self.k
            path[t + 1, 1] = self.k
        path[n_even - 1, 1] = self.k
        
        return path
        
    
