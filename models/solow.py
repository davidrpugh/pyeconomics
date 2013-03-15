import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import optimize, stats

class Model(object):
    
    def __init__(self, params=None, k=None, **kwargs):
        """A discrete time version of the classic Solow (1956) model of
        economic growth. 

        Summary of exogenous parameters of the Solow Model:
        
            g:     growth rate of technology
            n:     population growth rate
            delta: rate of capital depreciation
            alpha: capital share of output
            s:     savings rate
            rho:   (optional) persistence of technology shocks.
            sigma: (optional) standard deviation of technology shocks.
            
        Required attributes: 
        
            params: Either a dictionary of parameter values for the 
                    model, i.e. {'alpha':0.33,...}, or a valid ISO 
                    country code.  

                    If an ISO country code is specified, the model 
                    parameters will be calibrated as follows:

                        g:
                        n:
                        delta:
                        alpha:
                        s:
                        
            k:      a number representing the initial condition of 
                    the state variable (i.e., capital per effective 
                    worker.

        Optional keyword arguments:

            timing:     One of either 'discrete' or 'continuous'. 
                        Default is 'discrete'.

        """
        # initial value of the state variable
        self.k          = k
        
        # create the dictionary of parameter values (use copy!)
        self.param_dict = params.copy()

        # is the model discrete or continuous
        self.timing     = get.kwargs('timing', 'discrete')
        
        # dictionary of steady state values        
        self.SS_dict    = self.get_steadyStateValues() 
    
    def get_samplePath(self, T=None, tol=None, **kwargs):
        """Generate path of length T from current value of state;  or
        generate a sample path long enough such that the difference 
        between the steady state value and current value is less than 
        tol (useful for welfare comparisons).

        Optional kwargs: 

            'stochastic': If you wish to simulate a stochastic 
                          version of the Solow model, then you will
                          need set this flag to True.
            'seed':       a seed value for the random number 
                          generator (required if 'stochastic'=True).

        Inputs:
         
            T:   length of sample path.
            tol: desired tolerance for convergence to steady state.

        Returns: an array of shape (T,3) representing a sample path of 
                 the Solow economy.
        
        """
        ##### optional keyword arguments for stochastic model! #####

        # by default you are simulating a deterministic model...
        self.stochastic = kwargs.get('stochastic', False)
        self.seed       = kwargs.get('seed', None)

        # ...but want stochastic model to nest deterministic model.
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

        ##### generates the sample path #####

        # catch user errors
        if self.stochastic == True and tol != None:
            assert "tol must be 'None' if simulating stochastic model."
        elif T !=None and tol != None:
            assert "Either N or tol must be specified (but not both!)."

         # generate sample path of fixed length T 
        if T != None:
            path = np.zeros((int(T), 4))
            
            for t in xrange(int(T)):
                path[t, 0] = self.k # capital
                path[t, 1] = self.c # consumption
                path[t, 2] = self.z # technology shock
                path[t, 3] = self.e # disturbance
                self.update(policy, **kwargs)

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
                self.update(policy)
                dist = np.abs(self.k - k_star)   
                n_iter += 1
                
        return path

    def update(self, policy=None):
        """Update the state variable, k, the control, c, as well as 
        the technology shocks and their innovations (if stochastics 
        are turned on!)

        """
        # are we simulating a stochastic model?
        if self.stochastic == False:            
            # choose consumption 
            self.c = self.get_consumption()

            # update capital (note that z = zplus = 1!)
            self.k = self.get_newCapital(self.k, self.c, 1, 1)

        else:
            # need to update the innovation and technology
            self.e = self.eps.rvs()
            self.z = self.get_newTechShock(self.z, self.e)

    def get_consumption(self, policy):
        """Decision rule for choosing consumption.

        The following consumption policies are implemented:

            'solow': simple 'rule of thumb' decision rule. Consume a
                     fixed fraction of output.

        """
        # extract parameters
        s     = self.param_dict['s']
        alpha = self.param_dict['alpha']

        return (1 - s) * self.k**alpha

    def get_newTechShock(self, z, eplus):
        """Defines an AR(1) process representing technology shocks.

        Inputs:

            z: previous value of the technology shock, z.
            e: a draw from a N(0, sigma) RV.

        Returns: current value for the technology shock.
        
        """
        # extract params
        rho = self.param_dict['rho']

        # technology shock is log-normal
        zplus = self.z**rho * np.exp(eplus)

        return zplus
    
    def get_newCapital(self, k, c, z, zplus):
        """Function that maps the capital stock per effective worker 
        in period t, k_t, to the capital stock per effective worker in
        period k_{t+1}. In the deterministic case (i.e., when there are
        no technology shocks), the equation of motion is:

        k_{t+1} = \\frac{1}{(1 + g)(1 + n)}[(1 - \delta)k_t + k_t^{\alpha} - c_t]

        In the stochastic case, the equation of motion also needs to 
        take into account the values of technology shocks, z_t and 
        z_{t+1}:

        k_{t+1} = \\frac{1}{(1 + g)(1 + n)}\\frac{z_t}{z_{t+1}}[(1 - \delta)k_t + k_t^{\alpha} - c_t]

        The houshold's stock of captial per effective worker in period 
        t+1 depends on its choice of consumption in period t! Also, 
        note that the determinisitc model is a special case of the 
        stochastic model where z_t = z_{t+1} = 1.

        Inputs: 

            1. capital per effective worker, k.
            2. consumption per effective worker, c.
            3. technology shock, z.
            4. technology shock, zplus.

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
        kplus = (1 / ((1 + g) * (1 + n))) * ((1 - delta) * k + k**alpha - c) 

        return kplus
            
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
        
    
