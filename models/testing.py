import numpy as np
from scipy import integrate, linalg, optimize
import matplotlib.pyplot as plt

class SteadyState(object):
    """Abstract class representing the deterministic steady state of a model."""

    def __init__(self, model):
        """Initializes a SteadyState object with the following attributes:
        
            model: (object) An instance of the RamseyModel or SolowModel class.
            
        """
        self.model = model
        
    def set_functions(self, func_dict):
        """Modifies the model's steady_state_functions attribute.
        
        Arguments:
            
            func_dict: (dict) Dictionary of analytic function defining the 
                       model's steady state as functions of model parameters.
                       
        """  
        self.functions = func_dict
        
    def set_values(self):
        """Computes the steady state values using the dictionaries of analytic
        functions and model parameters.
                  
        """
        # store steady state values in a dictionary
        steady_state_values = {}
        
        # populate the dictionary of steady state values
        for key, func in self.functions.iteritems():
            steady_state_values[key] = func(self.model.params)
        
        self.values = steady_state_values
    
    def find_values(self, func, method, **kwargs):
        """Provides an interface to the various methods for finding the root of 
        a univariate or multivariate function in scipy.optimize.
    
        Arguments:
            
            method:   (str) For univariate functions method must be one of: 
                      'brentq', 'brenth', 'ridder', 'bisect', or 'newton'; for
                      multivariate function method must be one of: 'hybr', 'lm',
                      'broyden1', 'broyden2', 'anderson', 'linearmixing', 
                      'diagbroyden', 'excitingmixing', 'krylov'.
            
            **kwargs: (dict) Dictionary of method specific keyword arguments.
        
        """ 
        # list of valid multivariate root-finding methods
        multivariate_methods = ['hybr', 'lm', 'broyden1', 'broyden2', 
                                'anderson', 'linearmixing', 'diagbroyden', 
                                'excitingmixing', 'krylov']  
        # univariate methods     
        if method == 'brentq':
            res = optimize.brentq(func, args=(self.model.params,), **kwargs)
        elif method == 'brenth':
            res = optimize.brenth(func, args=(self.model.params,), **kwargs)
        elif method == 'ridder':
            res = optimize.ridder(func, args=(self.model.params,), **kwargs)
        elif method == 'bisect':
            res = optimize.bisect(func, args=(self.model.params,), **kwargs)
        elif method == 'newton':
            res = optimize.newton(func, args=(self.model.params,), **kwargs)
        
        # multivariate methods are handled by optimize.root
        elif method in multivariate_methods:
            res = optimize.root(func, args=(self.model.params,), **kwargs)
        else:
            raise ValueError, 'Unrecognized method!'
    
        return res
        
    def check_stability(self, eval_jacobian):
        """Computes the eigenvalues and eigenvectors of the Jacobian matrix
        evaluated at the deterministic steady state. In order for the steady
        state to be dynamically stable we need to have:

            1. same number of stable eigenvalues as pre-determined 
               variables (i.e., state variables)
            2. same number of unstable eigenvalues as control (i.e., 
               jump variables).

        Arguments:
            
            eval_jacobian: (array-like) Model jacobian evaluated at steady state.
            
        Returns: A list containing...

            eigenvalues:  The eigenvalues of the Jacobian matrix.
            eigenvectors: The eigenvectors of the Jacobian matrix.  
            
        """
        eigenvalues, eigenvectors = linalg.eig(eval_jacobian) 
        return [eigenvalues, eigenvectors]

class SolowModel(object):
    """Base class for working with Solow Growth models."""
    
    def __init__(self, params, capital, jacobian):
        """Initializes a SolowModel object with the following attributes:
        
            params:      (dict) Dictionary of model parameters.
            
            capital:     (callable) Equation of motion for capital (per 
                         person/effective person).
                                
            jacobian:    (callable) Function returning the model's jacobian 
                         matrix of partial derivatives.
            
        """
        # initialize model attributes
        self.params       = params
        self.capital      = capital
        self.jacobian     = jacobian

        # initialize a SteadyState object
        self.steady_state = SteadyState(self)
        
        # create and instance of the scipy.integrate.ode class      
        self.simulator    = integrate.ode(capital, jacobian)
    
    def update_model_parameters(self, new_params):
        """Updates the model's parameter dictionary."""
        self.params = new_params.copy()
        
    def simulate(self, k0, T, h, integrator, step=False, relax=False, **kwargs):
        """Simulates model trajectories.
        
        Arguments:
                
            k0:         (float) Initial conditions for capital (per 
                        person/effective person).
                      
            T:          (int) Length of desired trajectory. 
                                            
            h:          (float) Step-size for computing the solution.
            
            integrator: (str) Must be one of 'vode', 'lsoda', 'dopri5', 'dop85'. 
                        See documentation for integrate.ode for details.
            
            step:       (boolean) The following integrators support step: 
                        'vode', 'zvode', 'lsoda'. Default is False. 
                         
            relax:      (boolean) The following integrators support run_relax: 
                        'vode', 'zvode', 'lsoda'. Default is False. 
                     
            **kwargs:   (dict) Dictionary of method specific keyword args.
                
        Returns: 
                     
           solution: (array-like) Trajectory for capital (per person/effective
                     person).
               
        """        
        # select the integrator
        self.simulator.set_integrator(integrator, **kwargs)
        
        # pass the model parameters as additional args
        self.simulator.set_f_params(self.params)
        self.simulator.set_jac_params(self.params)
        
        # set the initial condition
        self.simulator.set_initial_value(k0, 0)
        
        # create a storage container for the trajectory
        solution = []
        
        # generate a trajectory of length T
        while self.simulator.successful() and self.simulator.t <= T: 
            solution.append([self.simulator.t, self.simulator.y])               
            self.simulator.integrate(self.simulator.t + h, step, relax)       
    
        return np.array(solution)
                               
    def plot_phase_diagram(self, N=1000):
        """Generates a plot of the phase diagram for the Solow model.
        
        Arguments:
            
            N: (int, optional) Number of points to plot . Default is 1000. 
            
        Returns: A list containing...
        
            ax:   (object) Axes object representing the plot.
            line: (object) Line2D object representing the k-dot locus.
        
        """
        # create a new figure and subplot
        ax     = plt.subplot(111)

        # use value of k_star to anchor the plot
        k_star = self.steady_state_values['k_star']
        grid   = np.linspace(0, 2 * k_star, N)

        # plot the evolution of capital
        data   = self.capital(grid, self.params)
        line   = ax.plot(grid, data)

        # adjust the ylims
        ax.set_ylim(-np.max(np.abs(data)), np.max(np.abs(data)))

        # remove the right and top spines
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')

        # hide the top and right ticks
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()

        # center the bottom spine on the data and demarcate the steady state
        ax.spines['bottom'].set_position(('data', 0))
        ax.set_xticks([k_star])
        ax.set_xticklabels([r'$k^*$'])

        # label the vertical axis
        ax.set_ylabel(r'$\dot{k}$', rotation='horizontal', fontsize=15)

        # provide a title
        ax.set_title('Phase diagram for $k$ in the Solow model', fontsize=15)

        return [ax, line]
        
    def plot_solow_diagram(self, gridmax, N=1000):
        """Generates the classic Solow diagram.
        
        Arguments:
            
            gridmax: (float) Maximum value for capital per effective worker.
            N:       (int, optional) Number of points to plot. Default is 1000. 
            
        Returns: A list containing...
        
            ax:                (object) Axes object representing the plot.
            output:            (object) Line2D object representing output.
            actual_invest:     (object) Line2D object representing actual 
                               investment.
            break_even_invest: (object) Line2D object representing break-even
                               investment.
        
        """
        # extract params
        n     = self.params['n']
        g     = self.params['g']
        s     = self.params['s']
        alpha = self.params['alpha']
        delta = self.params['delta']
        
        # create a new figure and subplot
        ax = plt.subplot(111)

        # grid of values for capital per effective worker
        grid   = np.linspace(0, gridmax, N)
          
        # plot output, actual and break even investment             
        output            = ax.plot(grid, grid**alpha, 'r-')[0]
        actual_invest     = ax.plot(grid, s * grid**alpha, 'g-')[0]
        break_even_invest = ax.plot(grid, (n + g + delta) * grid, 'b-')[0]

        # remove the right and top spines
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')

        # hide the top and right ticks
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()

        # axes, labels, title, legend, etc
        ax.set_xlabel('Capital per effective worker, $k$', fontsize=15)
        ax.set_xlim(0, gridmax)
        ax.set_ylabel('$y$, $i$, and $c$', rotation='horizontal', 
                      fontsize=15)
        ax.set_title('Classic Solow Diagram\n' + 
                     '$s=%g, n=%g, g=%g, \delta=%g$' %(s,n,g,delta), 
                     fontsize=20)
        
        return [ax, output, actual_invest, break_even_invest]    
  
class RamseyModel(object):
    """Base class for working with continuous-time growth models."""
    
    def __init__(self, params, capital, consumption, jacobian, utility):
        """Initializes a RamseyModel object with the following attributes:
        
            params:   (dict) Dictionary of model parameters.
            
            capital:  (callable) Equation of motion for capital (per 
                      person/effective person).
                                
            euler:    (callable) Equation describing consumption (per 
                      person/effective person).
                                
            jacobian: (callable) Function returning the model's jacobian matrix 
                      of partial derivatives.
                                
            utility:  (callable) Utility function representing household 
                      preferences over consumption.
            
        """
        self.params      = params
        self.capital     = capital
        self.consumption = consumption
        self.jacobian    = jacobian
        self.utility     = utility
                   
        # initialize model steady state
        self.steady_state = SteadyState(self)
        
        # create an instance of the scipy.integrate.ode class
        f = lambda t, y, params: [self.capital(y, params), 
                                  self.consumption(y, params)]
        jac = lambda t, y, params: self.jacobian(y, params)        
        self.simulator = integrate.ode(f, jac)
                
    def update_model_parameters(self, new_params):
        """Updates the model's parameter dictionary."""
        self.params = new_params.copy()
        
    def simulate(self, y0, T=None, tol=None, method=None, h=None, **kwargs):
        """Simulates model trajectories.
        
        Arguments:
                
            y0:       (array-like) Initial conditions for capital and
                      consumption (per person/effective person).
            T:        (int) Length of desired trajectory. Only one of T or tol
                      should be specified.
            tol:      (float) Desired convergence tolerance. Simulation will 
                      stop if the trajectory comes within tol of steady state.
                      Only one of T or tol should be specified.
            h:   
            method:
            **kwargs: (dict) Dictionary of method specific keyword args.
                
        Returns: 
                     
        
        """        
        # select the integrator
        self.simulator.set_integrator(method, **kwargs)
        
        # pass the model parameters as additional args
        self.simulator.set_f_params(self.params)
        self.simulator.set_jac_params(self.params)
        
        # set the initial condition
        self.simulator.set_initial_value(k0, 0)
        
        # create a storage container for the trajectory
        solution = [[self.simulator.t, self.simulator.y]]
        
        # generate a trajectory of length T
        if T != None:
            for i in range(T):
                self.simulator.integrate(self.simulator.t + h)
                if not self.simulator.successful():
                    break
                else:
                    solution.append([self.simulator.t, self.simulator.y])               
        
        # iterate until within some tolerance of steady state
        elif tol != None:
            # ordering should be k_star, c_star
            steady_state = [self.steady_state.values['k_star'],
                            self.steady_state.values['c_star']]
            dist = np.max(np.abs(self.simulator.y - steady_state))
            
            while self.simulator.successful() and dist < tol:
                solution.append([self.simulator.t, self.simulator.y]) 
                dist = np.max(np.abs(self.simulator.y - steady_state))
                
        else:
            raise ValueError, 'Either T or tol must be provided!'        
    
        return np.array(solution)
                
