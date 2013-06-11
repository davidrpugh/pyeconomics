import scipy.io as spio
import matplotlib.pyplot as plt

class Model(object):
    """Generic class for working with the output of a Dynare model file in Python.

    TODO: 

        1) Pass the model file to MatLab/Octave directly using mlabwrap/Octave 
        magic commands.

    """
    
    def __init__(self, output):
        """Initializes a Model object from a Dynare .mat output file."""
        
        # Loads the matlab file
        self.output = spio.loadmat(output)
        
        # Generic information about the MatLab file
        self.header   = self.output['__header__']
        self.globals  = self.output['__globals__']
        self.version  = self.output['__version__']
        
        # Contains various informations about the model.
        self.M_       = ModelInformation(self)
        
        # Contains values of the options used by Dynare. 
        self.options_ = ModelOptions(self)
        
        # Contains results of the various computations.
        self.oo_      = ModelResults(self)
        
        # Contains the model impulse response functions (IRFs)
        self.irfs     = ImpulseResponse(self.oo_)

class ModelInformation(object):
    
    def __init__(self, model):
        """Initializes a ModelInformation object using the M_ data structure 
        from the model output file.

        """
        self.model = model
        
        # Strips off extra dimensions from M_ structure
        self.info_dict = {key:model.output['M_'][0,0][key] for key in 
                          model.output['M_'].dtype.names}
        
    def __getitem__(self, key):
        """Called to implement evaluation of self[key].""" 
        return self.info_dict[key]
        
class ModelOptions(object):
    
    def __init__(self, model):
        """Initializes a ModeOptions object using the options_ data structure 
        from the model output file. Should prove useful when I am ready to 
        construct a .mod file and call Dynare/Octave from within Python.

        """
        self.model = model

        # strips off extra dimensions from options_ structure
        self.options_dict = {key:model.output['options_'][0,0][key] for key in 
                             model.output['options_'].dtype.names}

    def __getitem__(self, key):
        """Called to implement evaluation of self[key].""" 
        return self.options_dict[key]

class ModelResults(object):
    
    def __init__(self, model):
        """Initializes a ModelResults object using the oo_ data structure from 
        the model output file.

        """
        self.model = model
        
        # strips off extra dimensions
        self.results_dict = {key:model.output['oo_'][0,0][key] for key in 
                             model.output['oo_'].dtype.names}
        
        # replace steady_state array with a dictionary
        self.results_dict['steady_state'] = self.get_steady_states(model)
    
    def __getitem__(self, key):
        """Called to implement evaluation of self[key].""" 
        return self.results_dict[key]
    
    def get_steady_states(self, model):
        """Links variable names with steady-state values in a Python dictionary.

        """
        SS_dict = {}
        
        # make sure to strip white space from var names!
        for i, var in enumerate(model.M_['endo_names']): 
            SS_dict[var.strip()] = model.output['oo_'][0,0]['steady_state'][i, 0]
            
        return SS_dict
    
class DecisionRule(object):
    """
    TODO: 
        1) Parse the dr data structure to recover the policy and transistion 
        functions.
        2) Write a plotting method for policy functions (see Dynare reference 
        manual pg. 40)

    """
    def __init__(self, results):
        """Initializes a DecisionRule object as a Python dict from a ModelResults
        object.

        """
        self.model = results.model

        # strips off extra dimensions
        self.dr_dict = {key:results['dr'][key][0,0] for key in 
                        results['dr'].dtype.names}
        
    def __getitem__(self, key):
        """Called to implement evaluation of self[key].""" 
        return self.dr_dict[key]

class ImpulseResponse(object):
    
    def __init__(self, results):
        """Initializes an ImpulseResponseFunction object as a Python dict from a
        ModelResults object."""
        self.model = results.model
        
        # strips off extra dimensions
        self.irfs_dict = {var: results['irfs'][var][0,0].flatten() for var in 
                          results['irfs'].dtype.names}
        
    def __getitem__(self, key):
        """Called to implement evaluation of self[key].""" 
        return self.irfs_dict[key]
    
    def plot(self, ax, variable, percentage=True, formatting=True):
        """Plots the impulse response functions (IRFs) for your a variable.

        Arguments:

            ax:         (Axes) Matplotlib Axes object on which to plot the IRF.
            variable:   (tuple) Variable names whose IRFs you wish to plot.
            percentage: (boolean) Whether or not you want to plot IRFs in terms
                        of percentage deviations from steady state. Default is 
                        True.
            formatting: (boolean) Whether or not you want to to include
                        additional formatting for the Axes object. Default is 
                        True.

        Returns: A list containing...
        
            ax:       Matplotlib Axes object.
            irf_plot: Matplotlib Line2D object.

        TODO:
        
            1) Incorporate LaTex variable names to make plots look better.
            2) Incorporate error bands for IRFs with order = 2 and 3.
            
        """

        # always pad with some zeros
        padding = np.zeros(self.model.options_['irf'][0,0] / 5.0)
        
        if percentage == True:
            # extract steady state value
            ss_val = self.model.oo_['steady_state'][variable.split('_')[0]]
            
            # compute percentage deviations
            data = np.append(padding, 100 * (self[variable] / ss_val))
            
            # label the y-axis accordingly
            ax.set_ylabel('% Deviations')
                             
        else:
            # by default IRFs are absolute deviations
            data = np.append(padding, self[variable])
            
            # label the y-axis accordingly
            ax.set_ylabel('Absolute deviations')
        
        # plots the IRF
        irf_plot, = ax.plot(data, label=variable)
        
        # optional formatting for the Axes object
        if formatting == True:
            ax.set_xlabel('Periods')
            ax.set_title('IRF for %s' %variable)
            ax.legend(loc='best', frameon=False)
            ax.grid()

        return [ax, irf_plot]

    def grid_plot(self, variables, nrows, ncols, sharey=True, percentage=True, 
                  formatting=True):
        """Plots the impulse response functions (IRFs) for several model 
        variables.

        Arguments:

            variables:  (tuple) Variable names whose IRFs you wish to plot.
            nrows:
            ncols:
            sharey:     (boolean) Whether or not all subplots should share a
                        common y-scale. Default is True.
            percentage: (boolean) Whether or not you want to plot IRFs in terms
                        of percentage deviations from steady state. Default is 
                        True.
            formatting: (boolean) Whether or not you want to to include
                        additional formatting for the Axes objects. Default is 
                        True.

        """
        # sanity check for args
        if nrows * ncols < len(variables):
            raise Exception, ('Wrong number of subplots! nrows * ncols should ' +
                              'be >= number of variables.')

        # create new Figure
        fig, axarr = plt.subplots(nrows, ncols, sharey=sharey, squeeze=False, 
                                  figsize=(8,8))
        
        # initialize a variable counter (will break if extra subplots!)
        i = 0
        
        for j in range(nrows):
            for k in range(ncols):
                tmp_var = variables[i] 
                tmp_ax  = axarr[j][k]
            
                # plot the IRF
                tmp_ax, tmp_irf = self.plot(tmp_ax, tmp_var, percentage, False)
            
                # only label x-axis for plots in last row
                if j == nrows - 1:
                    tmp_ax.set_xlabel('Periods')
                else:
                    tmp_ax.set_xlabel('')
                    
                # only label y-axis for plots in first column
                if k == 0 and percentage == True:
                    tmp_ax.set_ylabel('% Deviations')
                elif k == 0 and percentage == False:
                    tmp_ax.set_ylabel('Absolute Deviations')
                else:
                    tmp_ax.set_ylabel('')
                    
                # additional, optional, formatting
                if formatting == True:
                    tmp_ax.grid()
                    tmp_ax.legend(loc='best', frameon=False)
        
                # increment the variable counter
                i += 1
                
        return [fig, axarr]
            
