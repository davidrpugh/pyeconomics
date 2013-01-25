import numpy as np

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
