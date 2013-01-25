import numpy as np

def CRRA(C):
    """Constant relative risk aversion utility function: 

    C = \frac{C^{1 - \\theta} - 1}{1 - \\theta}

    where C is consumption per worker/capita (depending). The parameter 
    \\theta plays two roles:
    
        1. coefficient of relative risk aversion
        2. inverse of the inter-temporal elasticity of substitution

    For low values of \\theta (i.e., those less than 1), the agent is 
    nearly risk neutral (i.e., CRRA is almost linear!) and is quite
    willing to shift consumption across time to take advantage of small
    changes to interest/discount rates. For high values of \\theta 
    (i.e., those greater than 1) the opposite is true: the agent is 
    risk averse and is unwilling to shift consumption across time. For 
    \\theta == 1, CRRA is isomorphic to log utility.
    
    """
    # extract the params
    theta = params['theta']
    
    if theta != 1:
        return (C**(1 - theta) - 1) / (1 - theta)
    else:
        return np.log(C)
    
