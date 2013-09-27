# dictionary of model parameters
params = {'DurationUnempG':1.5,
          'DurationUnempB':2.5,
          'corr':0.25,
          'UnempG':0.04,
          'UnempB':0.1,
          'DurationZG':8.0,
          'DurationZB':8.0,
          }

class Model(object):
    """Abstract class for solving and simulating macro models with heterogenous
    agents a la Krusell and Smith (1998, JPE).

    """
    def __init__(self, params):
        """Initializes a Model object with the following attributes:

        """
        # make sure to use copy!
        self.params = params.copy()

# Defining the transition matrix.
def get_transition_matrix(params):
    """Constructs the transition matrix."""

    # unemployment rates depend only on the aggregate productivity shock
    Unemp = np.array([[params['UnempG']], 
                      [params['UnempB']])

    # probability of remaining in 'Good/High' productivity state
    pZG = 1 - 1 / params['DurationZG']

    # probability of remaining in the 'Bad/Low' productivity state
    pZB = 1 - 1 / params['DurationZB']

    # matrix of transition probabilities for aggregate state
    PZ = np.array([[pZG, 1 - pZG], 
                   [1 - pZB, pZB]])

    # transition probabilities between employment states when aggregate productivity is high
    p22 = 1 - 1 / params['DurationUnempG']
    p21 = 1 - p22
    p11 = (((1 - params['UnempG']) - params['UnempG'] * p21) / 
           (1 - params['UnempG']))
    P11 = np.array([[p11, 1 - p11], 
                    [p21, p22]])

    # transition probabilities between employment states when aggregate productivity is low
    p22 = 1 - 1 / params['DurationUnempB']
    p21 = 1 - p22
    p11 = (((1 - params['UnempB']) - params['UnempB'] * p21) / 
           (1 - params['UnempB']))
    P00 = np.array([[p11, 1 - p11], 
                    [p21, p22]])

    # transition probabilities between employment states when aggregate productivity is high
    p22 = (1 + params['corr']) * p22
    p21 = 1 - p22
    p11 = (((1 - params['UnempB']) - params['UnempG'] * p21) / 
           (1 - params['UnempG']))
    P10 = np.array([[p11, 1 - p11], 
                    [p21, p22]])

    p22 = (1 - params['corr']) * (1 - 1 / params['DurationUnempG'])
    p21 = 1 - p22
    p11 = (((1 - params['UnempG']) - params['UnempB'] * p21) / 
           (1 - params['UnempB']))

    P01 = np.array([[p11,1 - p11], 
                    [p21, p22]])

    P = np.vstack((np.hstack((PZ[0, 0] * P11, PZ[0, 1] * P10)), 
                   np.hstack((PZ[1, 0] * P01, PZ[1, 1] * P00))))

    return P
