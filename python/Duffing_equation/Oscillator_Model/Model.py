#
# Model class
#

import numpy as np

class Model:
    """A Model for open and closing of ion channels
    Parameters
    ----------
    gamma: float
        damping parameter; default = 1
    driving_amplitude: float
        default = 1
    driving_freq: float
        default = 0
        if driving_freq > 0, the driving term is periodic
        if driving_freq == 0, the driving term is stochastic
            and uniformly distributed betwen plus/minus driving_amplitude
    initial_Condition: array of floats
        default = [0,0]
        initial conditions for position and velocity
    """
    def __init__(self, gamma=1, driving_amplitude=1, driving_freq=0,  initial_Condition=[0,0], potential_params=0):
        self.__gamma = gamma
        self.__IC = initial_Condition
        self.__dim = len(initial_Condition)
        self.__driving_amplitude = driving_amplitude
        self.__driving_freq = driving_freq
        self.__potential_params = potential_params

    def dV_potential(self, x):
        """
        Specify a potential 
        returns first derivative of potential
        """
        fac = 0.5
        if self.__potential_params == 0:
            return fac * (-x + x*x*x)
        else:
            # return func(x)
            return fac * (-x + x*x*x + self.__potential_params[2]) * (1 + self.__potential_params[0] * np.cos(self.__potential_params[1]*x))
        
    def drivingterm(self, t):
        """
        Specify a driving term
        """
        if self.__driving_freq > 0:
            #print('periodic driving term')
            return self.__driving_amplitude * np.cos(self.__driving_freq*t)
        else:
            #print('stochastic driving term')
            #return random.uniform(-1,1)
            print('Error in driving term')
            return False
    
     
    @property
    def gamma(self):
        """
        Friction coefficient gamma
        """
        return self.__gamma
    
    @property
    def driving_amplitude(self):
        """
        Driving amplitude
        """
        return self.__driving_amplitude
        
    @property
    def driving_freq(self):
        """
        Driving frequency
        """
        return self.__driving_freq
    
    @property
    def IC(self):
        """
        Initial condition
        """
        return self.__IC
    
    @property
    def dim(self):
        """
        Returns the number of peripheral compartments.
        """
        return self.__dim
    
    @property
    def potential_params(self):
        return self.__potential_params

