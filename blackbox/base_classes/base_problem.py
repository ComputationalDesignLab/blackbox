from abc import abstractmethod, ABC
import numpy as np
from ..msg import print_msg

class BaseProblem(ABC):
    """
        An abstract base problem class which is inherited by all other problem classes.
    """

    @abstractmethod
    def __call__(self):
        
        pass

    def check_input(self, x: np.ndarray):
        """
            Method to validate the input before evaluation

            Parameters
            ----------
            x: np.ndarray
                Input sample(s) to be evaluated, can be either a 1D or 2D numpy array
        """

        x = np.atleast_2d(x)

        try:
            assert x.shape[1] == self.in_dim, "Input dimension must match the problem dimension"
            assert np.all(x >= self.lb) and np.all(x <= self.ub), "Input values must be within the bounds defined by lb and ub"
            
        except AssertionError as e:
            print_msg(str(e))
