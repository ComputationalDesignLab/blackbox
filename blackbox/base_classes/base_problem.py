from abc import abstractmethod, ABC
import numpy as np

class BaseFunction(ABC):
    """
        An abstract base function class which is inherited by all other function classes.
    """

    @abstractmethod
    def __call__(self):
        
        pass

    def checkInput(self, x):
        """
            Method to validate the input before evaluation

            Parameters
            ----------
            x: np.ndarray
                Input sample(s) to be evaluated, can be either a 1D or 2D numpy array
        """

        x = np.atleast_2d(x)

        assert x.shape[1] == self.dim, "Input dimension must match the problem dimension"

        assert np.all(x >= self.lb) and np.all(x <= self.ub), "Input values must be within the bounds defined by lb and ub"
