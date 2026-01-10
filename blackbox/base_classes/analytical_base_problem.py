from abc import abstractmethod
import numpy as np
from .base_problem import BaseProblem

class SingleOutputAnalyticalProblem(BaseProblem):
    """
        Base class for all single output analytical problems
    """

    def __init__(self, num_inputs: int, negate: bool):
        
        assert isinstance(num_inputs, int) and num_inputs > 0, "num_inputs must be a positive integer"
        assert isinstance(negate, bool), "negate must be a boolean variable"

        self.num_inputs = num_inputs
        self.num_outputs = 1
        self.negate = negate

    def __call__(self, x: np.ndarray) -> np.ndarray:
        """
            Evalaute the function for given input x

            Parameters
            ----------
            x: np.ndarray
                1D/2D numpy array of shape (dim,) or (n_samples,dim)

            Returns
            -------
            y: np.ndarray
                1D/2D numpy array of shape (n_samples,1) or (1,) 
                containing the required output value for each input sample
        """

        self._check_input(x)

        y = self._evaluate(np.atleast_2d(x)) # evaluate the function
        
        y = y.reshape(-1,self.num_outputs) # ensure output is 2D

        if self.negate:
            y = -y

        if x.ndim == 1:
            y = y.reshape(-1,)
        
        return y
    
    def set_bounds(self, bounds: tuple):
        """
            Method to change the default bounds of the problem

            Parameters
            ----------
            bounds: tuple
                a tuple containing two 1D numpy array - first one is lower bound and second
                entry is upper bound
        """

        assert isinstance(bounds, tuple) and len(bounds) == 2, "bounds must be a tuple with two entries"
        for i, bound in enumerate(bounds):
            assert isinstance(bound, np.ndarray) and bound.ndim == 1, f"entry {i+1} in bounds must be a 1D numpy array"
            assert bound.shape[0] == self.num_inputs, f"size of entry {i+1} in bounds is not same as number of inputs"
        assert np.all(bounds[0] < bounds[1]), "Lower bound must be less than upper bound for each input"

        self.bounds = bounds

    def _check_input(self, x: np.ndarray):
        """
            Method to validate the input before evaluation

            Parameters
            ----------
            x: np.ndarray
                Input sample(s) to be evaluated, can be either a 1D or 2D numpy array
        """

        x = np.atleast_2d(x)

        assert x.shape[1] == self.bounds[0].shape[0], "Input dimension must match the problem dimension"
        assert np.all(x >= self.bounds[0]) and np.all(x <= self.bounds[1]), "Input values must be within the bounds"
    
    @abstractmethod
    def _evaluate(self, x: np.ndarray) -> np.ndarray:
        pass
