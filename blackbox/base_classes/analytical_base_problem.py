from abc import abstractmethod
import numpy as np
from .base_problem import BaseProblem

class SingleOutputAnalyticalProblem(BaseProblem):

    def __init__(self, num_inputs: int, negate: bool):
        """
            Base class for all single output analytical problems
        """

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
    
    def set_bounds(self, bounds: np.ndarray):
        """
            Method to set the bounds or change the default bounds of the problem

            Parameter
            ---------
            bounds: np.ndarray
                2D numpy array of shape (2,num_inputs) containing the new lower and upper bounds.
                First and second row correspond to lower and upper bounds, respectively
        """

        assert isinstance(bounds, np.ndarray), "bounds must be a numpy array"
        assert bounds.ndim == 2, "bounds must be a 2D numpy array"
        assert bounds.shape == (2,self.num_inputs), f"bounds must have shape (2,{self.num_inputs})"
        assert np.all(bounds[0,:] < bounds[1,:]), "Lower bounds must be less than upper bounds"

        self.lb = bounds[0,:]
        self.ub = bounds[1,:]

    def _check_input(self, x: np.ndarray):
        """
            Method to validate the input before evaluation

            Parameters
            ----------
            x: np.ndarray
                Input sample(s) to be evaluated, can be either a 1D or 2D numpy array
        """

        x = np.atleast_2d(x)

        assert x.shape[1] == self.lb.shape[0], "Input dimension must match the problem dimension"
        assert np.all(x >= self.lb) and np.all(x <= self.ub), "Input values must be within the bounds defined by lb and ub"
    
    @abstractmethod
    def _evaluate(self, x: np.ndarray) -> np.ndarray:
        pass
