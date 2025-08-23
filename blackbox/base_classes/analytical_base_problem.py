from abc import abstractmethod
import numpy as np
from .base_problem import BaseProblem

class AnalyticalProblem(BaseProblem):

    def __init__(self, lb, ub, negate):
        """
            Base class for all analytical problems
        """

        try:
            assert isinstance(lb, np.ndarray) and isinstance(ub, np.ndarray), "Lower and upper bounds must be numpy arrays"
            assert lb.ndim == 1 and ub.ndim == 1, "Lower and upper bounds must be one-dimensional numpy arrays"
            assert lb.shape == ub.shape, "Lower and upper bounds must have the same shape"
            assert np.all(lb < ub), "Lower bounds must be less than upper bounds"
            assert isinstance(negate, bool), "negate must be a boolean variable"
        
        except AssertionError as e:
            self._error(str(e))

        self.lb = lb
        self.ub = ub
        self.dim = lb.shape[0]
        self.negate = negate

    def __call__(self, x):
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
                containing the function value for each input sample
        """

        self.check_input(x)

        ndim = x.ndim

        if ndim == 1:
            x = x.reshape(1,-1)

        y = self.evaluate(x) 
        
        if self.negate:
            y = -y

        y = y.reshape(-1,1)
        
        if ndim == 1:
            y = y.reshape(-1,)
        
        return y

    @abstractmethod
    def evaluate(self):
        pass
