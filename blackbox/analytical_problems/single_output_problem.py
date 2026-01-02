import numpy as np
from ..msg import print_msg
from abc import abstractmethod

class SingleOutputAnalyticalProblem():

    def __init__(self, num_inputs: int, negate: bool):
        """
            Base class for all single output analytical problems
        """

        try:
            assert isinstance(num_inputs, int) and num_inputs > 0, "num_inputs must be a positive integer"
            assert isinstance(negate, bool), "negate must be a boolean variable"
        
        except AssertionError as e:
            print_msg(str(e))

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

        try:
            assert isinstance(bounds, np.ndarray), "bounds must be a numpy array"
            assert bounds.ndim == 2, "bounds must be a 2D numpy array"
            assert bounds.shape == (2,self.num_inputs), f"bounds must have shape (2,{self.num_inputs})"
            assert np.all(bounds[0,:] < bounds[1,:]), "Lower bounds must be less than upper bounds"
        
        except AssertionError as e:
            print_msg(str(e))

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

        try:
            assert x.shape[1] == self.lb.shape[0], "Input dimension must match the problem dimension"
            assert np.all(x >= self.lb) and np.all(x <= self.ub), "Input values must be within the bounds defined by lb and ub"
            
        except AssertionError as e:
            print_msg(str(e))
    
    @abstractmethod
    def _evaluate(self, x: np.ndarray) -> np.ndarray:
        pass


class Ackley(SingleOutputAnalyticalProblem):

    def __init__(self, num_inputs: int = 10, negate=False):
        """
            Class for defining the ackley function

            Default range of bounds for the design variable is [-32.768,32.768]^10

            Parameters
            ----------
            num_inputs: int
                number of input dimensions, default = 10
            negate: bool
                negate the values before returning
        """

        super().__init__(num_inputs, negate)

        self.lb = np.array([-32.768]*self.num_inputs)
        self.ub = np.array([32.768]*self.num_inputs)

    def _evaluate(self, x:np.ndarray) -> np.ndarray:

        a = 20.0
        b = 0.2
        c = 2.0 * np.pi

        return -a * np.exp( -b * np.linalg.norm(x, axis=-1) / np.sqrt(self.num_inputs) ) \
            - np.exp( np.mean(np.cos(c*x), axis=-1) ) + a + np.e


class Levy(SingleOutputAnalyticalProblem):

    def __init__(self, num_inputs: int = 10, negate=False):
        """
            Class for defining the levy function

            Parameters
            ----------
            num_inputs: int
                number of input dimensions, default = 10
            negate: bool
                negate the values before returning
        """

        super().__init__(num_inputs, negate)

        self.lb = np.array([-10]*num_inputs)
        self.ub = np.array([10]*num_inputs)

    def _evaluate(self, x:np.ndarray) -> np.ndarray:

        w = 1.0 + (x - 1.0) / 4.0

        part1 = np.sin(np.pi*w[:,0])**2
        
        part2 = np.sum( (w[:,:-1] - 1.0)**2 * (1.0 + 10.0*np.sin(np.pi*w[:,:-1] + 1.0)**2), axis=-1 )

        part3 = (w[:,-1] - 1.0)**2 * (1.0 + np.sin(2*np.pi*w[:,-1])**2 )

        return part1 + part2 + part3


class Rastrigin(SingleOutputAnalyticalProblem):

    def __init__(self, num_inputs: int = 10, negate=False):
        """
            Class for defining the rastrigin function

            Parameters
            ----------
            num_inputs: int
                number of input dimensions, default = 10
            negate: bool
                negate the values before returning
        """

        super().__init__(num_inputs,negate)

        self.lb = np.array([-5.12]*num_inputs)
        self.ub = np.array([5.12]*num_inputs)

    def _evaluate(self, x:np.ndarray) -> np.ndarray:

        return 10.0*self.in_dim + np.sum(x**2 - 10.0*np.cos(2.0*np.pi*x), axis=-1)


class Hartmann(SingleOutputAnalyticalProblem):

    def __init__(self, num_inputs: int = 6, negate: bool = False):
        """
            Class for defining the hartmann function

            Parameters
            ----------
            num_inputs: int
                number of input dimensions, default = 6
            negate: bool
                negate the values before returning
        """

        super().__init__(num_inputs, negate)

        self.lb = np.zeros(self.num_inputs)
        self.ub = np.ones(self.num_inputs)

        try:
            assert self.num_inputs in [3,4,6], "Hartmann function is only defined for 3, 4 or 6 dimensions"
        
        except AssertionError as e:
            print_msg(str(e))

        if self.num_inputs == 3:

            self.A = np.array([[3.0, 10, 30], [0.1, 10, 35], [3.0, 10, 30], [0.1, 10, 35]])

            self.P = np.array([[3689, 1170, 2673],
                               [4699, 4387, 7470],
                               [1091, 8732, 5547],
                               [381, 5743, 8828.0]])

        elif self.num_inputs == 4:

            self.A = np.array([[10, 3, 17, 3.5],
                               [0.05, 10, 17, 0.1],
                               [3, 3.5, 1.7, 10],
                               [17, 8, 0.05, 10],])
            
            self.P = np.array([[1312, 1696, 5569, 124],
                               [2329, 4135, 8307, 3736],
                               [2348, 1451, 3522, 2883],
                               [4047, 8828, 8732, 5743]])

        elif self.num_inputs == 6:

            self.A = np.array([[10, 3, 17, 3.5, 1.7, 8],
                               [0.05, 10, 17, 0.1, 8, 14],
                               [3, 3.5, 1.7, 10, 17, 8],
                               [17, 8, 0.05, 10, 0.1, 14]])

            self.P = np.array([[1312, 1696, 5569, 124, 8283, 5886],
                               [2329, 4135, 8307, 3736, 1004, 9991],
                               [2348, 1451, 3522, 2883, 3047, 6650],
                               [4047, 8828, 8732, 5743, 1091, 381]])

        self.alpha = np.array([[1], [1.2], [3.0], [3.2]])

    def _evaluate(self, x:np.ndarray) -> np.ndarray:

        innersum = np.sum(self.A * (np.expand_dims(x,-2) - 1e-4*self.P)**2, axis=-1)

        y = -np.matmul(np.exp(-innersum), self.alpha)

        if self.num_inputs == 4:
            y = (1.1 + y)/0.839

        return y


class Branin(SingleOutputAnalyticalProblem):

    def __init__(self, negate=False):
        """
            Class for defining the branin function

            Parameters
            ----------
            negate: bool
                negate the values before returning
        """

        super().__init__(num_inputs=2, negate=negate)

        self.lb = np.array([-5.0, 0.0])
        self.ub = np.array([10.0, 15.0])
    
    def _evaluate(self, x:np.ndarray) -> np.ndarray:

        a = 1.0
        b = 5.1 / (4.0 * np.pi**2)
        c = 5.0 / np.pi
        r = 6.0
        s = 10.0
        t = 1.0 / (8.0 * np.pi)

        x1 = x[:,0]
        x2 = x[:,1]

        return a * (x2 - b * x1**2 + c * x1 - r)**2 + s * (1 - t) * np.cos(x1) + s
