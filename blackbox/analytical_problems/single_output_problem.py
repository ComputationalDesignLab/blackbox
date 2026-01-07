import numpy as np
from ..base_classes.analytical_base_problem import SingleOutputAnalyticalProblem

class Ackley(SingleOutputAnalyticalProblem):

    def __init__(self, num_inputs: int = 10, negate: bool = False):
        """
            Class for defining the ackley function. More details about this
            problem can be found here: https://www.sfu.ca/~ssurjano/ackley.html

            Default range of bounds for the design variable is [-32.768,32.768]^10

            Parameters
            ----------
            num_inputs: int
                number of input dimensions, default = 10
            negate: bool
                negate the values before returning
        """

        super().__init__(num_inputs, negate)

        self.bounds = (np.array([-32.768]*self.num_inputs), np.array([32.768]*self.num_inputs))
        self._x_opt = np.array([0.0]*self.num_inputs)
        self._y_opt = 0.0

    def _evaluate(self, x:np.ndarray) -> np.ndarray:

        a = 20.0
        b = 0.2
        c = 2.0 * np.pi

        return -a * np.exp( -b * np.linalg.norm(x, axis=-1) / np.sqrt(self.num_inputs) ) \
            - np.exp( np.mean(np.cos(c*x), axis=-1) ) + a + np.e


class Levy(SingleOutputAnalyticalProblem):

    def __init__(self, num_inputs: int = 10, negate: bool = False):
        """
            Class for defining the levy function. More details about this
            problem can be found here: https://www.sfu.ca/~ssurjano/levy.html

            Parameters
            ----------
            num_inputs: int
                number of input dimensions, default = 10
            negate: bool
                negate the values before returning
        """

        super().__init__(num_inputs, negate)

        self.bounds = (np.array([-10]*num_inputs), np.array([10]*num_inputs))
        self._x_opt = np.array([1.0]*self.num_inputs)
        self._y_opt = 0.0

    def _evaluate(self, x:np.ndarray) -> np.ndarray:

        w = 1.0 + (x - 1.0) / 4.0

        part1 = np.sin(np.pi*w[:,0])**2
        
        part2 = np.sum( (w[:,:-1] - 1.0)**2 * (1.0 + 10.0*np.sin(np.pi*w[:,:-1] + 1.0)**2), axis=-1 )

        part3 = (w[:,-1] - 1.0)**2 * (1.0 + np.sin(2*np.pi*w[:,-1])**2 )

        return part1 + part2 + part3


class Rastrigin(SingleOutputAnalyticalProblem):

    def __init__(self, num_inputs: int = 10, negate: bool = False):
        """
            Class for defining the rastrigin function. More details about this
            problem can be found here: https://www.sfu.ca/~ssurjano/rastrigin.html

            Parameters
            ----------
            num_inputs: int
                number of input dimensions, default = 10
            negate: bool
                negate the values before returning
        """

        super().__init__(num_inputs, negate)

        self.bounds = (np.array([-5.12]*num_inputs), np.array([5.12]*num_inputs))
        self._x_opt = np.array([0.0]*self.num_inputs)
        self._y_opt = 0.0

    def _evaluate(self, x:np.ndarray) -> np.ndarray:

        return 10.0*self.num_inputs + np.sum(x**2 - 10.0*np.cos(2.0*np.pi*x), axis=-1)


class Hartmann(SingleOutputAnalyticalProblem):

    def __init__(self, num_inputs: int = 6, negate: bool = False):
        """
            Class for defining the hartmann function. More details about this
            problem can be found here:
            
            Hartmann 3D - https://www.sfu.ca/~ssurjano/hart3.html
            Hartmann 4D - https://www.sfu.ca/~ssurjano/hart4.html
            Hartmann 6D - https://www.sfu.ca/~ssurjano/hart6.html

            Parameters
            ----------
            num_inputs: int
                number of input dimensions, default = 6
            negate: bool
                negate the values before returning
        """

        super().__init__(num_inputs, negate)

        self.bounds = (np.zeros(self.num_inputs), np.ones(self.num_inputs))

        assert self.num_inputs in [3,4,6], "Hartmann function is only defined for 3, 4 or 6 dimensions"

        if self.num_inputs == 3:

            self.A = np.array([[3.0, 10, 30], [0.1, 10, 35], [3.0, 10, 30], [0.1, 10, 35]])

            self.P = np.array([[3689, 1170, 2673],
                               [4699, 4387, 7470],
                               [1091, 8732, 5547],
                               [381, 5743, 8828.0]])
            
            self._x_opt = np.array([0.114614, 0.555649, 0.852547])
            self._y_opt = -3.86278

        elif self.num_inputs == 4:

            self.A = np.array([[10, 3, 17, 3.5],
                               [0.05, 10, 17, 0.1],
                               [3, 3.5, 1.7, 10],
                               [17, 8, 0.05, 10],])
            
            self.P = np.array([[1312, 1696, 5569, 124.0],
                               [2329, 4135, 8307, 3736],
                               [2348, 1451, 3522, 2883],
                               [4047, 8828, 8732, 5743]])
            
            self._x_opt = np.array([0.1873, 0.1906, 0.5566, 0.2647])
            self._y_opt = -3.134353

        elif self.num_inputs == 6:

            self.A = np.array([[10, 3, 17, 3.5, 1.7, 8],
                               [0.05, 10, 17, 0.1, 8, 14],
                               [3, 3.5, 1.7, 10, 17, 8],
                               [17, 8, 0.05, 10, 0.1, 14]])

            self.P = np.array([[1312, 1696, 5569, 124, 8283, 5886],
                               [2329, 4135, 8307, 3736, 1004, 9991],
                               [2348, 1451, 3522, 2883, 3047, 6650],
                               [4047, 8828, 8732, 5743, 1091, 381]])
            
            self._x_opt = np.array([0.20169, 0.150011, 0.476874, 0.275332, 0.311652, 0.6573])
            self._y_opt = -3.322368

        self.alpha = np.array([[1.0], [1.2], [3.0], [3.2]])

    def _evaluate(self, x:np.ndarray) -> np.ndarray:

        innersum = np.sum(self.A * (np.expand_dims(x,-2) - 1e-4*self.P)**2, axis=-1)

        y = -np.matmul(np.exp(-innersum), self.alpha)

        if self.num_inputs == 4:
            y = (1.1 + y)/0.839

        return y


class Branin(SingleOutputAnalyticalProblem):

    def __init__(self, negate: bool = False):
        """
            Class for defining the branin function

            Parameters
            ----------
            negate: bool
                negate the values before returning
        """

        super().__init__(num_inputs=2, negate=negate)

        self.bounds = (np.array([-5.0, 0.0]), np.array([10.0, 15.0]))
        self._x_opt = np.array([-np.pi, 12.275]) # TO DO: can add other optima as well
        self._y_opt = 0.397887
    
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


class ModifiedBranin(SingleOutputAnalyticalProblem):

    def __init__(self, negate: bool = False):
        """
            Class for defining the modified branin function

            Parameters
            ----------
            negate: bool
                negate the values before returning
        """

        super().__init__(num_inputs=2, negate=negate)

        self.bounds = (np.array([-5.0, 0.0]), np.array([10.0, 15.0]))
        self._x_opt = np.array([-3.6893, 13.6301]) # TO DO: can add other optima as well
        self._y_opt = -16.64402

    def _evaluate(self, x:np.ndarray) -> np.ndarray:

        a = 1.0
        b = 5.1 / (4.0 * np.pi**2)
        c = 5.0 / np.pi
        r = 6.0
        s = 10.0
        t = 1.0 / (8.0 * np.pi)

        x1 = x[:,0]
        x2 = x[:,1]

        return a * (x2 - b * x1**2 + c * x1 - r)**2 + s * (1 - t) * np.cos(x1) + s + 5*x1
