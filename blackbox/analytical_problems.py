import numpy as np
from .base_classes import AnalyticalProblem

class Ackley(AnalyticalProblem):

    def __init__(self, lb, ub, negate=False, a=20.0, b=0.2, c=2.0*np.pi):
        """
            Class for defining the ackley function

            Default range of bounds for the design variable is [-32.768,32.768]^10

            Parameters
            ----------
            lb: np.ndarray
                Lower bounds of the function
            ub: np.ndarray
                Upper bounds of the function
            negate: bool
                negate the values before returning
            a: float, optional
                Parameter a of the ackley function, default is 20.0
            b: float, optional
                Parameter b of the ackley function, default is 0.2
            c: float, optional
                Parameter c of the ackley function, default is 2.0 * pi
        """

        super().__init__(lb,ub,negate)

        self.a = a
        self.b = b
        self.c = c

    def evaluate(self, x):

        self.check_input(x)

        return -self.a * np.exp( -self.b * np.linalg.norm(x, axis=-1) / np.sqrt(self.dim) ) \
            - np.exp( np.mean(np.cos(self.c*x), axis=-1) ) + self.a + np.e


class Levy(AnalyticalProblem):

    def __init__(self, lb, ub, negate=False):
        """
            Class for defining the levy function

            Parameters
            ----------
            lb: np.ndarray
                Lower bounds of the function
            ub: np.ndarray
                Upper bounds of the function
            negate: bool
                negate the values before returning
        """

        super().__init__(lb,ub,negate)

    def evaluate(self, x):

        self.check_input(x)

        w = 1.0 + (x - 1.0) / 4.0

        part1 = np.sin(np.pi*w[:,0])**2
        
        part2 = np.sum( (w[:,:-1] - 1.0)**2 * (1.0 + 10.0*np.sin(np.pi*w[:,:-1] + 1.0)**2), axis=-1 )

        part3 = (w[:,-1] - 1.0)**2 * (1.0 + np.sin(2*np.pi*w[:,-1])**2 )

        return part1 + part2 + part3


class Rastrigin(AnalyticalProblem):

    def __init__(self, lb, ub, negate=False):
        """
            Class for defining the rastrigin function

            Parameters
            ----------
            lb: np.ndarray
                Lower bounds of the function
            ub: np.ndarray
                Upper bounds of the function
            negate: bool
                negate the values before returning
        """

        super().__init__(lb,ub,negate)

    def evaluate(self, x):

        self.check_input(x)

        return 10.0*self.dim + np.sum(x**2 - 10.0*np.cos(2.0*np.pi*x), axis=-1)


class Hartmann(AnalyticalProblem):

    def __init__(self, lb, ub, negate=False):
        """
            Class for defining the hartmann function

            Parameters
            ----------
            lb: np.ndarray
                Lower bounds of the function
            ub: np.ndarray
                Upper bounds of the function
            negate: bool
                negate the values before returning
        """

        super().__init__(lb,ub,negate)

        assert self.dim in [3,4,6], "Hartmann function is only defined for 3, 4 or 6 dimensions"

        if self.dim == 3:

            self.A = np.array([[3.0, 10, 30], [0.1, 10, 35], [3.0, 10, 30], [0.1, 10, 35]])

            self.P = np.array([[3689, 1170, 2673],
                               [4699, 4387, 7470],
                               [1091, 8732, 5547],
                               [381, 5743, 8828.0]])

        elif self.dim == 4:

            self.A = np.array([[10, 3, 17, 3.5],
                               [0.05, 10, 17, 0.1],
                               [3, 3.5, 1.7, 10],
                               [17, 8, 0.05, 10],])
            
            self.P = np.array([[1312, 1696, 5569, 124],
                               [2329, 4135, 8307, 3736],
                               [2348, 1451, 3522, 2883],
                               [4047, 8828, 8732, 5743]])

        elif self.dim == 6:

            self.A = np.array([[10, 3, 17, 3.5, 1.7, 8],
                               [0.05, 10, 17, 0.1, 8, 14],
                               [3, 3.5, 1.7, 10, 17, 8],
                               [17, 8, 0.05, 10, 0.1, 14]])

            self.P = np.array([[1312, 1696, 5569, 124, 8283, 5886],
                               [2329, 4135, 8307, 3736, 1004, 9991],
                               [2348, 1451, 3522, 2883, 3047, 6650],
                               [4047, 8828, 8732, 5743, 1091, 381]])

        self.alpha = np.array([[1], [1.2], [3.0], [3.2]])

    def evaluate(self, x):

        self.check_input(x)

        innersum = np.sum(self.A * (np.expand_dims(x,-2) - 1e-4*self.P)**2, axis=-1)

        y = -np.matmul(np.exp(-innersum), self.alpha)

        if self.dim == 4:
            y = (1.1 + y)/0.839

        return y
