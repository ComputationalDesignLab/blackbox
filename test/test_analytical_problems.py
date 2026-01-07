from abc import abstractmethod
import unittest
from scipy.stats.qmc import LatinHypercube
from blackbox.analytical_problems import Ackley
import numpy as np

class BaseSingleOutputAnalyticalProblemsTest(unittest.TestCase):

    # def setUp(self):
    #     """
    #         special method from unittest class - runs before each test function
    #     """

    #     self.assertTrue(
    #         hasattr(self, "functions"),
    #         "Test class must define `functions`"
    #     )

    #     self.assertTrue(
    #         hasattr(self, "num_samples"),
    #         "Test class must define `num_samples`"
    #     )

    def test_output_shape(self):

        for function in self.functions:

            # test function values
            sampler = LatinHypercube(d=function.num_inputs)
            x = sampler.random(self.num_samples)
            x = function.bounds[0] + (function.bounds[1] - function.bounds[0]) * x
            y = self.function(x)

            self.assertTrue(
                y.shape == (self.num_samples,1),
                "output shape for multi sample input is not correct"
            )

            x = np.random.rand(function.num_inputs)
            x = function.bounds[0] + (function.bounds[1] - function.bounds[0]) * x
            y = self.function(x)

            self.assertTrue(
                y.shape == (self.num_samples,1),
                "output shape for a single sample input is not correct"
            )

    @property
    @abstractmethod
    def functions(self) -> list:

        pass

class OptimizationFunctionTest(unittest.TestCase):

    def test_known_optimum(self):

        for function in self.functions:

            # check attributes
            self.assertTrue(
                hasattr(function, "_x_opt"),
                "`_x_opt` is not implemented in the class"
            )

            self.assertTrue(
                hasattr(function, "_y_opt"),
                "`_y_opt` is not implemented in the class"
            )

            # check function optimum
            y_opt = -function._y_opt if function.negate else function._y_opt
            y_opt_comptued = function(function._x_opt)
            np.testing.assert_allclose(y_opt_comptued, y_opt, atol=1e-8)

            # test function values
            sampler = LatinHypercube(d=function.num_inputs)
            x = sampler.random(self.num_samples)
            x = function.bounds[0] + (function.bounds[1] - function.bounds[0]) * x
            y = function(x)
            self.assertTrue(
                np.all(y >= y_opt),
                "function implementation is not correct, some of the y values are better than function extremum"
            )

    @property
    @abstractmethod
    def functions(self) -> list:

        pass

class TestAckley(
    BaseSingleOutputAnalyticalProblemsTest,
    OptimizationFunctionTest
):
    
    num_samples = 10
    functions = [Ackley(), Ackley(num_inputs=5),  Ackley(negate=True), Ackley(num_inputs=15, negate=True)]

if __name__ == '__main__':
    unittest.main()
