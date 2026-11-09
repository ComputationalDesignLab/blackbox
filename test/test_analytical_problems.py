import unittest
from scipy.stats.qmc import LatinHypercube
from blackbox.analytical_problems import Ackley, Levy, Rastrigin, Hartmann, Branin, ModifiedBranin
import numpy as np

class SingleOutputAnalyticalProblemTestCaseMixin:

    def test_output_shape(self):

        for function in self.functions:

            # test function values
            sampler = LatinHypercube(d=function.num_inputs)
            x = sampler.random(self.num_samples)
            x = function.bounds[0] + (function.bounds[1] - function.bounds[0]) * x
            y = function(x)

            self.assertTrue(
                y.shape == (self.num_samples,1),
                "output shape for multi sample input is not correct"
            )

            x = np.random.rand(function.num_inputs)
            x = function.bounds[0] + (function.bounds[1] - function.bounds[0]) * x
            y = function(x)

            self.assertTrue(
                y.shape == (1,),
                "output shape for a single sample input is not correct"
            )

class AnalyticalOptimizationProblemTestCaseMixin:
    """
        Mixin class for adding test cases for analytical optimization
        problems that have known optimum value and location
    """

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
            np.testing.assert_allclose(y_opt_comptued, y_opt, atol=1e-6)

            # test function values
            sampler = LatinHypercube(d=function.num_inputs)
            x = sampler.random(self.num_samples)
            x = function.bounds[0] + (function.bounds[1] - function.bounds[0]) * x
            y = function(x)
            self.assertTrue(
                np.all(y <= y_opt) if function.negate else np.all(y >= y_opt),
                "function implementation is not correct, some of the y values are better than function extremum"
            )

class TestAckley(
    unittest.TestCase,
    SingleOutputAnalyticalProblemTestCaseMixin,
    AnalyticalOptimizationProblemTestCaseMixin
):
    
    num_samples = 10
    functions = [Ackley(), Ackley(num_inputs=5),  Ackley(negate=True), Ackley(num_inputs=15, negate=True)]

class TestLevy(
    unittest.TestCase,
    SingleOutputAnalyticalProblemTestCaseMixin,
    AnalyticalOptimizationProblemTestCaseMixin
):
    
    num_samples = 10
    functions = [Levy(), Levy(num_inputs=5),  Levy(negate=True), Levy(num_inputs=15, negate=True)]

class TestRastrigin(
    unittest.TestCase,
    SingleOutputAnalyticalProblemTestCaseMixin,
    AnalyticalOptimizationProblemTestCaseMixin
):
    
    num_samples = 10
    functions = [Rastrigin(), Rastrigin(num_inputs=5),  Rastrigin(negate=True), Rastrigin(num_inputs=15, negate=True)]

class TestHartmann(
    unittest.TestCase,
    SingleOutputAnalyticalProblemTestCaseMixin,
    AnalyticalOptimizationProblemTestCaseMixin
):
    
    num_samples = 10
    functions = [Hartmann(num_inputs=3), Hartmann(num_inputs=3, negate=True),
                 Hartmann(num_inputs=4), Hartmann(num_inputs=4, negate=True),
                 Hartmann(num_inputs=6), Hartmann(num_inputs=6, negate=True)]
    
class TestBranin(
    unittest.TestCase,
    SingleOutputAnalyticalProblemTestCaseMixin,
    AnalyticalOptimizationProblemTestCaseMixin
):
    
    num_samples = 10
    functions = [Branin(),  Branin(negate=True)]

class TestModifiedBranin(
    unittest.TestCase,
    SingleOutputAnalyticalProblemTestCaseMixin,
    AnalyticalOptimizationProblemTestCaseMixin
):
    
    num_samples = 10
    functions = [ModifiedBranin(),  ModifiedBranin(negate=True)]

if __name__ == '__main__':
    unittest.main()
