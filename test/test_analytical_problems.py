import unittest, os
from scipy.io import loadmat
from blackbox.analytical_problems import Ackley
import numpy as np

baseDir = os.path.dirname(os.path.abspath(__file__))

class TestAckley(unittest.TestCase):

    def test_ackley1(self):
        
        data = loadmat(os.path.join(baseDir, 'ref/ackley_values.mat'))

        x = data['x']
        y = data['y']

        dim = x.shape[1]
        lb = -32.768*np.ones(dim)
        ub = 32.768*np.ones(dim)

        problem = Ackley(lb, ub)

        np.testing.assert_allclose(problem(x), y, rtol=1e-6)

if __name__ == '__main__':
    unittest.main()