from blackbox.analytical_problems import Levy
from blackbox.doe import generate_lhs_samples
import numpy as np
from pyDOE3 import lhs

dim = 25

problem = Levy(dim, negate=False)

x = lhs(dim, 10, criterion="cm")

x = problem.lb + (problem.ub - problem.lb) * x

print(x)

y = problem(x)

print(y)

print(x[0,:])

print(problem(x[0,:]))
