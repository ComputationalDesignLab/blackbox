from blackbox.analytical_problems import Ackley
from blackbox.doe import generate_lhs_samples
import numpy as np
from scipy.io import savemat
from pyDOE3 import lhs

dim = 10

problem = Ackley(num_inputs=dim)

x = lhs(dim, 10, criterion="cm")

x = problem.lb + (problem.ub - problem.lb) * x

print(x)

y = problem(x)

print(y)

print(x[0,:])

print(problem(x[0,:]))

# data = {
#     'x': x,
#     'y': y
# }

# savemat('ackley_values.mat', data)
