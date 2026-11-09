# rst EXAMPLE start
from pyDOE3 import lhs
from scipy.io import savemat
from blackbox.analytical_problems import Ackley

dim = 10

# initialize problem
problem = Ackley(num_inputs=dim)

# samples
x = lhs(dim, 10, criterion="cm")
x = problem.lb + (problem.ub - problem.lb) * x

y = problem(x) # evaluate samples

# save data
savemat('ackley_values.mat', {'x': x, 'y': y})
# rst EXAMPLE end