from blackbox.analytical_problems import Ackley
from blackbox.doe import generate_lhs_samples
import numpy as np

dim = 10
lb = -32.768 * np.ones(dim)
ub = 32.768 * np.ones(dim)

problem = Ackley(lb, ub)

x = generate_lhs_samples(problem, 10)

print(x)

y = problem(x)

print(y)