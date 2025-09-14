from blackbox.analytical_problems import Levy
from blackbox.doe import generate_lhs_samples
import numpy as np

dim = 25
lb = -10 * np.ones(dim)
ub = 10 * np.ones(dim)

problem = Levy(lb, ub, negate=True)

x = generate_lhs_samples(problem, 20)

print(x)

y = problem(x)

print(y)

print(x[0,:])

print(problem(x[0,:]))
