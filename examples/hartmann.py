from blackbox.analytical_problems import Hartmann
from blackbox.doe import generate_lhs_samples
import numpy as np

dim = 3
lb = np.zeros(dim)
ub = np.ones(dim)

problem = Hartmann(lb, ub, negate=True)

x = generate_lhs_samples(problem, 20)

print(x)

y = problem(x)

print(y)

print(problem(np.array([0.114614, 0.555649, 0.852547])))

dim = 6
lb = np.zeros(dim)
ub = np.ones(dim)

problem = Hartmann(lb, ub, negate=True)

x = generate_lhs_samples(problem, 20)

print(x)

y = problem(x)

print(y)

print(problem(np.array([0.20169, 0.150011, 0.476874, 0.275332, 0.311652, 0.6573])))

dim = 4
lb = np.zeros(dim)
ub = np.ones(dim)

problem = Hartmann(lb, ub, negate=True)

x = generate_lhs_samples(problem, 20)

print(x)

y = problem(x)

print(y)
