from blackbox.mujoco_problems import HalfCheetah
from pyDOE3 import lhs

problem = HalfCheetah()

x = lhs(problem.dim, 10, criterion="cm")

x = problem.lb + (problem.ub - problem.lb) * x

print(x)

y = problem(x)

print(y)