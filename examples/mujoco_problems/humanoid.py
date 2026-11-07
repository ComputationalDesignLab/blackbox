from pyDOE3 import lhs
from blackbox.mujoco_problems import MujocoProblem

# initialize mujoco problem class
problem = MujocoProblem(env_name='Humanoid-v5', num_rollouts=5, seed=1)

# create sampling plan
x = lhs(problem.dim, 10, criterion="cm", random_state=1) 
x = problem.lb + (problem.ub - problem.lb) * x 

# evaluate sampling plan
y = problem(x)

print(f"Problem dimension: {problem.dim}")
print(f"\nSampling plan:\n{x}")
print(f"\nOutput:\n{y}")
