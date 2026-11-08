# rst start
from pyDOE3 import lhs
from blackbox.rover_trajectory_problem import Rover, plot_traj

dim = 80
force_start = True
force_goal = True

# create rover object
problem = Rover(dim=dim, force_start=force_start, force_goal=force_goal)

# create sampling plan
x = lhs(problem.dim, 10, criterion="cm", random_state=10)
x = problem.lb + (problem.ub - problem.lb) * x 

# evaluate sampling plan
y = problem(x)

print(f"Problem dimension: {problem.dim}")
print(f"\nSampling plan:\n{x}")
print(f"\nOutput:\n{y}")

# plot trajectory for first sample
plot_traj(control_points=x[0], rover=problem, show_figure=False)
# rst end