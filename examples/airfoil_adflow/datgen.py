from blackbox.airfoil_adflow_problem import AirfoilADflow, AirfoilADflowOptions
from baseclasses import AeroProblem
import numpy as np
from pyDOE3 import lhs

solver_options = {
    # Common Parameters
    "monitorvariables": ["cl", "cd", "yplus"],
    # Physics Parameters
    "equationType": "RANS",
    "smoother": "DADI",
    "MGCycle": "sg",
    "nsubiterturb": 10,
    "nCycles": 7000,
    # ANK Solver Parameters
    "useANKSolver": True,
    "ANKSubspaceSize": 400,
    "ANKASMOverlap": 3,
    "ANKPCILUFill": 4,
    "ANKJacobianLag": 5,
    "ANKOuterPreconIts": 3,
    "ANKInnerPreconIts": 3,
    # NK Solver Parameters
    "useNKSolver": True,
    "NKSwitchTol": 1e-6,
    "NKSubspaceSize": 400,
    "NKASMOverlap": 3,
    "NKPCILUFill": 4,
    "NKJacobianLag": 5,
    "NKOuterPreconIts": 3,
    "NKInnerPreconIts": 3,
    # Termination Criteria
    "L2Convergence": 1e-14
}

# Creating aeroproblem for adflow
aero_problem = AeroProblem(
    name="airfoil", alpha=2.0, mach=0.734, reynolds=6.5e6, reynoldsLength=1.0, T=288.15, 
    areaRef=1.0, chordRef=1.0, xRef = 0.25, yRef = 0.0, zRef = 0.0
)

options = AirfoilADflowOptions(
    airfoil_file="airfoil.dat",
    solver_options=solver_options,
    aero_problem = aero_problem,
    plot_airfoil = True,
    alpha = "explicit",
    write_surface_output = True,
    write_volume_output = True
)

# Example for generating samples
airfoil = AirfoilADflow(options=options)

# Add lower surface as a parameter
coeff = airfoil.parametrization.lower_cst # get the fitted CST coeff
lb = coeff - np.sign(coeff)*0.3*coeff
ub = coeff + np.sign(coeff)*0.3*coeff

airfoil.add_parameter("lower_cst", lb, ub)

# Add upper surface as a parameter
coeff = airfoil.parametrization.upper_cst # get the fitted CST coeff
lb = coeff - np.sign(coeff)*0.3*coeff
ub = coeff + np.sign(coeff)*0.3*coeff

airfoil.add_parameter("upper_cst", lb, ub)

# Add mach number as a parameter
airfoil.add_parameter("mach", 0.65, 0.7)

# Add alpha as a parameter
airfoil.add_parameter("alpha", 2.0, 2.5)

# Add reynolds number as a parameter
airfoil.add_parameter("reynolds", 6e6, 6.5e6)

# generate sampling for evaluation
samples = lhs(airfoil.bounds[0].shape[0], 3, "cm", 100)
samples = airfoil.bounds[0] + (airfoil.bounds[1] - airfoil.bounds[0]) * samples

# evaluate samples
airfoil(samples, return_results=False)
