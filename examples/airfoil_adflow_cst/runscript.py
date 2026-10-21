# rst INIT start
# options for adflow
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
    # Termination Criteria
    "L2Convergence": 1e-14
}

from baseclasses import AeroProblem

# creating aeroproblem for defining flow conditions
aero_problem = AeroProblem(
    name="airfoil", alpha=2.0, mach=0.734, reynolds=6.5e6, reynoldsLength=1.0, T=288.15, 
    areaRef=1.0, chordRef=1.0, xRef = 0.25, yRef = 0.0, zRef = 0.0
)

from blackbox.airfoil_adflow_cst_problem import AirfoilADflowCSTOptions

options = AirfoilADflowCSTOptions(
    airfoil_file="airfoil.dat",
    solver_options=solver_options,
    aero_problem=aero_problem,
    plot_airfoil=True,
    alpha="explicit",
    scalar_outputs=['cl', 'cd', 'cm'],
    surface_outputs=['cp', 'cf'],
    volume_outputs=['cp', 'mach'],
    write_surface_output=True,
    write_volume_output=True,
    num_processors=4
)

from blackbox.airfoil_adflow_cst_problem import AirfoilADflowCST

airfoil = AirfoilADflowCST(options=options)
# rst INIT end

# rst PARM start
import numpy as np

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
# rst PARM end

# rst RUN start
from pyDOE3 import lhs

# generate sampling plan for evaluation
samples = lhs(airfoil.bounds[0].shape[0], 10, "cm", 100)
samples = airfoil.bounds[0] + (airfoil.bounds[1] - airfoil.bounds[0]) * samples

# evaluate samples
airfoil(x=samples, return_results=False)
# rst RUN end