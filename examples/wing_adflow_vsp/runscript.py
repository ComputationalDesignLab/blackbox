# rst FLOW start
from baseclasses import AeroProblem

aero_problem = AeroProblem(name="wing", alpha=2.5, mach=0.85, altitude=10000, areaRef=45.5, chordRef=3.56)

aero_problem.addDV(key="mach", lower=0.75, upper=0.8)

aero_problem.addDV(key="alpha", lower=1, upper=4)

aero_problem.addDV(key="altitude", lower=8000, upper=12000)
# rst FLOW end

# rst SHAPE start
import numpy as np
from blackbox.wing_problem import WingVSP

wing_vsp = WingVSP("wing.vsp3")

# add wing planform variables
wing_vsp.add_shape_parameter(name="Dihedral", lower=0.0, upper=5.0, sections=[1])
wing_vsp.add_shape_parameter(name="Sweep", lower=25.0, upper=35.0, sections=[1])
wing_vsp.add_shape_parameter(name="Span", lower=12.0, upper=16.0, sections=[1])
wing_vsp.add_shape_parameter(name="Twist", lower=-5.0, upper=0.0, sections=[1])

# cst coefficients
upper_cst_coeff = np.array([0.12344046, 0.14809657, 0.14858145, 0.2168004,  0.17607825, 0.21018404])
lower_cst_coeff = np.array([-0.13198943, -0.11895939, -0.22056435, -0.12743513, -0.08232715, 0.05055414])

# add airfoil upper surface cst coefficients
wing_vsp.add_airfoil_cst_parameters(
    "upper", 
    upper_cst_coeff - np.sign(upper_cst_coeff)*0.3*upper_cst_coeff,
    upper_cst_coeff + np.sign(upper_cst_coeff)*0.3*upper_cst_coeff,
    sections=[0,1]
)

# add airfoil lower surface cst coefficients
wing_vsp.add_airfoil_cst_parameters(
    "lower", 
    lower_cst_coeff - np.sign(lower_cst_coeff)*0.3*lower_cst_coeff,
    lower_cst_coeff + np.sign(lower_cst_coeff)*0.3*lower_cst_coeff,
    sections=[0,1]
)
# rst SHAPE end

# rst INIT start
from blackbox.wing_problem import WingADflowVSPOptions, WingADflowVSP

solver_options = {
    # Common Parameters
    "monitorvariables": ["cl", "cd", "yplus"],
    # Physics Parameters
    "equationType": "RANS",
    "smoother": "DADI",
    "MGCycle": "sg",
    "nsubiterturb": 10,
    "nCycles": 7000,
    "liftIndex": 3,
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

options = WingADflowVSPOptions(
    wing_vsp=wing_vsp,
    mesh_file="wing_vol_mesh.cgns",
    solver_options=solver_options,
    aero_problem=aero_problem,
    write_surface_output=True,
    write_volume_output=True,
    num_processors=10,
    write_vsp_file=True,
    write_stl_file=True,
    write_slice_file=True,
    write_lift_distribution=True
)

wing = WingADflowVSP(options=options)
# rst INIT end

# rst RUN start
from pyDOE3 import lhs

# generate sampling plan for evaluation
samples = lhs(n=wing.bounds[0].shape[0], samples=10, criterion="cm", iterations=100)
samples = wing.bounds[0] + (wing.bounds[1] - wing.bounds[0]) * samples

wing(samples, return_results=False)
# rst RUN end