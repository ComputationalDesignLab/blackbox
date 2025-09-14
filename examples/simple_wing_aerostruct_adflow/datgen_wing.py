from blackbox import AeroStructFFD
from baseclasses import AeroProblem
import numpy as np
from tacs import constitutive, elements, functions

solverOptions = {
    # Common Parameters
    "monitorvariables": ["cl", "cd", "yplus"],
    "writeTecplotSurfaceSolution": True,
    "writeSurfaceSolution": False,
    "writeVolumeSolution": False,
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

# Callback function used to setup TACS element objects and DVs
def element_callback(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):

    # Material properties
    rho = 2500.0  # density kg/m^3
    E = 70.0e9  # Young's modulus (Pa)
    nu = 0.30  # Poisson's ratio
    ys = 350e6  # yield stress
    t = 0.01  # shell thickness, m

    # Setup (isotropic) property and constitutive objects
    prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
    # Set one thickness dv for every component
    con = constitutive.IsoShellConstitutive(prop, t=t, tNum=dvNum)

    # For each element type in this component, pass back the appropriate tacs element object
    transform = None
    elem = elements.Quad4Shell(transform, con)

    return elem

def problem_setup(scenario_name, fea_assembler, problem):
    """
    Helper function to add fixed forces and eval functions
    to structural problems used in tacs builder
    """

    # Add TACS Functions
    if scenario_name == "cruise":
        problem.addFunction("mass", functions.StructuralMass)
    problem.addFunction(
        "ks_vmfailure", functions.KSFailure, safetyFactor=1.0, ksWeight=100.0
    )

    # Add gravity load
    g = np.array([0.0, 0.0, -9.81])  # m/s^2
    if scenario_name == "maneuver":
        problem.addInertialLoad(2.5 * g)
    else:  # cruise
        problem.addInertialLoad(g)

# Flow details
ap = AeroProblem(name="wing", alpha=2.5, mach=0.85, altitude=10000, areaRef=45.5, chordRef=3.56, evalFuncs=["cl", "cd"])

options = {
    "aeroSolver": "adflow",
    "solverOptions": solverOptions,
    "tacsElementCallback": element_callback,
    "tacsProblemSetup": problem_setup,
    "gridFile": "wing_volMesh.cgns",
    "ffdFile": "wing_ffd.xyz",
    "structMeshFile": "wingbox.bdf",
    "liftIndex": 2, # Very important
    "aeroProblem": ap,
    "noOfProcessors": 8,
    "sliceLocation": [0.14, 3.22, 6.3, 9.38, 12.46, 13.86],
    "writeDeformedFFD": True,
    "samplingCriterion": "ese"
}

# Create the wing object
wing = AeroStructFFD(options=options)

# Add alpha as a design variable
wing.addDV("alpha", lowerBound=1.5, upperBound=4.5)

# Add the wing shape as a design variable
lowerBound = np.array([-0.03]*wing.nffd)
upperBound = np.array([0.03]*wing.nffd)
wing.addDV("shape", lowerBound=lowerBound, upperBound=upperBound)

# Add the wing twist as a design variable
lowerBound = np.array([-2.0]*wing.nTwist)
upperBound = np.array([2.0]*wing.nTwist)
wing.addDV("twist", lowerBound=lowerBound, upperBound=upperBound)

