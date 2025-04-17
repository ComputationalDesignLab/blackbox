from blackbox import AirfoilCST
from baseclasses import AeroProblem
import numpy as np

solverOptions = {
    # Common Parameters
    "monitorvariables": ["cl", "cd", "cmz", "yplus"],
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

meshingOptions = {
    # ---------------------------
    #        Input Parameters
    # ---------------------------
    "unattachedEdgesAreSymmetry": False,
    "outerFaceBC": "farfield",
    "autoConnect": True,
    "BC": {1: {"jLow": "zSymm", "jHigh": "zSymm"}},
    "families": "wall",
    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    "N": 129,
    "s0": 1e-6,
    "marchDist": 100.0,
}

# Creating aeroproblem for adflow
ap = AeroProblem(
    name="ap", alpha=2.0, mach=0.734, reynolds=6.5e6, reynoldsLength=1.0, T=288.15, 
    areaRef=1.0, chordRef=1.0, evalFuncs=["cl", "cd", "cmz"], xRef = 0.25, yRef = 0.0, zRef = 0.0
)

# Options for blackbox
options = {
    "solverOptions": solverOptions,
    "noOfProcessors": 8,
    "aeroProblem": ap,
    "airfoilFile": "rae2822.dat",
    "numCST": [6, 6],
    "meshingOptions": meshingOptions,
    "writeAirfoilCoordinates": True,
    "plotAirfoil": True,
    "writeSliceFile": True,
    "samplingCriterion": "ese"
}

# Example for generating samples
airfoil = AirfoilCST(options=options)

######### Multi Analysis

# Adding design variable
airfoil.addDV("alpha", 2.0, 3.0)

# Adding lower surface CST coeffs as DV
coeff = airfoil.DVGeo.lowerCST # get the fitted CST coeff
lb = coeff - np.sign(coeff)*0.3*coeff
ub = coeff + np.sign(coeff)*0.3*coeff
airfoil.addDV("lower", lowerBound=lb, upperBound=ub)

# Adding upper surface CST coeffs as DV
coeff = airfoil.DVGeo.upperCST # get the fitted CST coeff
lb = coeff - np.sign(coeff)*0.3*coeff
ub = coeff + np.sign(coeff)*0.3*coeff
airfoil.addDV("upper", lowerBound=lb, upperBound=ub)

# Generating the samples
airfoil.generateSamples(numSamples=5)
