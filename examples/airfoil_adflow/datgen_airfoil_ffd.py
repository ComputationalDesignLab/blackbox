from blackbox import AirfoilFFD
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

nffd = 20 # Number of FFD points

# Options for blackbox
options = {
    # Requried options
    "airfoilFile": "rae2822.dat",
    "nffd": 20,
    # FFD Box options
    "fitted": True,
    "ymarginl": 0.015,
    "ymarginu": 0.015,
    # Sampling options
    "samplingCriterion": "ese",
    # Fixing the LE/TE
    "fixLETE": True,
    # Smoothing options
    "smoothing": True,
    "smoothingTolerance": 5e-4,
    "smoothingTheta": 0.6,
    # Other options
    "noOfProcessors": 8,
    "aeroProblem": ap,
    "solverOptions": solverOptions,
    "meshingOptions": meshingOptions,
    "writeAirfoilCoordinates": True,
    "plotAirfoil": True,
    "writeDeformedFFD": True,
}

# Example for generating samples
airfoil = AirfoilFFD(options=options)

# Add alpha as an design variables
airfoil.addDV("alpha", lowerBound=1.5, upperBound=3.5)

############ Setting bounds for shape design variables ############

coeff = airfoil.DVGeo.origFFDCoef
index = airfoil.DVGeo.getLocalIndex(0)

# Initialize bounds
lowerBound = np.zeros(int(coeff.shape[0]/2))
upperBound = np.zeros(int(coeff.shape[0]/2))

# Lift direction
liftIndex = 1 # y

# Computing the y distance between upper and lower surface FFD points
# 0 at thrid index denotes section location - there are two sections
# (0,1) at second index denotes lower and upper surface of FFD points
dist = coeff[index[:,1,0], liftIndex] - coeff[index[:,0,0], liftIndex]

# Lower and upper surface FFD point index in DV
lowerIndex = np.linspace(0, nffd-2, dist.shape[0], dtype=int) # lower surface ffd point index in DV
upperIndex = np.linspace(1, nffd-1, dist.shape[0], dtype=int) # upper surface ffd point index in DV

# FFD local thickness
delta = 0.15

# Lower bounds
lowerBound[lowerIndex] = -delta * dist
lowerBound[upperIndex] = -delta * dist

# Upper bounds
upperBound[lowerIndex] = delta * dist
upperBound[upperIndex] = delta * dist

# Add shape as a design variables
# First and last entry is dropped since option for fixing LE/TE is set in blackbox
airfoil.addDV("shape", lowerBound=lowerBound[1:-1], upperBound=upperBound[1:-1])

# Generate samples
airfoil.generateSamples(5)
