# Script for generating airfoil samples 
from blackbox import AirfoilCST
from baseclasses import AeroProblem
import numpy as np

# options for DAFoam
solverOptions = {
    "designSurfaces": ["wing"],
    "solverName": "DARhoSimpleCFoam",
    "primalMinResTol": 1e-6,
    "primalMinResTolDiff": 1e4,
    "checkMeshThreshold": {
            "maxAspectRatio": 1000.0,
            "maxNonOrth": 70.0,
            "maxSkewness": 4.0,
            "maxIncorrectlyOrientedFaces": 0,
        }
}

# options for pyhyp
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
    "N": 33,
    "s0": 4e-3,
    "marchDist": 20.0,
    # ---------------------------
    #   Pseudo Grid Parameters
    # ---------------------------
    "ps0": -1.0,
    "pGridRatio": -1.0,
    "cMax": 1.0,
    # ---------------------------
    #   Smoothing parameters
    # ---------------------------
    "epsE": 2.0,
    "epsI": 4.0,
    "theta": 2.0,
    "volCoef": 0.20,
    "volBlend": 0.0005,
    "volSmoothIter": 20,
}

# Creating aeroproblem
ap = AeroProblem(
    name="ap", alpha=2.0, mach=0.734, reynolds=6.5e6, reynoldsLength=1.0, T=288.15, 
    areaRef=1.0, chordRef=1.0, evalFuncs=["cl", "cd", "cmz"], xRef = 0.25, yRef = 0.0, zRef = 0.0
)

# Options for blackbox
options = {
    "solver": "dafoam",
    "solverOptions": solverOptions,
    "meshingOptions": meshingOptions,
    "noOfProcessors": 8,
    "aeroProblem": ap,
    "airfoilFile": "naca0012.dat",
    "numCST": [6, 6],
    "writeAirfoilCoordinates": True,
    "plotAirfoil": True,
    "writeSliceFile": True,
    "samplingCriterion": "ese"
}

# Example for generating samples
airfoil = AirfoilCST(options=options)

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
