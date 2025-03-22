# Script for generating airfoil samples 

from blackbox import AirfoilCST
from baseclasses import AeroProblem
import numpy as np

# Input Parameters
U0 = 238.0
p0 = 101325.0
T0 = 300.0
nuTilda0 = 4.5e-5
CL_target = 0.5
aoa0 = 3.0
A0 = 1.0
# rho is used for normalizing CD and CL
rho0 = p0 / T0 / 287

# Shift most of the solver options inside
# function can be based on the evalfuncs from aeroproblem
# inputinfo can be based on operating condition design variables
# primalBC can be done inside
# constant variables U0, p0, T0, nutilda, rho0 can be done in runscript 

# options for DAFoam
solverOptions = {
    "designSurfaces": ["wing"], # internal
    "solverName": "DARhoSimpleCFoam", # user
    "primalMinResTol": 1e-8, # user
    "primalBC": { # internal
        "U0": {"variable": "U", "patches": ["inout"], "value": [U0, 0.0, 0.0]},
        "p0": {"variable": "p", "patches": ["inout"], "value": [p0]},
        "T0": {"variable": "T", "patches": ["inout"], "value": [T0]},
        "nuTilda0": {"variable": "nuTilda", "patches": ["inout"], "value": [nuTilda0]},
        "useWallFunction": True,
    },
    "function": { # internal
        "CD": {
            "type": "force",
            "source": "patchToFace",
            "patches": ["wing"],
            "directionMode": "parallelToFlow",
            "patchVelocityInputName": "patchV",
            "scale": 1.0 / (0.5 * U0 * U0 * A0 * rho0),
        },
        "CL": {
            "type": "force",
            "source": "patchToFace",
            "patches": ["wing"],
            "directionMode": "normalToFlow",
            "patchVelocityInputName": "patchV",
            "scale": 1.0 / (0.5 * U0 * U0 * A0 * rho0),
        },
    },
    "normalizeStates": { # internal
        "U": U0,
        "p": p0,
        "T": T0,
        "nuTilda": nuTilda0 * 10.0,
        "phi": 1.0,
    },
    "inputInfo": { # internal
        "patchV": {
            "type": "patchVelocity",
            "patches": ["inout"],
            "flowAxis": "x",
            "normalAxis": "y",
            "components": ["solver", "function"],
        },
    },
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
airfoil.generateSamples(numSamples=2)
