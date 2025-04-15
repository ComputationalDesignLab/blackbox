from blackbox import WingFFD
from baseclasses import AeroProblem
import numpy as np

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

ap = AeroProblem(name="wing", alpha=2.5, mach=0.85, altitude=10000, areaRef=45.5, chordRef=3.56, evalFuncs=["cl", "cd"])

options = {
    "solver": "adflow",
    "solverOptions": solverOptions,
    "gridFile": "wing_volMesh.cgns",
    "ffdFile": "wing_ffd.xyz",
    "liftIndex": 2, # Very important
    "aeroProblem": ap,
    "noOfProcessors": 8,
    "sliceLocation": [0.14, 3.22, 6.3, 9.38, 12.46, 13.86],
    "writeDeformedFFD": True,
    # "alpha": "implicit",
    # "targetCLTol": 1e-4,
    # "startingAlpha": 3.0,
    "samplingCriterion": "ese"
}

# Create the wing object
wing = WingFFD(options=options)

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

# Generate samples
wing.generateSamples(5)
