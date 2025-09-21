from blackbox import AeroStructFFD
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

# Flow details
ap = AeroProblem(name="aerostruct", alpha=2.5, mach=0.85, altitude=10000, areaRef=45.5, chordRef=3.56, evalFuncs=["cl", "cd", "cmz"])

options = {
    "aeroSolver": "adflow",
    "aeroSolverOptions": solverOptions,
    "structSolverConfigFile": "tacs_setup.py",
    "gridFile": "wing_volMesh.cgns",
    "ffdFile": "wing_ffd.xyz",
    "structMeshFile": "wingbox.bdf",
    "liftIndex": 2, # Very important
    "aeroProblem": ap,
    "noOfProcessors": 4,
    "sliceLocation": [0.14, 3.22, 6.3, 9.38, 12.46, 13.86],
    "writeDeformedFFD": True,
    "writeForceField": True,
    "writeDisplacementField": True,
    "samplingCriterion": "ese"
}

# Create the wing object
wing = AeroStructFFD(options=options)

# Add alpha as a design variable
wing.addDV("mach", lowerBound=0.7, upperBound=0.9)

# Add alpha as a design variable
wing.addDV("alpha", lowerBound=1.5, upperBound=4.5)

# Add the wing shape as a design variable
lowerBound = np.array([-0.1]*wing.nffd)
upperBound = np.array([0.1]*wing.nffd)
wing.addDV("shape", lowerBound=lowerBound, upperBound=upperBound)

# Add the wing twist as a design variable
lowerBound = np.array([-5.0]*wing.nTwist)
upperBound = np.array([5.0]*wing.nTwist)
wing.addDV("twist", lowerBound=lowerBound, upperBound=upperBound)

wing.generateSamples(5)
