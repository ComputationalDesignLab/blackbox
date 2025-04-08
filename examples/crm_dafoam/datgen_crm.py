from blackbox import WingFFD
from baseclasses import AeroProblem
import numpy as np

# options for DAFoam
solverOptions = {
    "gridFile": "crm_volMesh.cgns",
    "liftIndex": 3,
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

# Creating aeroproblem for adflow
ap = AeroProblem(
    name="crm", alpha=2.2, mach=0.85, reynolds=5e6, reynoldsLength=1.0, T=298.15,
    areaRef=3.407014, chordRef=1.0, evalFuncs=["cl", "cd", "cmy"], xRef = 1.2077, yRef = 0.0, 
    zRef = 0.007669
)

options = {
    "solver": "dafoam",
    "solverOptions": solverOptions,
    "ffdFile": "crm_ffd.xyz",
    "aeroProblem": ap,
    "noOfProcessors": 8,
    "sliceLocation": [0.883, 1.003, 2.093, 2.612, 3.112, 3.548],
    "writeDeformedFFD": True
}

# Create the wing object
wing = WingFFD(options=options)

# Add alpha as a design variable
wing.addDV("alpha", lowerBound=1.5, upperBound=3.5)

# Add mach as a design variable
wing.addDV("mach", lowerBound=0.7, upperBound=0.8)

# Add the wing shape as a design variable
# lowerBound = np.array([-0.01]*wing.nffd)
# upperBound = np.array([0.01]*wing.nffd)
# wing.addDV("shape", lowerBound=lowerBound, upperBound=upperBound)

# Add the wing twist as a design variable
# lowerBound = np.array([-2.0]*wing.nTwist)
# upperBound = np.array([2.0]*wing.nTwist)
# wing.addDV("twist", lowerBound=lowerBound, upperBound=upperBound)

# Generate samples
wing.generateSamples(2)
