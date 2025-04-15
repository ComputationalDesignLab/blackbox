# Runscript to get forces, only execute this script after running the analysis

import numpy as np
from mpi4py import MPI
from dafoam import PYDAFOAM
import pickle
from scipy.io import savemat

# Reading input file
filehandler = open("input.pickle", 'rb') 
input = pickle.load(filehandler)
filehandler.close()

# Getting aero problem from input file
ap = input["aeroProblem"]

# Assigning non-shape DVs to get correct pref
if "alpha" in input.keys():
    ap.alpha = input["alpha"][0]

if "mach" in input.keys():
    ap.mach = input["mach"][0]

if "altitude" in input.keys():
    ap.altitude = input["altitude"][0]

daOptions = {
    "designSurfaces": ["wing"],
    "solverName": "DARhoSimpleCFoam",
    "outputInfo": {
        "f_aero": {
            "type": "forceCouplingOutput",
            "patches": ["wing"],
            "components": ["forceCoupling"],
            "pRef": float(ap.P),
        },
    },
    "printDAOptions": False
}

DASolver = PYDAFOAM(options=daOptions, comm=MPI.COMM_WORLD)

outputSize = DASolver.solver.getOutputSize("f_aero", "forceCouplingOutput")

forces = np.zeros(outputSize)

DASolver.solver.calcOutput("f_aero", "forceCouplingOutput", forces)

forces = forces.reshape((-1, 3))

surfaceCoords = DASolver.getSurfaceCoordinates()

data = {
    "forces": forces,
    "surfaceCoords": surfaceCoords
}

# Save the data file
savemat("surfaceForces.mat", data)
