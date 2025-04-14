# Runscript to get forces, only execute this script after running the analysis

import numpy as np
from mpi4py import MPI
from dafoam import PYDAFOAM
import argparse
from scipy.io import savemat

parser = argparse.ArgumentParser(description="Get forces from DAFoam")
parser.add_argument("--pRef")
args = parser.parse_args()

daOptions = {
    "designSurfaces": ["wing"],
    "solverName": "DARhoSimpleCFoam",
    "outputInfo": {
        "f_aero": {
            "type": "forceCouplingOutput",
            "patches": ["wing"],
            "components": ["forceCoupling"],
            "pRef": float(args.pRef),
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
