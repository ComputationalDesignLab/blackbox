# Runscript to get forces, only execute this script after running the analysis

import numpy as np
from mpi4py import MPI
from dafoam import PYDAFOAM
import pickle, os
from scipy.io import savemat

# Getting MPI comm
comm = MPI.COMM_WORLD
parent_comm = comm.Get_parent()

# Redirecting the stdout - only root processor does printing
if comm.rank == 0:
    stdout = os.dup(1)
    log = open("log_forces.txt", "a")
    os.dup2(log.fileno(), 1)

# Send the processor
parent_comm.send(os.getpid(), dest=0, tag=comm.rank)

try: 

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

except Exception as e:
    if comm.rank == 0:
        print(e)

finally:
    # close the file
    if comm.rank == 0:
        # Redirecting to original stdout
        os.dup2(stdout, 1)
        os.close(stdout)
        log.close()

    # Important to wait for all processors to finish before calling disconnect
    # Otherwise, program will enter deadlock
    comm.barrier()
        
    # Getting intercomm and disconnecting
    # Otherwise, program will enter deadlock
    parent_comm.Disconnect()
