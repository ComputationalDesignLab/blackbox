############## Script file for running airfoil analysis.
# Imports
import pickle, os
from mpi4py import MPI
from adflow import ADFLOW

# Getting MPI comm
comm = MPI.COMM_WORLD
parent_comm = comm.Get_parent()

# Redirecting the stdout - only root processor does printing
if comm.rank == 0:
    stdout = os.dup(1)
    log = open("log.txt", "a")
    os.dup2(log.fileno(), 1)

# Send the processor
parent_comm.send(os.getpid(), dest=0, tag=comm.rank)

try:
    ############## Reading input file for the analysis

    # Reading input file
    filehandler = open("input.pickle", 'rb') 
    input = pickle.load(filehandler)
    filehandler.close()

    # Getting aero problem from input file
    ap = input["aeroProblem"]
    slice = input["sliceLocation"]
    writeLiftDistribution = input["writeLiftDistribution"]
    CL_target = input["targetCL"]
    tol = input["targetCLTol"]
    alpha0 = input["startingAlpha"]

    # Assigning non-shape DVs
    if "mach" in input.keys():
        ap.mach = input["mach"][0]

    if "altitude" in input.keys():
        ap.altitude = input["altitude"][0]

    # Getting solver and meshing options from input file
    solverOptions = input["solverOptions"]
    solverOptions["gridFile"] = "volMesh.cgns"
    
    ############## Settign up adflow

    if comm.rank == 0:
        print("")
        print("#" + "-"*129 + "#")
        print(" "*59 + "Analysis Log" + ""*59)
        print("#" + "-"*129 + "#")
        print("")

    # Creating adflow object
    CFDSolver = ADFLOW(options=solverOptions, comm=comm)

    # Adding pressure distribution output
    if solverOptions["liftIndex"] == 2: # y
        sliceDirection = "z"
    elif solverOptions["liftIndex"] == 3: # z
        sliceDirection = "y"

    for loc in slice: 
        CFDSolver.addSlices(sliceDirection, loc, sliceType="absolute")

    if writeLiftDistribution:
        CFDSolver.addLiftDistribution(nSegments=200, direction=sliceDirection)

    ############## Solving for the CL
    solveCLResults = CFDSolver.solveCL(ap, CLStar=CL_target, alpha0=alpha0, tol=tol, autoReset=False, maxIter=10, writeSolution=True)

    # Evaluating objectives
    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    CFDSolver.checkSolutionFailure(ap, funcs)

    # printing the result
    if comm.rank == 0:
        print("")
        print("#" + "-"*129 + "#")
        print(" "*59 + "Result" + ""*59)
        print("#" + "-"*129 + "#")
        print("")

        output = {}

        # Printing and storing results based on evalFuncs in aero problem
        for obj in ap.evalFuncs:
            print("{} = ".format(obj), funcs["{}_{}".format(ap.name, obj)])
            output["{}".format(obj)] = funcs["{}_{}".format(ap.name, obj)]

        # Other mandatory outputs
        print("fail = ", not solveCLResults["converged"])
        output["fail"] = not solveCLResults["converged"]

        # Storing the results in output file
        filehandler = open("output.pickle", "xb")
        pickle.dump(output, filehandler)
        filehandler.close()

except Exception as e:
    if comm.rank == 0:
        print(e)

finally:
    # Redirecting to original stdout
    if comm.rank == 0:
        os.dup2(stdout, 1)

        # close the file and stdout
        log.close()
        os.close(stdout)

    # Getting intercomm and disconnecting
    # Otherwise, program will enter deadlock
    parent_comm.Disconnect()
