import pickle, os, time
from mpi4py import MPI
from baseclasses import AeroProblem
from adflow import ADFLOW
from idwarp import USMesh
from pyhyp import pyHyp
from cgnsutilities.cgnsutilities import readGrid

# Getting MPI comm
comm = MPI.COMM_WORLD
parent_comm = comm.Get_parent()

# Send the processor
parent_comm.send(os.getpid(), dest=0, tag=comm.rank)

# Redirecting the stdout
stdout = os.dup(1)
log = open("log.txt", "a")
os.dup2(log.fileno(), 1)

try:
    ############## Reading input file for the analysis

    # Reading input file
    filehandler = open("input.pickle", 'rb') 
    input = pickle.load(filehandler)
    filehandler.close()

    # Getting aero problem from input file
    ap = input["aeroProblem"]
    refine = input["refine"]
    slice = input["writeSliceFile"]
    CL_target = input["targetCL"]

    # Assigning non-shape DVs
    if "alpha" in input.keys():
        ap.alpha = input["alpha"][0]

    if "mach" in input.keys():
        ap.mach = input["mach"][0]

    if "altitude" in input.keys():
        ap.altitude = input["altitude"][0]

    # Getting solver and meshing options from input file
    solverOptions = input["solverOptions"]
    solverOptions["gridFile"] = "volMesh.cgns"
    solverOptions["liftindex"] = 2 # Always 2 since meshing is done internally

    meshingOptions = input["meshingOptions"]
    meshingOptions["inputFile"] = "surfMesh.xyz"

    ############## Generating mesh

    if comm.rank == 0:
        print("#" + "-"*129 + "#")
        print(" "*59 + "Meshing Log" + ""*59)
        print("#" + "-"*129 + "#")
        print("")

    hyp = pyHyp(options=meshingOptions, comm=comm)
    hyp.run()
    hyp.writeCGNS("volMesh.cgns")

    ############## Refining the mesh

    if comm.rank == 0:
        # Read the grid
        grid = readGrid("volMesh.cgns")

        # Only one processor has to do this
        if refine == 1:
            grid.refine(['i', 'k'])
        if refine == 2:
            grid.refine(['i', 'k'])
            grid.refine(['i', 'k'])
        if refine == -1:
            grid.coarsen()
        if refine == -2:
            grid.coarsen()
            grid.coarsen()

        grid.writeToCGNS("volMesh.cgns")
        
    # Wait till root is done with refining/coarse of mesh
    comm.barrier()

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
    if slice:
        CFDSolver.addSlices("z", 0.5, sliceType="absolute")

    ############## Solving for the CL
    CFDSolver.solveCL(ap, CLStar=CL_target, alpha0=2.8, delta=0.2, tol=0.0001, autoReset=False, maxIter=8)

    ############# Post-processing

    # Run CFD at obtained alpha - it will be quick
    CFDSolver(ap)

    # Evaluating objectives
    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    CFDSolver.checkSolutionFailure(ap, funcs)

    # printing the result
    if MPI.COMM_WORLD.rank == 0:
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
        print("fail = ", funcs["fail"])
        output["fail"] = funcs["fail"]

        # Storing the results in output file
        filehandler = open("output.pickle", "xb")
        pickle.dump(output, filehandler)
        filehandler.close()

except Exception as e:
    if comm.rank == 0:
        print(e)

finally:
    # Redirecting to original stdout
    os.dup2(stdout, 1)

    # close the file and stdout
    log.close()
    os.close(stdout)

    # Getting intercomm and disconnecting
    # Otherwise, program will enter deadlock
    parent_comm.Disconnect()