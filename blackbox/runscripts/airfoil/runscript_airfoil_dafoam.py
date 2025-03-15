############## Script file for running airfoil analysis using DAfoam
# Imports
import pickle, os
from mpi4py import MPI
from pyhyp import pyHyp
from cgnsutilities.cgnsutilities import readGrid
import openmdao.api as om
from dafoam.mphys import DAFoamBuilder

# Getting MPI comm
comm = MPI.COMM_WORLD
parent_comm = comm.Get_parent()

# Send the processor
parent_comm.send(os.getpid(), dest=0, tag=comm.rank)

try:
    # Redirecting the stdout - only root processor does printing
    if comm.rank == 0:
        log = open("log.txt", "a")
        stdout = os.dup(1)
        os.dup2(log.fileno(), 1)

    ############## Reading input file for the analysis

    # Reading input file
    filehandler = open("input.pickle", 'rb') 
    input = pickle.load(filehandler)
    filehandler.close()

    # Getting aero problem from input file
    ap = input["aeroProblem"]
    refine = input["refine"]
    slice = input["writeSliceFile"]

    # Assigning non-shape DVs
    if "alpha" in input.keys():
        ap.alpha = input["alpha"][0]

    if "mach" in input.keys():
        ap.mach = input["mach"][0]

    if "altitude" in input.keys():
        ap.altitude = input["altitude"][0]

    # Getting solver and meshing options from input file
    solverOptions = input["solverOptions"]

    # solverOptions["gridFile"] = "volMesh.cgns"
    # solverOptions["liftindex"] = 2 # Always 2 since meshing is done internally

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
    # hyp.writeCGNS("volMesh.cgns")
    hyp.writePlot3D("volMesh.xyz")

    ############## Refining the mesh

    # Only one processor has to do this
    if comm.rank == 0:

    #     # Read the grid
    #     grid = readGrid("volMesh.cgns")

    #     if refine == 1:
    #         grid.refine(['i', 'k'])
    #     if refine == 2:
    #         grid.refine(['i', 'k'])
    #         grid.refine(['i', 'k'])
    #     if refine == -1:
    #         grid.coarsen()
    #     if refine == -2:
    #         grid.coarsen()
    #         grid.coarsen()

    #     grid.writePlot3d("volMesh.xyz")
    #     os.system("rm volMesh.cgns")

        # Converting plot3d vol mesh to openfoam friendly vol mesh
        os.system("plot3dToFoam -noBlank volumeMesh.xyz")
        os.system("autoPatch 30 -overwrite")
        os.system("createPatch -overwrite")
        os.system("renumberMesh -overwrite")

    # Wait till root is done with refining/coarse of mesh
    comm.barrier()
    
    ############## Settign up dafoam

    if comm.rank == 0:
        print("")
        print("#" + "-"*129 + "#")
        print(" "*59 + "Analysis Log" + ""*59)
        print("#" + "-"*129 + "#")
        print("")

    # Creating adflow object
    # CFDSolver = ADFLOW(options=solverOptions, comm=comm)

    # Adding pressure distribution output
    # if slice:
    #     CFDSolver.addSlices("z", 0.5, sliceType="absolute")

    ############## Run CFD
    # CFDSolver(ap)

    ############## Evaluating objectives
    # funcs = {}
    # CFDSolver.evalFunctions(ap, funcs)
    # CFDSolver.checkSolutionFailure(ap, funcs)

    ############# Post-processing

    # printing the result
    if comm.rank == 0:
        # print("")
        # print("#" + "-"*129 + "#")
        # print(" "*59 + "Result" + ""*59)
        # print("#" + "-"*129 + "#")
        # print("")

        # output = {}

        # # Printing and storing results based on evalFuncs in aero problem
        # for obj in ap.evalFuncs:
        #     print("{} = ".format(obj), funcs["{}_{}".format(ap.name, obj)])
        #     output["{}".format(obj)] = funcs["{}_{}".format(ap.name, obj)]

        # # Other mandatory outputs
        # print("fail = ", funcs["fail"])
        # output["fail"] = funcs["fail"]

        # # Storing the results in output file
        # filehandler = open("output.pickle", "xb")
        # pickle.dump(output, filehandler)
        # filehandler.close()

        # Redirecting to original stdout
        os.dup2(stdout, 1)
        os.close(stdout)

except Exception as e:
    if comm.rank == 0:
        print(e)

finally:
    # close the file
    if comm.rank == 0:
        log.close()

    # Getting intercomm and disconnecting
    # Otherwise, program will enter deadlock
    parent_comm.Disconnect()
