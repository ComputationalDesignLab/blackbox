############# Script file for running airfoil analysis.

# Imports
import pickle, os, json
from mpi4py import MPI
from adflow import ADFLOW
from pyhyp import pyHyp

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

    # Getting some options
    ap = input["aero_problem"]
    refine_volume_mesh = input["refine_volume_mesh"]
    scalar_outputs = input["scalar_outputs"]
    solverOptions = input["solver_options"]
    meshingOptions = input["meshing_options"]

    # implicit/explicit alpha options
    alpha_type = input["alpha_type"]
    CL_target = input["target_CL"]
    target_CL_tol = input["target_CL_tol"]
    starting_alpha = input["starting_alpha"]

    # Assigning non-shape DVs
    if "alpha" in input.keys():
        ap.alpha = input["alpha"][0]

    if "mach" in input.keys():
        ap.mach = input["mach"][0]

    if "reynolds" in input.keys():
        ap.reynolds = input["reynolds"][0]

    ############## Generating mesh

    if comm.rank == 0:
        print("#" + "-"*129 + "#")
        print(" "*59 + "Meshing Log" + ""*59)
        print("#" + "-"*129 + "#")
        print("")

    hyp = pyHyp(options=meshingOptions, comm=comm)
    hyp.run()
    hyp.writeCGNS(solverOptions["gridFile"])

    ############## Refining the mesh

    # Only one processor has to do this
    if comm.rank == 0:

        if refine_volume_mesh != 0:
            from cgnsutilities.cgnsutilities import readGrid

            # Read the grid
            grid = readGrid(solverOptions["gridFile"])

            if refine_volume_mesh > 0:
                for _ in range(refine_volume_mesh):
                    grid.refine(["i", "k"])

            else:
                for _ in range(-refine_volume_mesh):
                    grid.coarsen()

            grid.writeToCGNS(solverOptions["gridFile"])

    # Wait till root is done with refining/coarse of mesh
    comm.barrier()
    
    ############## Settign up adflow

    if comm.rank == 0:
        print("")
        print("#" + "-"*129 + "#")
        print(" "*59 + "Solver Log" + ""*59)
        print("#" + "-"*129 + "#")
        print("")

    # Creating adflow object
    CFDSolver = ADFLOW(options=solverOptions, comm=comm)

    # Adding pressure distribution output
    if input["write_slice_file"]:
        CFDSolver.addSlices("z", 0.5, sliceType="absolute")

    ############## Run CFD

    if input["alpha_type"] == "explicit":

        CFDSolver(ap)

        # Evaluating objectives
        funcs = {}
        CFDSolver.evalFunctions(ap, funcs, evalFuncs=scalar_outputs)
        CFDSolver.checkSolutionFailure(ap, funcs)
        
    elif input["alpha_type"] == "implicit":

        # Run CFD
        itr_results = CFDSolver.solveCL(ap, CLStar=CL_target, alpha0=input["starting_alpha"], delta=0.2, tol=input["target_CL_tol"], autoReset=False, maxIter=8, writeSolution=True)

        # Evaluating objectives
        funcs = {}
        CFDSolver.evalFunctions(ap, funcs, evalFuncs=scalar_outputs)

    ############## post-processing
    
    if comm.rank == 0:

        if os.path.exists(f"{ap.name}_surf.cgns"):
            os.rename(f"{ap.name}_surf.cgns", "surface.cgns")

        if os.path.exists(f"{ap.name}_vol.cgns"):
            os.rename(f"{ap.name}_vol.cgns", "volume.cgns")

        if alpha_type == "implicit":
            funcs["fail"] = not itr_results["converged"]

        # rename the pitching moment and change the sign
        if f"{ap.name}_cmz" in funcs.keys():
            funcs[f"{ap.name}_cm"] = -funcs.pop(f"{ap.name}_cmz")

        # remove ap name from keys
        funcs = {
            k[len(f"{ap.name}_"):] if k.startswith(f"{ap.name}_") else k: v
            for k, v in funcs.items()
        }

        # dump the scalar outputs to json file
        with open("scalar_outputs.json", "w") as fp:
            json.dump(funcs, fp, indent=4)
        fp.close()

        print("")
        print("#" + "-"*129 + "#")
        print(" "*59 + "Result" + ""*59)
        print("#" + "-"*129 + "#")
        print("")

        # Printing and storing results based on evalFuncs in aero problem
        for key, value in funcs.items():
            print(f"{key} = {value}")

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
