############## Script file for running airfoil analysis.
# Imports
import pickle, os, h5py
from mpi4py import MPI
from adflow import ADFLOW
from pyhyp import pyHyp
from cgnsutilities.cgnsutilities import readGrid
import numpy as np
import pyvista

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
    ap = input["aero_problem"]
    refine = input["refine"]
    slice = input["write_slice_file"]
    get_flowfield_data = input["get_flowfield_data"]
    CL_target = input["target_CL"]
    target_CL_tol = input["target_CL_tol"]
    starting_alpha = input["starting_alpha"]

    # Assigning non-shape DVs
    if "alpha" in input.keys():
        ap.alpha = input["alpha"][0]

    if "mach" in input.keys():
        ap.mach = input["mach"][0]

    if "altitude" in input.keys():
        ap.altitude = input["altitude"][0]

    # Getting solver and meshing options from input file
    solverOptions = input["solver_options"]
    solverOptions["gridFile"] = "vol_mesh.cgns"
    solverOptions["liftindex"] = 2 # Always 2 since meshing is done internally

    meshingOptions = input["meshing_options"]
    meshingOptions["inputFile"] = "surf_mesh.xyz"

    ############## Generating mesh

    if comm.rank == 0:
        print("#" + "-"*129 + "#")
        print(" "*59 + "Meshing Log" + ""*59)
        print("#" + "-"*129 + "#")
        print("")

    hyp = pyHyp(options=meshingOptions, comm=comm)
    hyp.run()
    hyp.writeCGNS("vol_mesh.cgns")

    ############## Refining the mesh

    # Only one processor has to do this
    if comm.rank == 0:

        # Read the grid
        grid = readGrid("vol_mesh.cgns")

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

        grid.writeToCGNS("vol_mesh.cgns")

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

    ############## Run CFD
    itr_results = CFDSolver.solveCL(ap, CLStar=CL_target, alpha0=starting_alpha, delta=0.2, tol=target_CL_tol, autoReset=False, maxIter=8, writeSolution=True)

    ############## Evaluating objectives
    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)

    ############# Post-processing

    # printing the result
    if comm.rank == 0:
        print("")
        print("#" + "-"*129 + "#")
        print(" "*59 + "Result" + ""*59)
        print("#" + "-"*129 + "#")
        print("")

        # Storing the results in output file
        f = h5py.File('output.hdf5','w')

        scalars = f.create_group("scalars")

        scalars.attrs["fail"] = not itr_results["converged"]

        # Printing and storing results based on evalFuncs in aero problem
        for obj in ap.evalFuncs:
            
            print("{} = ".format(obj), funcs["{}_{}".format(ap.name, obj)])

            scalars.attrs[f"{obj}"] = funcs["{}_{}".format(ap.name, obj)]
            
        if get_flowfield_data:

            field_group = f.create_group("fields")

            reader = pyvista.CGNSReader(f"{ap.name}_surf.cgns")
            
            reader.load_boundary_patch = False

            ds = reader.read() # read the mesh

            str_grid = ds[0][0] # get the base-block

            for var_name in set(ds[0][0].array_names):
                if var_name != "Base/Zone":
                    field_group.create_dataset(var_name.lower(), data=np.asarray(ds[0][0][var_name]))

            if solverOptions["writeSurfaceSolution"]:
                os.system(f"rm {ap.name}_surf.cgns")

        f.close()

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
