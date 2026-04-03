############# Script file for running airfoil analysis.

# Imports
import pickle, os, json
from mpi4py import MPI
from adflow import ADFLOW
from idwarp import USMesh
from pygeo import DVGeometryVSP
import warnings

warnings.filterwarnings(
    "ignore",
    message="Using internally generated IDWarp surfaces.*"
)

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
    scalar_outputs = input["scalar_outputs"]
    solver_options = input["solver_options"]
    vsp_file = input["vsp_file"]

    # implicit/explicit alpha options
    alpha_type = input["alpha_type"]
    CL_target = input["target_CL"]
    target_CL_tol = input["target_CL_tol"]
    starting_alpha = input["starting_alpha"]
    initial_delta_alpha = input["initial_delta_alpha"]
    max_iterations = input["max_iterations"]

    ############## Read and set flow variables

    # reading parameters
    with open("parameters.json") as fp:
        params = json.load(fp)
    fp.close()

    ap.setDesignVars(params)

    ############## Deform the volume mesh

    if vsp_file is not None:

        options = {
            'gridFile': solver_options["gridFile"],
            # 'fileType':'CGNS',
            # 'specifiedSurfaces':None,
            # 'symmetrySurfaces':None,
            # 'symmetryPlanes':[],
            # 'aExp': 3.0,
            # 'bExp': 5.0,
            'LdefFact': 75.0, # needed for large changes in mesh
            # 'alpha': 0.25,
            'errTol': 1e-5,
            # 'evalMode': 'fast',
            # 'useRotations': True,
            # 'zeroCornerRotations': True,
            # 'cornerAngle': 30.0,
            # 'bucketSize': 8,
        }

        # Create the mesh object
        mesh_deformer = USMesh(options=options, comm=comm)

        # Extract all original surface mesh coordinates
        orig_surface_mesh_coords = mesh_deformer.getSurfaceCoordinates()

        # create pygeo object
        geo_vsp = DVGeometryVSP(vsp_file, comm=comm)

        # Adding surface mesh co-ordinates as a pointset
        geo_vsp.addPointSet(orig_surface_mesh_coords, "wing_surface_mesh")

        # reading pygeo parameters
        with open("pygeo_parameters.json") as fp:
            geo_params = json.load(fp)
        fp.close()

        # add pygeo parameters
        for key in geo_params.keys():
            component, group, parm = key.split(":")
            geo_vsp.addVariable(component=component, group=group, parm=parm, scaledStep=False)

        # set value of parameters
        geo_vsp.setDesignVars(geo_params)

        # Update the wing surface mesh
        new_surface_mesh_coords = geo_vsp.update("wing_surface_mesh")

        # set updated surface mesh
        mesh_deformer.setSurfaceCoordinates(new_surface_mesh_coords)

        # deform volume mesh
        mesh_deformer.warpMesh()

        # write deformed vol mesh
        mesh_deformer.writeGrid('vol_mesh.cgns')

        # set deformed volume mesh for analysis
        solver_options["gridFile"] = 'vol_mesh.cgns'

        if input["write_vsp_file"]:
            geo_vsp.writeVSPFile("updated_model.vsp3")

        if input["write_stl_file"]:
            geo_vsp.vspModel.ExportFile("updated_model.stl", geo_vsp.vspModel.SET_ALL, geo_vsp.vspModel.EXPORT_STL)

    ############## Settign up adflow

    if comm.rank == 0:
        print("")
        print("#" + "-"*129 + "#")
        print(" "*59 + "Solver Log" + ""*59)
        print("#" + "-"*129 + "#")
        print("")

    # Creating adflow object
    CFDSolver = ADFLOW(options=solver_options, comm=comm)
    
    # Direction along wing span
    if CFDSolver.options["liftIndex"] == 2:
        direction = "z"
    elif CFDSolver.options["liftIndex"] == 3:
        direction = "y"

    # Adding lift distribution
    if input["write_lift_distribution"]:
        CFDSolver.addLiftDistribution(nSegments=input["num_segments"], direction=direction)

    # Adding wing slices
    if input["write_slice_file"]:
        CFDSolver.addSlices(positions=input["slice_location"], direction=direction)

    ############## Run CFD

    if input["alpha_type"] == "explicit":

        CFDSolver(ap)

        # Evaluating objectives
        funcs = {}
        CFDSolver.evalFunctions(ap, funcs, evalFuncs=scalar_outputs)
        CFDSolver.checkSolutionFailure(ap, funcs)
        
    elif input["alpha_type"] == "implicit":

        # Run CFD
        itr_results = CFDSolver.solveCL(
            ap, 
            CLStar=CL_target,
            alpha0=starting_alpha,
            delta=initial_delta_alpha,
            tol=target_CL_tol, 
            autoReset=False,
            maxIter=max_iterations,
            writeSolution=True
        )

        # Evaluating objectives
        funcs = {}
        CFDSolver.evalFunctions(ap, funcs, evalFuncs=scalar_outputs)

    ############## post-processing
    
    if comm.rank == 0:

        if os.path.exists(f"{ap.name}_surf.cgns"):
            os.rename(f"{ap.name}_surf.cgns", "surface_solution.cgns")

        if os.path.exists(f"{ap.name}_vol.cgns"):
            os.rename(f"{ap.name}_vol.cgns", "volume_solution.cgns")

        if input["write_lift_distribution"]:
            os.rename(f"{ap.name}_lift.dat", "lift_distribution.dat")

        if input["write_slice_file"]:
            os.rename(f"{ap.name}_slices.dat", "slices.dat")

        if alpha_type == "implicit":
            funcs["fail"] = not itr_results["converged"]

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
