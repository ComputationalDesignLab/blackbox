############## Script file for running airfoil analysis using DAfoam
# Imports
import pickle, os
from mpi4py import MPI
from cgnsutilities.cgnsutilities import readGrid
import openmdao.api as om
from dafoam.mphys import DAFoamBuilder
from mphys import Multipoint
from mphys.scenario_aerodynamic import ScenarioAerodynamic
import numpy as np

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
    getSurfaceForces = input["getSurfaceForces"]
    normalAxis = input["liftIndex"]

    # Assigning non-shape DVs
    if "alpha" in input.keys():
        ap.alpha = input["alpha"][0]

    if "mach" in input.keys():
        ap.mach = input["mach"][0]

    if "altitude" in input.keys():
        ap.altitude = input["altitude"][0]

    # Getting solver and meshing options from input file
    solverOptions = input["solverOptions"]

    ############## Refining the mesh

    # Only one processor has to do this
    if comm.rank == 0:

        # Read the cgns file and write plot3d file
        grid = readGrid("volMesh.cgns")
        grid.writePlot3d("volMesh.xyz")

        if comm.rank == 0:
            print("#" + "-"*129 + "#")
            print(" "*59 + "OpenFoam meshing Log" + ""*59)
            print("#" + "-"*129 + "#")
            print("")

        # Converting plot3d vol mesh to openfoam friendly vol mesh
        os.system("plot3dToFoam -noBlank volMesh.xyz")
        os.system("autoPatch 45 -overwrite")
        os.system("createPatch -overwrite")
        os.system("renumberMesh -overwrite")

    # Wait till root is done with refining/coarse of mesh
    comm.barrier()

    ############## Setting up dafoam

    # Parameters for DAfoam options
    U0 = float(ap.V)
    rho0 = float(ap.rho)
    p0 = float(ap.P)
    T0 = float(ap.T)
    nuTilda0 = 4.5e-5
    A0 = ap.areaRef
    L0 = ap.chordRef

    # Some DAfoam options
    solverOptions["printDAOptions"] = False
    solverOptions["function"] = {}
    
    for obj in ap.evalFuncs:

        if obj == "cd":
            
            solverOptions["function"]["CD"] = {}
            solverOptions["function"]["CD"]["type"] = "force"
            solverOptions["function"]["CD"]["source"] = "patchToFace"
            solverOptions["function"]["CD"]["patches"] = solverOptions["designSurfaces"]
            solverOptions["function"]["CD"]["directionMode"] = "parallelToFlow"
            solverOptions["function"]["CD"]["patchVelocityInputName"] = "patchV"
            solverOptions["function"]["CD"]["scale"] = 1.0 / (0.5 * U0 * U0 * A0 * rho0)

        elif obj == "cl":
            
            solverOptions["function"]["CL"] = {}
            solverOptions["function"]["CL"]["type"] = "force"
            solverOptions["function"]["CL"]["source"] = "patchToFace"
            solverOptions["function"]["CL"]["patches"] = solverOptions["designSurfaces"]
            solverOptions["function"]["CL"]["directionMode"] = "normalToFlow"
            solverOptions["function"]["CL"]["patchVelocityInputName"] = "patchV"
            solverOptions["function"]["CL"]["scale"] = 1.0 / (0.5 * U0 * U0 * A0 * rho0)

        elif obj == "cmz":

            solverOptions["function"]["CMZ"] = {}
            solverOptions["function"]["CMZ"]["type"] = "moment"
            solverOptions["function"]["CMZ"]["source"] = "patchToFace"
            solverOptions["function"]["CMZ"]["patches"] = solverOptions["designSurfaces"]
            solverOptions["function"]["CMZ"]["axis"] = [0.0, 0.0, 1.0]
            solverOptions["function"]["CMZ"]["center"] = [ap.xRef, ap.yRef, ap.zRef]
            solverOptions["function"]["CMZ"]["scale"] = 1.0 / (0.5 * U0 * U0 * A0 * L0 * rho0)

        elif obj == "cmy":

            solverOptions["function"]["CMY"] = {}
            solverOptions["function"]["CMY"]["type"] = "moment"
            solverOptions["function"]["CMY"]["source"] = "patchToFace"
            solverOptions["function"]["CMY"]["patches"] = solverOptions["designSurfaces"]
            solverOptions["function"]["CMY"]["axis"] = [0.0, 1.0, 0.0]
            solverOptions["function"]["CMY"]["center"] = [ap.xRef, ap.yRef, ap.zRef]
            solverOptions["function"]["CMY"]["scale"] = 1.0 / (0.5 * U0 * U0 * A0 * L0 * rho0)

    solverOptions["normalizeStates"] = { # internal
        "U": U0,
        "p": p0,
        "T": T0,
        "nuTilda": nuTilda0 * 10.0,
        "phi": 1.0,
    }

    solverOptions["primalBC"] = { # internal
        "U0": {"variable": "U", "patches": ["inout"], "value": [U0, 0.0, 0.0]},
        "p0": {"variable": "p", "patches": ["inout"], "value": [p0]},
        "T0": {"variable": "T", "patches": ["inout"], "value": [T0]},
        "nuTilda0": {"variable": "nuTilda", "patches": ["inout"], "value": [nuTilda0]},
        "useWallFunction": True,
    }

    if "alpha" in input.keys():

        solverOptions["inputInfo"] = { # internal
            "patchV": {
                "type": "patchVelocity",
                "patches": ["inout"],
                "flowAxis": "x",
                "normalAxis": normalAxis,
                "components": ["solver", "function"],
            },
        }

    # Top class to setup the optimization problem
    class Top(Multipoint):

        def setup(self):

            # create the builder to initialize the DASolver
            dafoam_builder = DAFoamBuilder(solverOptions, scenario="aerodynamic")
            dafoam_builder.initialize(comm)

            # add the design variable component to keep the top level design variables
            self.add_subsystem("dvs", om.IndepVarComp(), promotes=["*"])

            # add a scenario (flow condition) for 
            self.mphys_add_scenario("scenario1", ScenarioAerodynamic(aero_builder=dafoam_builder))

            # add the design variables to the dvs component's output
            self.dvs.add_output("patchV", val=np.array([U0, ap.alpha]))

            # manually connect the dvs output to the geometry and scenario1
            self.connect("patchV", "scenario1.patchV")
    
    # OpenMDAO setup
    prob = om.Problem(comm=comm)
    prob.model = Top()
    prob.setup(mode="rev")
    om.n2(prob, show_browser=False, outfile="mphys.html")

    ############## Run CFD
    prob.run_model()

    ############# Post-processing

    # printing the result
    if comm.rank == 0:

        # Reconstruct the field
        os.system("reconstructPar")

        print("")
        print("#" + "-"*129 + "#")
        print(" "*59 + "Result" + ""*59)
        print("#" + "-"*129 + "#")
        print("")

        output = {}

        # Printing and storing results based on evalFuncs in aero problem
        for obj in solverOptions["function"].keys():
            print(f"{obj.lower()} = ", prob.get_val(f"scenario1.aero_post.{obj}"))
            output[f"{obj.lower()}"] = prob.get_val(f"scenario1.aero_post.{obj}")

        # Other mandatory outputs - TO DO
        print("fail = ", prob.model.scenario1.coupling.solver.DASolver.primalFail)
        output["fail"] = prob.model.scenario1.coupling.solver.DASolver.primalFail

        # Storing the results in output file
        filehandler = open("output.pickle", "xb")
        pickle.dump(output, filehandler)
        filehandler.close()

except Exception as e:
    if comm.rank == 0:
        print(e)

finally:
    # close the file
    if comm.rank == 0:

        # os.system("rm -rf processor*")
        # os.system("rm -rf runscript_out")
        # os.system("rm volMesh.xyz")

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
