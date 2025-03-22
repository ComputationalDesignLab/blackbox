############## Script file for running airfoil analysis using DAfoam
# Imports
import pickle, os
from mpi4py import MPI
from pyhyp import pyHyp
from cgnsutilities.cgnsutilities import readGrid
import openmdao.api as om
from dafoam.mphys import DAFoamBuilder
from mphys import Multipoint
from mphys.scenario_aerodynamic import ScenarioAerodynamic
import numpy as np

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
    hyp.writeCGNS("volMesh.cgns")
    # hyp.writePlot3D("volMesh.xyz")

    ############## Refining the mesh

    # Only one processor has to do this
    if comm.rank == 0:

        # Read the grid
        grid = readGrid("volMesh.cgns")

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

        grid.writePlot3d("volMesh.xyz")
        os.system("rm volMesh.cgns")

        if comm.rank == 0:
            print("#" + "-"*129 + "#")
            print(" "*59 + "OpenFoam meshing Log" + ""*59)
            print("#" + "-"*129 + "#")
            print("")

        # Converting plot3d vol mesh to openfoam friendly vol mesh
        os.system("plot3dToFoam -noBlank volMesh.xyz")
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

    # Top class to setup the optimization problem
    class Top(Multipoint):

        def setup(self):

            # create the builder to initialize the DASolvers
            dafoam_builder = DAFoamBuilder(solverOptions, scenario="aerodynamic")
            dafoam_builder.initialize(self.comm)

            # # add the design variable component to keep the top level design variables
            self.add_subsystem("dvs", om.IndepVarComp(), promotes=["*"])

            # add the mesh component
            # self.add_subsystem("mesh", dafoam_builder.get_mesh_coordinate_subsystem())

            # # add the geometry component (FFD)
            # self.add_subsystem("geometry", OM_DVGEOCOMP(file="FFD/wingFFD.xyz", type="ffd"))

            # add a scenario (flow condition) for optimization, we pass the builder
            # to the scenario to actually run the flow and adjoint
            self.mphys_add_scenario("scenario1", ScenarioAerodynamic(aero_builder=dafoam_builder))

            # # need to manually connect the x_aero0 between the mesh and geometry components
            # # here x_aero0 means the surface coordinates of structurally undeformed mesh
            # self.connect("mesh.x_aero0", "geometry.x_aero_in")
            # # need to manually connect the x_aero0 between the geometry component and the scenario1
            # # scenario group
            # self.connect("geometry.x_aero0", "scenario1.x_aero")

            # Connecting mesh and scenario
            # self.connect("mesh.x_aero0", "scenario1.x_aero")

            # add the design variables to the dvs component's output
            self.dvs.add_output("patchV", val=np.array([238.0, ap.alpha]))

            # manually connect the dvs output to the geometry and scenario1
            self.connect("patchV", "scenario1.patchV")

        # def configure(self):

            # define the design variables to the top level
            # self.add_design_var("shape", lower=-1.0, upper=1.0, scaler=10.0)

            # here we fix the U0 magnitude and allows the aoa to change
            # self.add_design_var("patchV", lower=[U0, 0.0], upper=[U0, 10.0], scaler=0.1)

            # # add objective and constraints to the top level
            # self.add_objective("scenario1.aero_post.CD", scaler=1.0)
            # self.add_constraint("scenario1.aero_post.CL", equals=CL_target, scaler=1.0)
            # self.add_constraint("geometry.thickcon", lower=0.5, upper=3.0, scaler=1.0)
            # self.add_constraint("geometry.volcon", lower=1.0, scaler=1.0)
            # self.add_constraint("geometry.rcon", lower=0.8, scaler=1.0)
    
    # OpenMDAO setup
    prob = om.Problem()
    prob.model = Top()
    prob.setup()
    om.n2(prob, show_browser=False, outfile="mphys.html")

    prob.run_model()

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
