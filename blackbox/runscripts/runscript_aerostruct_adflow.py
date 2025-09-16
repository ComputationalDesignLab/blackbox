# Importing classes requried for openmdao and mphys
import openmdao.api as om
from mphys.multipoint import Multipoint
from  mphys.scenario_aerostructural import ScenarioAeroStructural
from tacs.mphys import TacsBuilder
from funtofem.mphys import MeldBuilder
from adflow.mphys import ADflowBuilder

# Importing other python packages
from baseclasses import AeroProblem
from mpi4py import MPI
import numpy as np
import pickle, os

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
    aeroSolverOptions = input["aeroSolverOptions"]
    tacsProblemSetup = input["tacsProblemSetup"]
    tacsElementCallback = input["tacsElementCallback"]
    storeFieldData = input["storeFieldData"]
    spanIndex = input["spanIndex"]

    # Assigning non-shape DVs
    if "alpha" in input.keys():
        ap.alpha = input["alpha"][0]

    if "mach" in input.keys():
        ap.mach = input["mach"][0]

    if "altitude" in input.keys():
        ap.altitude = input["altitude"][0]

    ############## Setting up the openmdao model

    class Top(Multipoint):

        def setup(self):

            # aero builder
            aero_builder = ADflowBuilder(aeroSolverOptions, scenario="Aerostructural")
            aero_builder.initialize(comm)
            self.add_subsystem("mesh_aero", aero_builder.get_mesh_coordinate_subsystem())

            # struct builder
            struct_builder = TacsBuilder(mesh_file="wingbox.bdf",
                element_callback=tacsElementCallback,
                problem_setup=tacsProblemSetup
            )
            self.add_subsystem("mesh_struct", struct_builder.get_mesh_coordinate_subsystem())

            # transfer builder
            if spanIndex == "j":
                isym = 1 # y-axis
            elif spanIndex == "k":
                isym = 2 # z-axis

            ldxfer_builder = MeldBuilder(aero_builder, struct_builder, isym=isym)
            ldxfer_builder.initialize(self.comm)

            # coupled aerostructural scenario
            nonlinear_solver = om.NonlinearBlockGS(maxiter=25, iprint=2, use_aitken=True, rtol=1e-12, atol=1e-12)
            linear_solver = om.LinearBlockGS(maxiter=25, iprint=2, use_aitken=True, rtol=1e-12, atol=1e-12)
            self.mphys_add_scenario(
                "scenario",
                ScenarioAeroStructural(
                    aero_builder=aero_builder,
                    struct_builder=struct_builder,
                    ldxfer_builder=ldxfer_builder,
                ),
                nonlinear_solver,
                linear_solver,
            )

        # def configure(self):
            
        #     # Adds linear and nonlinear solver defined in setup
        #     super().configure()

        #     # You need to add name while adding dv. It is the same name used for output/input.
        #     if "aoa" in input.keys():
        #         ap.addDV("alpha", name="aoa", units="deg")
        #     if "mach" in input.keys():
        #         ap.addDV("mach", name="mach")
        #     if "altitude" in input.keys():
        #         ap.addDV("altitude", name="altitude", units="m")

        #     self.scenario.coupling.aero.mphys_set_ap(ap)
        #     self.scenario.aero_post.mphys_set_ap(ap)

        #     if "aoa" in input.keys():
        #         self.connect("aoa", ["scenario.coupling.aero.aoa", "scenario.aero_post.aoa"])

        #     if "mach" in input.keys():
        #         self.connect("mach", ["scenario.coupling.aero.mach", "scenario.aero_post.mach"])

        #     if "altitude" in input.keys():
        #         self.connect("altitude", ["scenario.coupling.aero.altitude", "scenario.aero_post.altitude"])


    # OpenMDAO setup
    prob = om.Problem(comm=comm)
    prob.model = Top()
    prob.setup(mode="rev")
    om.n2(prob, show_browser=False, outfile="mphys.html")
    
    ############## Run CFD
    prob.run_model()

    # printing the result
    if comm.rank == 0:

        print("")
        print("#" + "-"*129 + "#")
        print(" "*59 + "Result" + ""*59)
        print("#" + "-"*129 + "#")
        print("")

        output = {}
    
        for value in prob.model.list_outputs(val=False, out_stream=None):
            if value[0] == "scenario.struct_post.eval_funcs.ks_vmfailure":
                prob.model.objectives.append("failure")
            if value[0] == "scenario.struct_post.mass_funcs.mass":
                prob.model.objectives.append("mass")

        print(prob.get_val("mesh_struct.fea_mesh.x_struct0", get_remote=False))

        for value in prob.model.objectives:
            if "cl" == value:
                print("cl = ", prob["scenario.aero_post.cl"])
                output["cl"] = prob["scenario.aero_post.cl"]
            if "cd" == value:
                print("cd = ", prob["scenario.aero_post.cd"])
                output["cd"] = prob["scenario.aero_post.cd"]
            if "failure" == value:
                print("failure = ", prob["scenario.struct_post.eval_funcs.ks_vmfailure"])
                output["failure"] = prob["scenario.struct_post.eval_funcs.ks_vmfailure"]
            if "mass" == value:
                print("mass = ", prob["scenario.struct_post.mass_funcs.mass"])
                output["mass"] = prob["scenario.struct_post.mass_funcs.mass"]

        
        
        filehandler = open("output.pickle", "xb")
        pickle.dump(output, filehandler)
        filehandler.close()

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
