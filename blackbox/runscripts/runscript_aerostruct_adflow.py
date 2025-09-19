# Importing classes requried for openmdao and mphys
import openmdao.api as om
from mphys.multipoint import Multipoint
from  mphys.scenario_aerostructural import ScenarioAeroStructural
from tacs.mphys import TacsBuilder
from funtofem.mphys import MeldBuilder
from adflow.mphys import ADflowBuilder
from mpi4py import MPI
import pickle, os
from struct_setup_file import * # importing from structural definitions

os.environ['OPENMDAO_REPORTS'] = "0" # disable report generation

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
    sliceLocation = input["sliceLocation"]
    storeFieldData = input["storeFieldData"]
    spanIndex = input["spanIndex"]

    # Assigning non-shape DVs
    if "alpha" in input.keys():
        ap.alpha = input["alpha"][0]

    if "mach" in input.keys():
        ap.mach = input["mach"][0]

    if "altitude" in input.keys():
        ap.altitude = input["altitude"][0]

    structSolverOptions = {
        "numberSolutions": False,
        "outputDir": ".",
        "printTiming": False,
        "printLevel": 1
    }

    ############## Setting up the openmdao model

    class Top(Multipoint):

        def setup(self):

            # aero builder
            aero_builder = ADflowBuilder(aeroSolverOptions, scenario="Aerostructural", write_solution=False)
            aero_builder.initialize(comm)
            self.add_subsystem("mesh_aero", aero_builder.get_mesh_coordinate_subsystem())

            # struct builder
            struct_builder = TacsBuilder(
                mesh_file="wingbox.bdf",
                element_callback=element_callback,
                problem_setup=problem_setup,
                write_solution=False,
                pytacs_options=structSolverOptions
            )
            struct_builder.initialize(self.comm)
            self.add_subsystem("mesh_struct", struct_builder.get_mesh_coordinate_subsystem())

            # transfer builder
            if spanIndex == "j":
                isym = 1 # y-axis
            elif spanIndex == "k":
                isym = 2 # z-axis

            ldxfer_builder = MeldBuilder(aero_builder, struct_builder, isym=isym)
            ldxfer_builder.initialize(self.comm)

            # coupled aerostructural scenario
            nonlinear_solver = om.NonlinearBlockGS(maxiter=25, iprint=2, use_aitken=True, rtol=1e-4, atol=1e-4, err_on_non_converge=True)
            linear_solver = om.LinearBlockGS(maxiter=25, iprint=2, use_aitken=True, rtol=1e-4, atol=1e-4, err_on_non_converge=True)
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

            for discipline in ["aero", "struct"]:
                self.mphys_connect_scenario_coordinate_source(f"mesh_{discipline}", "scenario", discipline)

        def configure(self):

            super().configure()

            self.scenario.coupling.aero.mphys_set_ap(ap)
            self.scenario.aero_post.mphys_set_ap(ap)

    # OpenMDAO setup
    prob = om.Problem(comm=comm)
    prob.model = Top()
    prob.setup(mode="rev")
    om.n2(prob, show_browser=False, outfile="mphys.html")

    # Adding pressure distribution output
    if spanIndex == "j":
        sliceDirection = "y" # y-axis
    elif spanIndex == "k":
        sliceDirection = "z" # z-axis

    for loc in sliceLocation: 
        prob.model.scenario.coupling.aero.solver.solver.addSlices(sliceDirection, loc, sliceType="absolute")

    ############## Run the model
    try:
        prob.run_model()
    except om.AnalysisError:
        fail = True
    else:
        fail = False

        # Write aero solution files
        prob.model.scenario.coupling.aero.solver.solver.writeSolution(baseName="aero_output")

        # Important to wait for all processors to finish before calling disconnect
        # Otherwise, program will enter deadlock
        comm.barrier()

        # Write struct solution files
        prob.model.scenario.coupling.struct.sp.writeSolution(baseName="struct_output")
    
    # Important to wait for all processors to finish before calling disconnect
    # Otherwise, program will enter deadlock
    comm.barrier()

    # printing the result
    if comm.rank == 0:

        print("")
        print("#" + "-"*129 + "#")
        print(" "*59 + "Result" + ""*59)
        print("#" + "-"*129 + "#")
        print("")

        output = {}

        for obj in ap.evalFuncs :
            print(f"{obj} = ", prob[f"scenario.aero_post.{obj}"].item())
            output[f"{obj}"] = prob[f"scenario.aero_post.{obj}"].item()

        for name in prob.model.scenario.struct_post.system_iter():
            for abs_name, meta in name.list_outputs(out_stream=None, prom_name=True):
                print(f"{abs_name} = ", meta["val"].item())
                output[f"{abs_name}"] = meta["val"].item()

        output["fail"] = fail

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
