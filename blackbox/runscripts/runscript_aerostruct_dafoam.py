# Importing classes requried for openmdao and mphys
import pickle, os
import openmdao.api as om
import numpy as np
from mphys import Multipoint
from  mphys.scenario_aerostructural import ScenarioAeroStructural
from tacs.mphys import TacsBuilder
from funtofem.mphys import MeldBuilder
from dafoam.mphys import DAFoamBuilder
from mpi4py import MPI
from scipy.io import savemat
from cgnsutilities.cgnsutilities import readGrid
import pyvista as pv
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
    writeForceField = input["writeForceField"]
    writeDisplacementField = input["writeDisplacementField"]
    spanIndex = input["spanIndex"]
    normalAxis = input["liftIndex"]

    # Assigning non-shape DVs
    if "alpha" in input.keys():
        ap.alpha = input["alpha"][0]

    if "mach" in input.keys():
        ap.mach = input["mach"][0]

    if "altitude" in input.keys():
        ap.altitude = input["altitude"][0]

    ############## Convert cgns mesh to openfoam friendly mesh

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

    # Wait till root is done with mesh conversion
    comm.barrier()

    ############## Setting up dafoam

    # Parameters for DAfoam options
    U0 = float(ap.V)
    rho0 = float(ap.rho)
    p0 = float(ap.P)
    T0 = float(ap.T)
    mu0 = float(ap.mu)
    nuTilda0 = 4.5e-5
    A0 = ap.areaRef
    L0 = ap.chordRef

    # Some DAfoam options
    aeroSolverOptions["printDAOptions"] = False
    aeroSolverOptions["function"] = {}
    
    for obj in ap.evalFuncs:

        if obj == "cd":
            
            aeroSolverOptions["function"]["CD"] = {}
            aeroSolverOptions["function"]["CD"]["type"] = "force"
            aeroSolverOptions["function"]["CD"]["source"] = "patchToFace"
            aeroSolverOptions["function"]["CD"]["patches"] = aeroSolverOptions["designSurfaces"]
            aeroSolverOptions["function"]["CD"]["directionMode"] = "parallelToFlow"
            aeroSolverOptions["function"]["CD"]["patchVelocityInputName"] = "patchV"
            aeroSolverOptions["function"]["CD"]["scale"] = 1.0 / (0.5 * U0 * U0 * A0 * rho0)

        elif obj == "cl":
            
            aeroSolverOptions["function"]["CL"] = {}
            aeroSolverOptions["function"]["CL"]["type"] = "force"
            aeroSolverOptions["function"]["CL"]["source"] = "patchToFace"
            aeroSolverOptions["function"]["CL"]["patches"] = aeroSolverOptions["designSurfaces"]
            aeroSolverOptions["function"]["CL"]["directionMode"] = "normalToFlow"
            aeroSolverOptions["function"]["CL"]["patchVelocityInputName"] = "patchV"
            aeroSolverOptions["function"]["CL"]["scale"] = 1.0 / (0.5 * U0 * U0 * A0 * rho0)

        elif obj == "cmz":

            aeroSolverOptions["function"]["CMZ"] = {}
            aeroSolverOptions["function"]["CMZ"]["type"] = "moment"
            aeroSolverOptions["function"]["CMZ"]["source"] = "patchToFace"
            aeroSolverOptions["function"]["CMZ"]["patches"] = aeroSolverOptions["designSurfaces"]
            aeroSolverOptions["function"]["CMZ"]["axis"] = [0.0, 0.0, 1.0]
            aeroSolverOptions["function"]["CMZ"]["center"] = [ap.xRef, ap.yRef, ap.zRef]
            aeroSolverOptions["function"]["CMZ"]["scale"] = 1.0 / (0.5 * U0 * U0 * A0 * L0 * rho0)

        elif obj == "cmy":

            aeroSolverOptions["function"]["CMY"] = {}
            aeroSolverOptions["function"]["CMY"]["type"] = "moment"
            aeroSolverOptions["function"]["CMY"]["source"] = "patchToFace"
            aeroSolverOptions["function"]["CMY"]["patches"] = aeroSolverOptions["designSurfaces"]
            aeroSolverOptions["function"]["CMY"]["axis"] = [0.0, 1.0, 0.0]
            aeroSolverOptions["function"]["CMY"]["center"] = [ap.xRef, ap.yRef, ap.zRef]
            aeroSolverOptions["function"]["CMY"]["scale"] = 1.0 / (0.5 * U0 * U0 * A0 * L0 * rho0)

    aeroSolverOptions["normalizeStates"] = { # internal, only required for adjoint solver
        "U": U0,
        "p": p0,
        "T": T0,
        "nuTilda": nuTilda0 * 10.0,
        "phi": 1.0,
    }

    aeroSolverOptions["primalBC"] = { # internal
        "U0": {"variable": "U", "patches": ["inout"], "value": [U0, 0.0, 0.0]},
        "p0": {"variable": "p", "patches": ["inout"], "value": [p0]},
        "T0": {"variable": "T", "patches": ["inout"], "value": [T0]},
        "nuTilda0": {"variable": "nuTilda", "patches": ["inout"], "value": [nuTilda0]},
        "useWallFunction": True,
        "thermo:mu": mu0,
    }

    if "alpha" in input.keys():

        aeroSolverOptions["inputInfo"] = { # internal - used in mphys
            "aero_vol_coords": {"type": "volCoord", "components": ["solver", "function"]},
            "patchV": {
                "type": "patchVelocity",
                "patches": ["inout"],
                "flowAxis": "x",
                "normalAxis": normalAxis,
                "components": ["solver", "function"],
            },
        }

    aeroSolverOptions["outputInfo"] = { # internal - used in  mphys
        "f_aero": {
            "type": "forceCouplingOutput",
            "patches": aeroSolverOptions["designSurfaces"],
            "components": ["forceCoupling"],
            "pRef": p0
        },
    }

    # point and normal for the symmetry plane
    if spanIndex == "j":
        symmetryPlanes = [[[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]]
        isym = 1 # y-axis
    elif spanIndex == "k":
        symmetryPlanes = [[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]]
        isym = 2 # z-axis

    mesh_options = {
        "gridFile": os.getcwd(),
        "fileType": "OpenFOAM",
        "symmetryPlanes": symmetryPlanes
    }

    ############## Setting up the openmdao model

    class Top(Multipoint):

        def setup(self):

            # add the design variable component to keep the top level design variables
            self.dvs = self.add_subsystem("dvs", om.IndepVarComp(), promotes=["*"])

            # add the design variables to the dvs component's output
            self.dvs.add_output("patchV", val=np.array([U0, ap.alpha]))

            # aero builder
            aero_builder = DAFoamBuilder(aeroSolverOptions, scenario="Aerostructural", mesh_options=mesh_options)
            aero_builder.initialize(comm)
            self.add_subsystem("mesh_aero", aero_builder.get_mesh_coordinate_subsystem())

            # struct builder
            struct_builder = TacsBuilder(
                mesh_file="wingbox.bdf",
                element_callback=element_callback,
                problem_setup=problem_setup,
                write_solution=False
            )
            struct_builder.initialize(self.comm)
            self.add_subsystem("mesh_struct", struct_builder.get_mesh_coordinate_subsystem())

            ldxfer_builder = MeldBuilder(aero_builder, struct_builder, isym=isym)
            ldxfer_builder.initialize(self.comm)

            # coupled aerostructural scenario
            nonlinear_solver = om.NonlinearBlockGS(maxiter=25, iprint=2, use_aitken=True, rtol=1e-6, atol=1e-6, err_on_non_converge=True)
            linear_solver = om.LinearBlockGS(maxiter=25, iprint=2, use_aitken=True, rtol=1e-6, atol=1e-6, err_on_non_converge=True)
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

            # manually connect the dvs output to the geometry and scenario1
            self.connect("patchV", "scenario.patchV")

    # OpenMDAO setup
    prob = om.Problem(comm=comm)
    prob.model = Top()
    prob.setup(mode="rev")
    om.n2(prob, show_browser=False, outfile="mphys.html")

    ############## Run the model
    try:
        prob.run_model()

    except om.AnalysisError:
        fail = True

    else:

        if prob.model.scenario.coupling.aero.solver.DASolver.primalFail:
            fail = True
        else:
            fail = False

            # Write struct solution files
            prob.model.scenario.coupling.struct.sp.setOption("numbersolutions", False)
            prob.model.scenario.coupling.struct.sp.writeSolution(baseName="struct_output")

            if writeDisplacementField:
                displacement = prob.get_val("scenario.coupling.struct.masker.u_struct_masked", get_remote=True)
                struct_coords = prob.get_val("mesh_struct.fea_mesh.x_struct0", get_remote=True)
                if comm.rank == 0:
                    data = {
                        "displacement": displacement.reshape(-1,6),
                        "struct_coords": struct_coords.reshape(-1,3)
                    }
                    savemat("displacement_field.mat", data)

    # printing the result
    if comm.rank == 0:

        # Reconstruct the field
        os.system("reconstructPar")

        os.system("rm -rf 0")
        # os.system("rm -rf constant")
        os.system("rm -rf system")

        os.system("rm -rf processor*")
        os.system("rm volMesh.xyz")

        # Write force field
        if writeForceField and not fail:

            os.system("touch case.foam") # create a dummy file so that pyvista knows this a openfoam  directory

            reader = pv.OpenFOAMReader("case.foam") # read openfoam results
            
            reader.set_active_time_value(reader.time_values[-1]) # set the latest time

            mesh = reader.read() # read results

            boundaries = mesh['boundary'] # reading only surface data for now

            aero_coords = boundaries[aeroSolverOptions["designSurfaces"][0]].points # extract mesh points

            force_point = boundaries[aeroSolverOptions["designSurfaces"][0]].point_data["volumeForceField"] # extract force field data
            pressure_point = boundaries[aeroSolverOptions["designSurfaces"][0]].point_data["p"].reshape(-1,1) # extract pressure
            cp_point = (pressure_point - p0) / 0.5 / rho0 / U0**2

            force_cell = boundaries[aeroSolverOptions["designSurfaces"][0]].cell_data["volumeForceField"] # extract force field data
            pressure_cell = boundaries[aeroSolverOptions["designSurfaces"][0]].cell_data["p"].reshape(-1,1) # extract pressure
            cp_cell = (pressure_cell - p0) / 0.5 / rho0 / U0**2

            data = {
                "force_point": force_point,
                "cp_point": cp_point,
                "force_cell": force_cell,
                "cp_cell": cp_cell,
                "aero_coords": aero_coords
            }
            savemat("aero_field_data.mat", data)

            os.system("rm case.foam")

        print("")
        print("#" + "-"*129 + "#")
        print(" "*59 + "Result" + ""*59)
        print("#" + "-"*129 + "#")
        print("")

        output = {}

        for obj in aeroSolverOptions["function"].keys():
            print(f"{obj.lower()} = ", prob[f"scenario.aero_post.{obj}"].item())
            output[f"{obj.lower()}"] = prob[f"scenario.aero_post.{obj}"].item()

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
