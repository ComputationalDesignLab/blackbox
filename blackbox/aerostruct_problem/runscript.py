############# Script file for running aerostruct analysis

# imports
import numpy as np
from scipy.io import savemat
import pickle, os, json, warnings

from mpi4py import MPI
import openmdao.api as om
from idwarp import USMesh
from mphys import Multipoint
from pygeo import DVGeometryVSP
from tacs.mphys import TacsBuilder
from funtofem.mphys import MeldBuilder
from adflow.mphys import ADflowBuilder
from tacs.pymeshloader import pyMeshLoader
from tacs import constitutive, elements, functions
from  mphys.scenario_aerostructural import ScenarioAeroStructural

warnings.filterwarnings(
    "ignore",
    message="Using internally generated IDWarp surfaces.*"
)

def element_callback(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
    """
        Callback function used to setup TACS element objects and DVs
    """

    # Material properties
    rho = 2500.0  # density kg/m^3
    E = 70.0e9  # Young's modulus (Pa)
    nu = 0.30  # Poisson's ratio
    ys = 350e6  # yield stress
    t = 0.01  # shell thickness, m

    # Setup (isotropic) property and constitutive objects
    prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
    # Set one thickness dv for every component
    con = constitutive.IsoShellConstitutive(prop, t=t, tNum=dvNum)

    # For each element type in this component, pass back the appropriate tacs element object
    transform = None
    elem = elements.Quad4Shell(transform, con)

    return elem

# getting MPI comm
comm = MPI.COMM_WORLD
parent_comm = comm.Get_parent()

# send the processor
parent_comm.send(os.getpid(), dest=0, tag=comm.rank)

try:

    # redirecting the stdout
    log = open("log.txt", "a")
    stdout = os.dup(1)
    os.dup2(log.fileno(), 1)

    ############## Reading input file for the analysis

    # reading input file
    filehandler = open("input.pickle", 'rb') 
    input = pickle.load(filehandler)
    filehandler.close()

    # getting some options
    ap = input["aero_problem"]
    scalar_outputs = input["scalar_outputs"]
    aero_solver_options = input["aero_solver_options"]
    vsp_file = input["vsp_file"]
    struct_mesh_file = input["struct_mesh_file"]

    ############## Read and set flow variables

    # reading parameters
    with open("parameters.json") as fp:
        params = json.load(fp)
    fp.close()

    ap.setDesignVars(params)

    ############## Deform the volume mesh

    if vsp_file is not None:

        options = {
            'gridFile': aero_solver_options["gridFile"],
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

        # create the mesh object
        mesh_deformer = USMesh(options=options, comm=comm)

        # extract all original surface mesh coordinates
        orig_surface_mesh_coords = mesh_deformer.getSurfaceCoordinates()

        # initialize the bdf mesh reader and read the bdf file
        struct_mesh = pyMeshLoader(comm, False)
        struct_mesh.scanBdfFile(struct_mesh_file)

        # get nastran object
        nastran_obj = struct_mesh.getBDFInfo()

        # collect structural nodes from BDF to create a 3D point cloud
        node_ids = []
        orig_struct_mesh_coords = np.zeros((0,3))

        for nid, node in nastran_obj.nodes.items():
            node_ids.append(nid)
            orig_struct_mesh_coords = np.vstack((orig_struct_mesh_coords, node.get_position()))

        # create pygeo object
        geo_vsp = DVGeometryVSP(vsp_file, comm=comm)

        # adding co-ordinates as a pointset
        geo_vsp.addPointSet(orig_surface_mesh_coords, "wing_surface_mesh")
        geo_vsp.addPointSet(orig_struct_mesh_coords, "struct_mesh")

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

        # update the wing surface mesh
        new_surface_mesh_coords = geo_vsp.update("wing_surface_mesh")
        new_struct_mesh_coords = geo_vsp.update("struct_mesh")

        # set updated surface mesh
        mesh_deformer.setSurfaceCoordinates(new_surface_mesh_coords)

        # deform volume mesh
        mesh_deformer.warpMesh()

        # write deformed vol mesh
        mesh_deformer.writeGrid('vol_mesh.cgns')

        # set deformed volume mesh for analysis
        aero_solver_options["gridFile"] = 'vol_mesh.cgns'

        # update the node location in the nastran object
        for (nid, node), xyz in zip(nastran_obj.nodes.items(), new_struct_mesh_coords):
                node.set_position(nastran_obj, xyz)

        # write the new wingbox file
        struct_mesh_file = "wingbox.bdf"
        nastran_obj.write_bdf(struct_mesh_file)

        if input["write_vsp_file"]:
            geo_vsp.writeVSPFile("updated_model.vsp3")

        if input["write_stl_file"]:
            geo_vsp.vspModel.ExportFile("updated_model.stl", geo_vsp.vspModel.SET_ALL, geo_vsp.vspModel.EXPORT_STL)

    ############## Setting up openmdao model

    # direction along wing span
    if "liftIndex" in aero_solver_options.keys():
        if aero_solver_options["liftIndex"] == 2:
            direction = "z"
            isym = 2
        elif aero_solver_options["liftIndex"] == 3:
            direction = "y"
            isym = 1
    else:
        direction = "z"
        isym = 2

    def problem_setup(scenario_name, fea_assembler, problem):
        """
            Function to add fixed forces and eval functions to structural problems used in tacs builder
        """

        # Add TACS Functions
        problem.addFunction("mass", functions.StructuralMass)
        problem.addFunction("ks_vmfailure", functions.KSFailure, safetyFactor=1.0, ksWeight=100.0)

        # Add 1g load
        if direction == "z":
            g = np.array([0.0, -9.81, 0.0])  # gravity is along negative y direction m/s^2
        elif direction == "y":
            g = np.array([0.0, 0.0, -9.81])  # gravity is along negative z direction m/s^2
        
        # Multiply by load factor if it is a parameter
        if "load_factor" in params.keys():
            g = params["load_factor"] * g
            
        problem.addInertialLoad(g)

    class Top(Multipoint):

        def setup(self):

            # aero builder
            aero_builder = ADflowBuilder(aero_solver_options, scenario="Aerostructural", write_solution=False)
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

            # add the design variable component to keep the top level design variables
            flow_parameters = self.add_subsystem("flow_parameters", om.IndepVarComp(), promotes=["*"])

            for name in ap.DVs:
                flow_parameters.add_output(name, np.array([params[name]]))

            # coupled aerostructural scenario
            nonlinear_solver = om.NonlinearBlockGS(maxiter=25, iprint=2, use_aitken=True, rtol=1e-12, atol=1e-12, err_on_non_converge=True)
            linear_solver = om.LinearBlockGS(maxiter=25, iprint=2, use_aitken=True, rtol=1e-12, atol=1e-12, err_on_non_converge=True)
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

            for name in ap.DVs:
                self.connect(name, f"scenario.coupling.aero.{name}")
                self.connect(name, f"scenario.aero_post.{name}")

    ############## OpenMDAO setup

    prob = om.Problem(comm=comm)
    prob.model = Top()
    prob.setup(mode="rev")
    om.n2(prob, show_browser=False, outfile="mphys.html")

    # adding wing slices
    if input["write_slice_file"]:
        prob.model.scenario.coupling.aero.solver.solver.addSlices(positions=input["slice_location"], direction=direction)

    # adding lift distribution
    if input["write_lift_distribution"]:
        prob.model.scenario.coupling.aero.solver.solver.addLiftDistribution(nSegments=input["num_segments"], direction=direction)

    ############## Run the model

    try:
        prob.run_model()

    except:
        funcs = {}
        funcs["fail"] = True
        
    else:
        # evaluating objectives
        funcs = {}
        prob.model.scenario.coupling.aero.solver.solver.evalFunctions(ap, funcs, evalFuncs=scalar_outputs)
        prob.model.scenario.coupling.aero.solver.solver.checkSolutionFailure(ap, funcs)

        # write aero solution files
        prob.model.scenario.coupling.aero.solver.solver.writeSolution()

        # write struct solution files
        prob.model.scenario.coupling.struct.sp.setOption("numbersolutions", False)
        prob.model.scenario.coupling.struct.sp.writeSolution(baseName="struct_output")

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

except Exception as e:
    if comm.rank == 0:
        print(e)

finally:
    # Redirecting to original stdout
    os.dup2(stdout, 1)
    os.close(stdout)

    # close the file
    log.close()

    # Getting intercomm and disconnecting
    # Otherwise, program will enter deadlock
    parent_comm.Disconnect()
