# Importing classes requried for openmdao and mphys
import openmdao.api as om
from mphys import Multipoint
# from mphys.scenarios import ScenarioAeroStructural
from mphys.scenario_aerostructural import ScenarioAeroStructural
from tacs.mphys import TacsBuilder
from funtofem.mphys import MeldBuilder
from adflow.mphys import ADflowBuilder
from tacs import constitutive, elements, functions

# Importing other python packages
from baseclasses import AeroProblem
from mpi4py import MPI
import numpy as np
import os

os.environ['OPENMDAO_REPORTS'] = "0" # disable report generation

# Getting MPI comm
comm = MPI.COMM_WORLD

# Callback function used to setup TACS element objects and DVs
def element_callback(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):

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

def problem_setup(scenario_name, fea_assembler, problem):
    """
    Helper function to add fixed forces and eval functions
    to structural problems used in tacs builder
    """

    problem.addFunction("mass", functions.StructuralMass) # wingbox mass
    problem.addFunction("ks_vmfailure", functions.KSFailure, safetyFactor=1.0, ksWeight=100.0) # ks aggregated failure criteria

    # Add inertial load
    g = np.array([0.0, -9.81, 0.0])  # m/s^2, -ve y points into ground
    problem.addInertialLoad(g) # cruise

solverOptions = {
    "gridfile": "wing_volMesh.cgns",
    "printAllOptions": False,
    "printIntro": False,
    "printTiming": False,
    "forcesAsTractions": False,
    # Common Parameters
    "monitorvariables": ["cl", "cd", "yplus"],
    "writeTecplotSurfaceSolution": True,
    "writeSurfaceSolution": False,
    "writeVolumeSolution": False,
    "numberSolutions": False,
    # Physics Parameters
    "equationType": "RANS",
    "smoother": "DADI",
    "MGCycle": "sg",
    "nsubiterturb": 10,
    "nCycles": 7000,
    "liftIndex": 2,
    # ANK Solver Parameters
    "useANKSolver": True,
    "ANKSubspaceSize": 400,
    "ANKASMOverlap": 3,
    "ANKPCILUFill": 4,
    "ANKJacobianLag": 5,
    "ANKOuterPreconIts": 3,
    "ANKInnerPreconIts": 3,
    # NK Solver Parameters
    "useNKSolver": True,
    "NKSwitchTol": 1e-6,
    "NKSubspaceSize": 400,
    "NKASMOverlap": 3,
    "NKPCILUFill": 4,
    "NKJacobianLag": 5,
    "NKOuterPreconIts": 3,
    "NKInnerPreconIts": 3,
    # Termination Criteria
    "L2Convergence": 1e-14
}

ap = AeroProblem(name="wing", alpha=2.5, mach=0.85, altitude=10000, areaRef=45.5, chordRef=3.56, evalFuncs=["cl", "cd", "cmz"])

class Top(Multipoint):

    def setup(self):

        # aero builder
        aero_builder = ADflowBuilder(solverOptions, scenario="Aerostructural", write_solution=False)
        aero_builder.initialize(self.comm)
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

        ldxfer_builder = MeldBuilder(aero_builder, struct_builder, isym=2)
        ldxfer_builder.initialize(self.comm)

        # coupled aerostructural scenario
        nonlinear_solver = om.NonlinearBlockGS(maxiter=25, iprint=2, use_aitken=True, rtol=1e-0, atol=1e-0, err_on_non_converge=True)
        linear_solver = om.LinearBlockGS(maxiter=25, iprint=2, use_aitken=True, rtol=1e-0, atol=1e-0, err_on_non_converge=True)
        self.mphys_add_scenario(
            "scenario",
            ScenarioAeroStructural(
                aero_builder=aero_builder,
                struct_builder=struct_builder,
                ldxfer_builder=ldxfer_builder
            ),
            nonlinear_solver,
            linear_solver,
        )

        for discipline in ["aero", "struct"]:
            self.mphys_connect_scenario_coordinate_source(f"mesh_{discipline}", "scenario", discipline)

        # self.connect("dv_struct", f"scenario.dv_struct")

    def configure(self):

        super().configure()

        self.scenario.coupling.aero.mphys_set_ap(ap)
        self.scenario.aero_post.mphys_set_ap(ap)

# OpenMDAO setup
prob = om.Problem(comm=comm)
prob.model = Top()
prob.setup(mode="rev")
om.n2(prob, show_browser=False, outfile="mphys.html")

# ############## Run CFD
try:
    prob.run_model()

except om.AnalysisError:
    fail = True
    
else:
    fail = False

    # Write aero solution files
    prob.model.scenario.coupling.aero.solver.solver.writeSolution(baseName="aero_output", number=None)

    # Write struct solution files
    prob.model.scenario.coupling.struct.sp.setOption("numbersolutions", False)
    prob.model.scenario.coupling.struct.sp.writeSolution(baseName="struct_output")

    prob.model.scenario.coupling.aero.solver.solver.writeMeshFile("adflow_mesh.cgns")

    force = prob.get_val("scenario.coupling.aero.force.f_aero", get_remote=True).reshape(-1,3) 

    coords = prob.get_val("scenario.coupling.geo_disp.x_aero", get_remote=True).reshape(-1,3)

    disp = prob.get_val("scenario.coupling.struct.masker.u_struct_masked", get_remote=True).reshape(-1,6)

    prob.model.scenario.coupling.aero.solver.solver.writeMeshFile("adflow_mesh.cgns")

    if comm.rank == 0:

        print(force.shape)

        print(coords.shape)

        print(np.unique(force, axis=0).shape)

        print(np.unique(coords, axis=0).shape)

        print(disp.shape)

# if comm.rank == 0:

#     for obj in ap.evalFuncs :
#         print(f"{obj} = ", prob[f"scenario.aero_post.{obj}"])

#     for name in prob.model.scenario.struct_post.system_iter():
#         for abs_name, meta in name.list_outputs(out_stream=None, prom_name=True):
#             print(abs_name, meta['val'])

#     print(prob.model.scenario.coupling.aero.aero_solver.curAP.fatalFail)

#     print(fail)
