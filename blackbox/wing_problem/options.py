import os, psutil
from baseclasses import AeroProblem
from pydantic import BaseModel, ConfigDict, Field, model_validator
from typing import Literal
from .wing_vsp import WingVSP

class WingADflowVSPOptions(BaseModel):
    """
        This class is used to define various settings for running ADflow simulations on wing geometries parameterized using OpenVSP
    """

    model_config = ConfigDict(
        extra='forbid',
        arbitrary_types_allowed=True
    )

    # --------------------
    # Mandatory arguments
    # --------------------

    mesh_file: str = Field(
        description="Path to the cgns volume mesh file for the wing geometry"
    )

    solver_options: dict = Field(
        description="Dictionary of ADflow solver options passed directly to the solver"
    )

    aero_problem: AeroProblem = Field(
        exclude=True,
        description="AeroProblem instance defining flow conditions for the ADflow solver, this field is excluded from serialization"
    )

    wing_vsp: WingVSP = Field(
        exclude=True,
        default=None,
        description="an instance of `WingVSP` class defining various shape parameters"
    )

    # ----------------------------
    # Simulation output arguments
    # ----------------------------

    scalar_outputs: list[Literal["cl", "clp", "clv", "cd", "cdp", "cdv", "cmx", "cmy", "cmz"]] = Field(
        default_factory=lambda: ["cl", "clp", "clv", "cd", "cdp", "cdv", "cmx", "cmy", "cmz"],
        description="List of scalar output identifiers to be written during the simulation, the provided list should contain atleast one valid identifier",
        min_length=1
    )

    surface_outputs: list[Literal["rho", "p", "temp", "cp", "vx", "vy", "vz", "cf", "cfx", "cfy", "cfz", "yplus"]] = Field(
        default_factory=lambda: ["rho", "p", "temp", "cp", "vx", "vy", "vz", "cf", "cfx", "cfy", "cfz", "yplus"],
        description="List of surface output identifiers to be written during the simulation, the provided list should contain atleast one valid identifier",
        min_length=1
    )

    volume_outputs: list[Literal["cp", "mach", "temp"]] = Field(
        default_factory=lambda: ["cp", "mach", "temp"],
        description="List of volume output identifiers to be written during the simulation, the provided list should contain atleast one valid identifier",
        min_length=1
    )

    # -------------------------------
    # Writing and plotting arguments
    # -------------------------------

    write_surface_output: bool = Field(
        default=False,
        description="If True, write surface output from adflow in cgns format. A ready-to-use HDF5 file containing surface data extracted from cgns file is also written"
    )

    write_volume_output: bool = Field(
        default=False,
        description="If True, write volume output from adflow in cgns format. A ready-to-use HDF5 file containing volume data extracted from cgns file is also written"
    )

    write_vsp_file: bool = Field(
        default=False,
        description="If true, write the updated VSP model in vsp3 file. This option will only work when shape variables are added"
    )

    write_stl_file: bool = Field(
        default=False,
        description="If true, write the updated VSP model in a stl file. This option will only work when shape variables are added"
    )

    write_lift_distribution: bool = Field(
        default=False,
        description="If true, add a lift distribution along the wing span and write the distribution file after analysis"
    )

    num_segments: int = Field(
        default=200,
        description="number of points to use along the wing span while writing lift distribution file, only relevant when `write_lift_distribution` is set to True",
        ge=1
    )

    write_slice_file: bool = Field(
        default=False,
        description="If true, add slices along the wing span and write the slice file after analysis"
    )

    slice_location: list = Field(
        default_factory=lambda: [0.05, 0.2, 0.4, 0.6, 0.8, 0.95],
        description="a list describing relative slice location along the wing span, only relevant when `write_slice_file` is set to True",
        min_length=1
    )

    # ------------------------
    # Alpha related arguments
    # ------------------------

    alpha: Literal["explicit", "implicit"] = Field(
        default="explicit",
        description="Angle-of-attack control mode. Use `implicit` if you want solver to solve for alpha based on a target lift coefficient. `NOTE`: If you use `implicit` mode, then `alpha` cannot be set as a parameter"
    )

    target_CL: float = Field(
        default=0.5,
        description="Target lift coefficient used when alpha is solved implicitly"
    )

    target_CL_tol: float = Field(
        default=1e-3,
        description="Tolerance on the lift coefficient when solving for implicit alpha"
    )

    starting_alpha: float = Field(
        default=2.5,
        description="Initial guess for the angle of attack (in degrees) when using implicit alpha mode"
    )

    initial_delta_alpha: float = Field(
        default=0.5,
        description="Initial guess for the change in alpha for secant search"
    )

    max_iterations: int = Field(
        default=8,
        description="number of iterations for secant search"
    )

    # ----------------
    # Other arguments
    # ----------------

    directory: str = Field(
        default="output",
        description="Directory where all simulation outputs will be written"
    )

    num_processors: int = Field(
        default=4,
        description="Number of processors to use for the ADflow simulation",
        ge=1,
        le=psutil.cpu_count(False)-1
    )

    @model_validator(mode="after")
    def post_init_validation(self):
        """
            Method for validating and setting some options after initialization
        """

        # set paths
        self.mesh_file = os.path.abspath(self.mesh_file)
        self.directory = os.path.abspath(self.directory)
        
        # check if mesh file path exists
        if not os.path.exists(self.mesh_file):
            raise ValueError("Provided mesh file path does not exist")
        
        for val in self.slice_location:
            assert isinstance(val, float) and 0.0 <= val <= 1.0, "entries in slice location list must be a float between 0 and 1"

        # set some solver and mesh options
        self.solver_options["volumeVariables"] = self.volume_outputs
        self.solver_options["surfaceVariables"] = self.surface_outputs
        self.solver_options["writeSurfaceSolution"] = self.write_surface_output
        self.solver_options["writeVolumeSolution"] = self.write_volume_output
        self.solver_options["printAllOptions"] = False
        self.solver_options["printIntro"] = False
        self.solver_options["outputDirectory"] = "."
        self.solver_options["numberSolutions"] = False
        self.solver_options["printTiming"] = False
        self.solver_options["gridFile"] = self.mesh_file

        return self
