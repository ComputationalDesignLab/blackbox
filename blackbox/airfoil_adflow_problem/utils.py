import os, psutil
from baseclasses import AeroProblem
from pydantic import BaseModel, ConfigDict, Field, model_validator, computed_field
from typing import Literal

class AirfoilADflowOptions(BaseModel):
    """
        This class is used to define various settings for running ADflow simulations on airfoil geometries parameterized using CST
    """

    model_config = ConfigDict(
        extra='forbid',
        arbitrary_types_allowed=True
    )

    # --------------------
    # Mandatory arguments
    # --------------------
    
    airfoil_file: str = Field(
        description="Path to the coordinate file used to define the airfoil geometry"
    )

    solver_options: dict = Field(
        description="Dictionary of ADflow solver options passed directly to the solver"
    )

    aero_problem: AeroProblem = Field(
        exclude=True,
        description="AeroProblem instance defining flow conditions for the ADflow solver, this field is excluded from serialization"
    )

    # --------------
    # CST arguments
    # --------------
    
    num_cst_upper: int = Field(
        default=6,
        description="Number of CST coefficients used to parameterize the upper surface of the airfoil",
        ge=1
    )

    num_cst_lower: int = Field(
        default=6,
        description="Number of CST coefficients used to parameterize the lower surface of the airfoil",
        ge=1
    )

    # ------------------
    # Meshing arguments
    # ------------------

    mesh_levels: int = Field(
        default=129,
        description="Number of levels to march in the normal direction from the airfoil surface",
        ge=1
    )

    initial_offwall_spacing: float = Field(
        default=1e-6,
        description="Height of the first mesh layer adjacent to the airfoil surface",
        ge=0.0
    )

    constant_offwall_layers: int = Field(
        default=1,
        description="Number of initial layers with initial offwall spacing before geometric growth is applied",
        ge=1
    )

    marching_distance: float = Field(
        default=100.0,
        description="Total distance to march along the normal distance for mesh generation, usually expressed as a multiple of airfoil chord length",
        ge=0.0
    )

    refine_volume_mesh: Literal[-2, -1, 0, 1, 2] = Field(
        default=0,
        description="Flag indicating whether volume mesh should be refined/coarsed after mesh generation. The possible values are:"
                    "-2: coarsen the mesh twice"
                    "-2: coarsen the mesh once"
                    "0: no change in generated mesh"
                    "1: refine the mesh once"
                    "2: refine the mesh twice"
    )

    # ----------------------------
    # Simulation output arguments
    # ----------------------------

    scalar_outputs: list[Literal["cl", "clp", "clv", "cd", "cdp", "cdv", "cm"]] = Field(
        default_factory=lambda: ["cl", "clp", "clv", "cd", "cdp", "cdv", "cm"],
        description="List of scalar output identifiers to be written during the simulation, the provided list should contain atleast one valid identifier",
        min_length=1
    )

    surface_outputs: list[Literal["cp", "vx", "vy", "cf", "cfx", "cfy", "yplus"]] = Field(
        default_factory=lambda: ["cp", "vx", "vy", "cf", "cfx", "cfy", "yplus"],
        description="List of surface output identifiers to be written during the simulation, the provided list should contain atleast one valid identifier",
        min_length=1
    )

    volume_outputs: list[Literal["cp", "mach"]] = Field(
        default_factory=lambda: ["cp", "mach"],
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

    write_slice_file: bool = Field(
        default=False,
        description="If True, writes a slice file for post-processing and visualization"
    )

    write_airfoil_coordinates: bool = Field(
        default=False,
        description="If True, writes a deformed airfoil coordinates in a dat file"
    )

    plot_airfoil: bool = Field(
        default=False,
        description="If True, generate a plot of the deformed airfoil geometry vs original file"
    )

    # ------------------------
    # Alpha related arguments
    # ------------------------

    alpha: Literal["explicit", "implicit"] = Field(
        default="explicit",
        description="Angle-of-attack control mode. Use `implicit` if you want solver to solve for alpha based on a target lift coefficient. `NOTE`: If you use `implicit` mode, then `alpha` cannot be set as a parameter"
    )

    target_CL: float = Field(
        default=0.824,
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

    compute_area: bool = Field(
        default=False,
        description="Flag to compute area and include it in scalar outputs"
    )

    @model_validator(mode="after")
    def post_init_validation(self):
        """
            Method for validating and setting some options after initialization
        """

        # set paths
        self.airfoil_file = os.path.abspath(self.airfoil_file)
        self.directory = os.path.abspath(self.directory)

        # check if airfoil file path exists
        if not os.path.exists(self.airfoil_file):
            raise ValueError("Provided airfoil file path does not exist")

        # replace "cm" in scalar outputs list with "cmz"
        if "cm" in self.scalar_outputs:
            self.scalar_outputs.remove("cm")
            self.scalar_outputs.append("cmz")

        # set some solver and mesh options
        self.solver_options["volumeVariables"] = self.volume_outputs
        self.solver_options["surfaceVariables"] = self.surface_outputs
        self.solver_options["writeSurfaceSolution"] = self.write_surface_output
        self.solver_options["writeVolumeSolution"] = self.write_volume_output
        self.solver_options["liftIndex"] = 2 # y-axis
        self.solver_options["printAllOptions"] = False
        self.solver_options["printIntro"] = False
        self.solver_options["outputDirectory"] = "."
        self.solver_options["numberSolutions"] = False
        self.solver_options["printTiming"] = False
        self.solver_options["gridFile"] = "vol_mesh.cgns"

        return self

    @computed_field
    @property
    def meshing_options(self) -> dict:
        """
            Meshing options dictionary derived from mesh-related parameters
        """

        return {
            "inputFile": "surf_mesh.xyz",
            "unattachedEdgesAreSymmetry": False,
            "outerFaceBC": "farfield",
            "BC": {1: {"jLow": "zSymm", "jHigh": "zSymm"}},
            "families": "wall",
            "N": self.mesh_levels,
            "s0": self.initial_offwall_spacing,
            "marchDist": self.marching_distance,
            "nConstantStart": self.constant_offwall_layers
        }
