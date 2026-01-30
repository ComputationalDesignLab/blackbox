import os, psutil
from baseclasses import AeroProblem
from dataclasses import field
from pydantic.dataclasses import dataclass
from pydantic import ConfigDict, Field

valid_scalar_outputs = ["cl", "clp", "clv", "cd", "cdp", "cdv", "cm"]
valid_surface_outputs = ["cp", "vx", "vy", "cf", "cfx", "cfy", "yplus"]
valid_volume_outputs = ["cp", "mach"]

@dataclass(config=ConfigDict(arbitrary_types_allowed=True))
class AirfoilADflowOptions:

    # Parameterization, solver, and meshing options
    airfoil_file: str
    solver_options: dict
    meshing_options: dict
    aero_problem: AeroProblem = Field(exclude=True)
    refine: int = 0
    num_cst_upper: int = 6
    num_cst_lower: int = 6

    # Other options
    directory: str = "output"
    num_processors: int = 4
    scalar_outputs: list[str] = field(default_factory=lambda: valid_scalar_outputs)
    surface_outputs: list[str] = field(default_factory=lambda: valid_surface_outputs)
    volume_outputs: list[str] = field(default_factory=lambda: valid_volume_outputs)

    # Writing and plotting options
    write_slice_file: bool = False
    write_airfoil_coordinates: bool = False
    write_deformed_ffd: bool = False
    write_surface_output: bool = False
    write_volume_output: bool = False
    plot_airfoil: bool = False

    # Implicit alpha options
    alpha: str = "explicit"
    target_CL: float = 0.824
    target_CL_tol: float = 1e-3
    starting_alpha: float = 2.5

    model_config = {"arbitrary_types_allowed": True}

    def __post_init__(self):
        """
            Method for validating and setting some options after initialization
        """

        assert isinstance(self.aero_problem, AeroProblem), "'aero_problem' argument must be an instance of AeroProblem class"
        assert os.path.exists(os.path.abspath(self.airfoil_file)), "provided file path for airfoil does not exists"
        assert self.num_cst_upper > 0, "number of CST coefficients for upper surface should be more than 0"
        assert self.num_cst_lower > 0, "number of CST coefficients for lower surface should be more than 0"
        assert psutil.cpu_count(False) >= self.num_processors + 1, "requested number of processors is more than available processors"

        assert self.refine in [-2,-1,0,1,2], "meshing refine options should be from -2, -1, 0, 1, and 2"
        assert self.alpha in ["explicit", "implicit"], "option 'alpha' should be 'explicit' or 'implicit'"

        # check scalar and surface output list
        for name in self.scalar_outputs:
            assert name in valid_scalar_outputs, f"{name} is not a valid scalar output"
            if name == "cm":
                self.scalar_outputs.remove("cm")
                self.scalar_outputs.append("cmz")

        for name in self.surface_outputs:
            assert name in valid_surface_outputs, f"{name} is not a valid surface output"

        for name in self.volume_outputs:
            assert name in valid_volume_outputs, f"{name} is not a valid volume output"

        if self.plot_airfoil:
            try:
                import matplotlib.pyplot as plt
            except ImportError:
                raise ImportError("'matplotlib' is required for plotting the airfoil")

        # set paths
        self.airfoil_file = os.path.abspath(self.airfoil_file)
        self.directory = os.path.abspath(self.directory)

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
        self.meshing_options["inputFile"] = "surf_mesh.xyz"
