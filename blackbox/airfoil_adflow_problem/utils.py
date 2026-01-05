import os
from typing import Optional
# from dataclasses import dataclass, field
# from baseclasses import AeroProblem
from pydantic.dataclasses import dataclass

@dataclass
class AirfoilOptions:

    # Parameterization, solver, and meshing options
    airfoil_file: str
    solver_options: dict
    meshing_options: dict
    # aero_problem: Optional[AeroProblem]
    refine: int = 0
    num_cst_upper: int = 6
    num_cst_lower: int = 6

    # Other options
    directory: str = "output"
    num_processors: int = 4
    get_flowfield_data: bool = False
    region: str = "surface"

    # Writing and plotting options
    write_slice_file: bool = False
    write_airfoil_coordinates: bool = False
    write_deformed_ffd: bool = False
    plot_airfoil: bool = False

    # Implicit alpha options
    alpha: str = "explicit"
    target_CL: float = 0.824
    target_CL_tol: float = 1e-3
    starting_alpha: float = 2.5

    # # FFD options
    # fitted: bool = False
    # xmargin: bool = 0.001
    # ymarginu: bool = 0.02
    # ymarginl: bool = 0.02
    # fix_LETE: bool = True

    # # Smoothing options
    # smoothing: bool = False
    # smoothing_theta: float = 0.75
    # smoothing_max_iter: int = 100
    # smoothing_tolerance: float = 5e-4

    def __post_init__(self):
        """
            Method for validating options after initialization
        """

        assert os.path.exists(os.path.abspath(self.airfoil_file)), "provided file path for airfoil does not exists"
        assert self.num_cst_upper > 0, "number of CST coefficients for upper surface should be more than 0"
        assert self.num_cst_lower > 0, "number of CST coefficients for lower surface should be more than 0"
        # assert psutil.cpu_count(False) >= self.num_processors + 1, "requested number of processors is more than available processors"

        assert self.refine in [-2,-1,0,1,2], "meshing refining options should be from -2, -1, 0, 1, and 2"

        assert self.region in ["surface", "field"], "region for field extraction should be either 'surface' or 'field'"

        assert self.alpha in ["explicit", "implicit"], "option 'alpha' should be 'explicit' or 'implicit'"

        

        self.airfoil_file = os.path.abspath(self.airfoil_file)

        self.directory = os.path.abspath(self.directory)
