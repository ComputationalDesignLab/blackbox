import os
from typing import Optional
from dataclasses import dataclass, field
from baseclasses import AeroProblem

@dataclass
class AirfoilOptions:

    # Parameterization, solver, and meshing options
    solver: str = "adflow"
    solver_options: dict = field(default_factory=dict)
    meshing_options: dict = field(default_factory=dict)
    aero_problem: Optional[AeroProblem] = None
    openfoam_directory: str = "."
    refine: int = 0
    num_cst_upper: int = 6
    num_cst_lower: int = 6
    num_ffd: int = 20
    get_flowfield_data: bool = False
    region: str = "surface"

    # Other options
    airfoil_file: str = "airfoil.dat"
    directory: str = "output"
    num_processors: int = 4

    # Writing and plotting options
    write_slice_file: bool = False
    write_airfoil_coordinates: bool = False
    write_deformed_ffd: bool = False
    plot_airfoil: bool = False

    # Implicit alpha options
    alpha: str = "explicit"
    target_CL: float = 0.824
    target_CL_tol: float = 1e-4
    starting_alpha: float = 2.5

    # FFD options
    fitted: bool = False
    xmargin: bool = 0.001
    ymarginu: bool = 0.02
    ymarginl: bool = 0.02
    fix_LETE: bool = True

    # Smoothing options
    smoothing: bool = False
    smoothing_theta: float = 0.75
    smoothing_max_iter: int = 100
    smoothing_tolerance: float = 5e-4

    def __post_init__(self):
        """
            Method for validating options after initialization
        """

        try:
            assert os.path.exists(os.path.abspath(self.airfoil_file)), "Hello"

        except Exception as e:
            print(e)
            exit()

        else:
            self.airfoil_file = os.path.abspath(self.airfoil_file)

        
