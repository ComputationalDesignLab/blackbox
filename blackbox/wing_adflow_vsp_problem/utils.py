import os, psutil
from baseclasses import AeroProblem
from pydantic import BaseModel, ConfigDict, Field, model_validator, computed_field
from typing import Literal

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
    
    vsp_file: str = Field(
        description="Path to the vsp file containing the wing geometry"
    )

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

    @model_validator(mode="after")
    def post_init_validation(self):
        """
            Method for validating and setting some options after initialization
        """

        # set paths
        self.vsp_file = os.path.abspath(self.vsp_file)
        self.mesh_file = os.path.abspath(self.mesh_file)
        self.directory = os.path.abspath(self.directory)

        # check if vsp file path exists
        if not os.path.exists(self.vsp_file):
            raise ValueError("Provided vsp file path does not exist")
        
        # check if mesh file path exists
        if not os.path.exists(self.mesh_file):
            raise ValueError("Provided mesh file path does not exist")

        # set some solver and mesh options
        self.solver_options["volumeVariables"] = self.volume_outputs
        self.solver_options["surfaceVariables"] = self.surface_outputs
        self.solver_options["writeSurfaceSolution"] = self.write_surface_output
        self.solver_options["writeVolumeSolution"] = self.write_volume_output
        # self.solver_options["liftIndex"] = 2 # y-axis
        self.solver_options["printAllOptions"] = False
        self.solver_options["printIntro"] = False
        self.solver_options["outputDirectory"] = "."
        self.solver_options["numberSolutions"] = False
        self.solver_options["printTiming"] = False
        self.solver_options["gridFile"] = self.mesh_file

        return self
