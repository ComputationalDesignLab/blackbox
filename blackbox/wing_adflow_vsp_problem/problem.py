import os, sys, pickle, json, h5py
import numpy as np
import pyvista as pv
from time import time
from packaging import version

from .utils import WingADflowVSPOptions
from ..base_classes.base_problem import BaseProblem
from pygeo import DVGeometryVSP
from cgnsutilities.cgnsutilities import readGrid
from idwarp import USMesh
from mpi4py import MPI
        
comm = MPI.COMM_WORLD

class WingADflowVSP(BaseProblem):

    def __init__(self, options: WingADflowVSPOptions):
        """
            Class for performing wing analysis using ADflow solver
            
            OpenVSP is used for parametrizing the airfoil in this problem

            idwarp is used for deforming the given volume mesh based on the wing geometry changes

            Refer to the documentation for more details about the analysis pipeline

            Parameters
            ----------
            options: WingADflowVSPOptions
                wing vsp option dataclass object
        """

        assert isinstance(options, WingADflowVSPOptions), "options argument should be an object of AirfoilOptions class"

        self.options = options

        # check if required packages are available with specific versions
        try:
            import adflow
        except:
            raise RuntimeError(
                "ADFLOW solver is not installed or can not be imported\n\n"
                "You can follow the installation guide: https://mdolab-adflow.readthedocs-hosted.com/en/latest/install.html"
            )
        else:
            MIN_VERSION = version.parse("2.11.0")
            CURR_VERSION = version.parse(adflow.__version__)
            assert CURR_VERSION >= MIN_VERSION, f"ADFlow version is {CURR_VERSION} but minimum {MIN_VERSION} is required"

        if self.options.write_surface_output or self.options.write_volume_output:
            try:
                import pyvista, h5py
            except ImportError as e:
                raise ValueError(
                    "`pyvista` and `h5py` are not installed, it is required when surface or volume output is set to `True`"
                    "Run this command to install both these pacages: pip install pyvista h5py"
                ) from e

        # Creating directory for storing the results
        if not os.path.isdir(self.options.directory):
            os.system("mkdir {}".format(self.options.directory))
        else:
            os.system("rm -r {}".format(self.options.directory))
            os.system("mkdir {}".format(self.options.directory))

        # Initializing the parametrization object
        self.parametrization = DVGeometryVSP(self.options.vsp_file, comm=comm)

        # Some checks
        geoms = self.parametrization.vspModel.FindGeoms()
        assert(len(geoms) == 1) and self.parametrization.vspModel.GetGeomTypeName(geoms[0]) == "Wing", "Your OpenVSP model should contain only one component of type `WING`"

        # Some initializations which will be used later
        self.parameters = []
        self.mask = np.array([])
        self.bounds = (np.array([]), np.array([]))
        self.samples_generated = 0
        self._perform_init_shape = True

    def add_wing_planform_parameter(self, name: str, lower_bound: float | np.ndarray, upper_bound: float | np.ndarray) -> None:
        """
            Method for adding wing planform parameters such as span, dihedral, sweep

            Parameters
            ----------
            name: str
                name of the parameter to be added. It can only be:
                    * `alpha`: angle of attack of the flow
                    * `mach`: mach number of the flow
                    * `altitude`: altitude number for the flow

            lower_bound: float or np.ndarray
                lower bound for the parameter

            upper_bound: float or np.ndarray
                upper bound for the parameter
        """

        # Perform some steps when this method is called first time
        if self._perform_init_shape:

            check_pygeo_version()

            self._perform_init_shape = False

        pass

    def add_flow_parameter(self, name: str, lower_bound: float, upper_bound: float) -> None:
        """
            Method for adding a flow parameter for the wing problem

            Parameters
            ----------
            name: str
                name of the parameter to be added. It can only be:
                    * `alpha`: angle of attack of the flow
                    * `mach`: mach number of the flow
                    * `altitude`: altitude number for the flow

            lower_bound: float
                lower bound for the parameter

            upper_bound: float
                upper bound for the parameter
        """

        assert isinstance(name, str) and name.lower() in ["alpha", "mach", "altitude"], "`name` must be one of `alpha`, `mach`, and `alititude`"
        assert isinstance(lower_bound, float) and isinstance(upper_bound, float), "Lower and upper bound must be float"
        assert upper_bound > lower_bound, "upper bound must be greater than lower bound"

        mask = [f"{name.lower()}"]
        lb = np.append(self.bounds[0], np.array([lower_bound]))
        ub = np.append(self.bounds[1], np.array([upper_bound]))
            
        self.bounds = (lb, ub)
        self.mask = np.append(self.mask, mask)
        self.parameters.append(name.lower())

    def __call__(self):

        pass

    ###########################################################
    ######### Below methods are for internal use only #########
    ###########################################################

def check_pygeo_version():

    try:
        import pygeo
    except:
        raise RuntimeError(
            "pyGeo is not installed or can not be imported\n\n"
            "You can follow the installation guide: https://mdolab-pygeo.readthedocs-hosted.com/en/latest/install.html"
        )
    else:
        MIN_VERSION = version.parse("1.17.0")
        CURR_VERSION = version.parse(pygeo.__version__)
        assert CURR_VERSION >= MIN_VERSION, f"pygeo version is {CURR_VERSION} but minimum {MIN_VERSION} is required"