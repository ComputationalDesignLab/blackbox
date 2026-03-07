import os, sys, pickle, json, h5py
import numpy as np
import pyvista as pv
from time import time
from packaging import version
from .options import WingADflowVSPOptions
from ..base_classes.base_problem import BaseProblem

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

        # Some initializations which will be used later
        self.parameters = []
        self.mask = np.array([])
        self.bounds = (np.array([]), np.array([]))
        self.samples_generated = 0
        self._perform_init_shape = True

    def __call__(self):

        pass

    ###########################################################
    ######### Below methods are for internal use only #########
    ###########################################################
