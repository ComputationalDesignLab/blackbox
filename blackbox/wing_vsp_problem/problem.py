import os, sys, pickle, json, h5py
import numpy as np
import pyvista as pv
from time import time
from packaging import version

from .utils import WingVSPOptions
from ..base_classes.base_problem import BaseProblem
from pygeo import DVGeometryVSP
from cgnsutilities.cgnsutilities import readGrid
from idwarp import USMesh
from mpi4py import MPI
        
comm = MPI.COMM_WORLD

class WingVSPADflow(BaseProblem):

    def __init__(self, options: WingVSPOptions):
        """
            Class for performing wing analysis using ADflow solver
            
            OpenVSP is used for parametrizing the airfoil in this problem

            idwarp is used for deforming the given volume mesh based on the wing geometry changes

            Refer to the documentation for more details about the analysis pipeline

            Parameters
            ----------
            options: WingVSPOptions
                wing vsp option dataclass object
        """

        assert isinstance(options, WingVSPOptions), "options argument should be an object of AirfoilOptions class"

        self.options = options

        # Creating directory for storing the results
        if not os.path.isdir(self.options.directory):
            os.system("mkdir {}".format(self.options.directory))
        else:
            os.system("rm -r {}".format(self.options.directory))
            os.system("mkdir {}".format(self.options.directory))

        # Initializing the parametrization object
        self.parametrization = DVGeometryVSP(self.options.vsp_file, comm=comm)

        



