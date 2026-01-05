import os
from .cst import CST
from .utils import AirfoilOptions
from ..base_classes.base_problem import BaseProblem

class AirfoilADflow(BaseProblem):

    def __init__(self, options: AirfoilOptions):
        """
            Class for performing airfoil analysis using ADflow solver
            
            CST is used for parametrizing the airfoil in this problem

            pyHyp is used for creating mesh around the airfoil

            Refer to the documentation for more details about the analysis pipeline
        """

        assert isinstance(options, AirfoilOptions), "options argument should be an object of AirfoilOptions class"

        self.options = options

        # Getting abs path for the storage directory
        self.options.directory = os.path.abspath(self.options.directory)

        # Creating directory for storing the results
        if not os.path.isdir(self.options.directory):
            os.system("mkdir {}".format(self.options.directory))
        else:
            os.system("rm -r {}".format(self.options.directory))
            os.system("mkdir {}".format(self.options.directory))

        # Initializing the parametrization object
        self.parametrization = CST(self.options.airfoil_file, num_cst=[self.options.num_cst_upper, self.options.num_cst_lower])

        # Some initializations which will be used later
        self.design = []
        self.samples_generated = 0

    