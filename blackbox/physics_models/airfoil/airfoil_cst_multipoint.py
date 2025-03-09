# Imports
import os
import numpy as np
from mpi4py import MPI
from pygeo import DVGeometryCST
from prefoil.utils import readCoordFile
from smt.sampling_methods import LHS
from ...base import AirfoilBaseClass, DefaultOptions

# Trying to import matplotlib
try:
    import matplotlib.pyplot as plt
except ImportError:
    msg_matplotlib = "Matplotlib is not installed"
else:
    msg_matplotlib = None

comm = MPI.COMM_WORLD

class AirfoilCSTMultipoint(AirfoilBaseClass):
    """
        This class provides methods for generating samples for a general airfoil
        using CST parameterization.
    """

    def __init__(self, options):

        # Partial checking of options argument provided by the user.
        if not isinstance(options, dict):
            self._error("The 'options' argument provided is not a dictionary.")
        elif options == {}:
            self._error("The 'options' argument provided is an empty dictionary.")

        # Creating an empty options dictionary
        self.options = {}

        # Setting up default options
        self._getDefaultOptions(DefaultOptions())

        # Setting up the required options list
        defaultOptions = list(self.options.keys())
        requiredOptions = ["airfoilFile", "numCST"]

        # Validating user provided options
        self._checkOptions(defaultOptions, requiredOptions, options)

        # Updating/Appending the default option list with user provided options
        self._setOptions(options)

        # Overiding/set some solver options
        self.options["solverOptions"]["printAllOptions"] = False
        self.options["solverOptions"]["printIntro"] = False
        self.options["solverOptions"]["outputDirectory"] = "."
        self.options["solverOptions"]["numberSolutions"] = False

        # Raise an error if matplotlib is not installed
        if self.options["plotAirfoil"]:
            if msg_matplotlib != None:
                self._error(msg_matplotlib)
            else:
                self.plt = plt

        # Getting abs path for the storage directory
        self.options["directory"] = os.path.abspath(self.options["directory"])

        # Setting up the folder for saving the results
        directory = self.options["directory"]

        # Creating directory for storing the results
        if not os.path.isdir(directory):
            os.system("mkdir {}".format(directory))
        else:
            os.system("rm -r {}".format(directory))
            os.system("mkdir {}".format(directory))

        # Read the coordinate file
        self.origCoords = readCoordFile(self.options["airfoilFile"])

        # Some validation for coordinate file
        if self.origCoords[0,0] != self.origCoords[-1,0]:
            self._error("The X coordinate of airfoil doesn't start and end at same point.")
        elif self.origCoords[0,1] != self.origCoords[-1,1]:
            self._error("The Y coordinate of airfoil doesn't start and end at same point.")

        # Initializing the parametrization object
        self.DVGeo = DVGeometryCST(self.options["airfoilFile"], numCST=self.options["numCST"], comm=comm)
        self.parametrization = "CST"

        # Adding pointset to the parametrization
        self.origCoords = np.hstack(( self.origCoords, np.zeros((self.origCoords.shape[0], 1)) ))
        self.DVGeo.addPointSet(self.origCoords, "airfoil")

        # Checking the number of points at trailing edge for blunt TE
        # Only two are allowed for CST. Otherwise, meshing will have problem.
        if not self.DVGeo.sharp:
            if len(np.where(self.origCoords[1:-1,0] == self.origCoords[0,0])[0]) > 1:
                self._error("There are more than two points in the trailing edge.")

        # Some initializations which will be used later
        self.DV = []
        self.genSamples = 0

        # Adding aeroproblem DVs
        for ap in self.options["aeroProblem"]:

            # Adding DVs
            for dv in ap.DVs.keys():

                # Extracting the lower and upper bound
                lowerBound = ap.DVs[dv].lower
                upperBound = ap.DVs[dv].upper

                # Setting up the locator
                locator = np.array(["{}".format(dv)])

                # Appending the DV details
                if len(self.DV) == 0:
                    self.upperBound = np.array([upperBound])
                    self.lowerBound = np.array([lowerBound])
                    self.locator = np.array([locator])
                else:
                    self.upperBound = np.append(self.upperBound, upperBound)
                    self.lowerBound = np.append(self.lowerBound, lowerBound)
                    self.locator = np.append(self.locator, locator)

                self.DV.append(dv)

        # Limits for sampler
        xlimits = np.hstack((self.lowerBound.reshape(-1,1), self.upperBound.reshape(-1,1)))

        # Creating the sampler
        self.sampler = LHS(xlimits=xlimits, criterion=self.options["samplingCriterion"], random_state=self.options["randomState"])

    # ----------------------------------------------------------------------------
    #                       Design Variable related methods
    # ----------------------------------------------------------------------------

    def addDV(self, name: str, lowerBound: list, upperBound: list) -> None:
        """
            Method for adding ONLY CST related DVs.

            Parameters
            ----------
            name : str
                Name of the DV. It can be "upper" or "lower".

            lowerBound : float or np.ndarray
                Lower bound of the DV.

            upperBound : float or np.ndarray
                Upper bound of the DV.
        """

        # Checking
        self._checkDV(name, lowerBound, upperBound)

        # Adding DV into DVGeo
        if name in ["alpha", "mach", "altitude"]:
            self._error("Cannot use \"addDV\" method for adding {} as design variable".format(name))

        if name == "upper" or name == "lower":
            locator = np.array(["{}".format(name)]*len(lowerBound))

            if len(self.DV) == 0:
                self.upperBound = upperBound
                self.lowerBound = lowerBound
                self.locator = locator
            else:
                self.upperBound = np.append(self.upperBound, upperBound)
                self.lowerBound = np.append(self.lowerBound, lowerBound)
                self.locator = np.append(self.locator, locator)
        else:
            locator = np.array(["{}".format(name)])

            if len(self.DV) == 0:
                self.upperBound = np.array([upperBound])
                self.lowerBound = np.array([lowerBound])
                self.locator = np.array([locator])
            else:
                self.upperBound = np.append(self.upperBound, upperBound)
                self.lowerBound = np.append(self.lowerBound, lowerBound)
                self.locator = np.append(self.locator, locator)    

        self.DVGeo.addDV("{}".format(name), "{}".format(name))

        # Adding the DV to the list
        self.DV.append(name)

        # Limits for sampler
        xlimits = np.hstack((self.lowerBound.reshape(-1,1), self.upperBound.reshape(-1,1)))

        # Creating the sampler
        self.sampler = LHS(xlimits=xlimits, criterion=self.options["samplingCriterion"], random_state=self.options["randomState"])

    def _checkDV(self, name: str, lb, ub) -> None:
        """
            Method for validating DV.

            Parameters
            ----------
            name : str
                Name of the DV.

            lb : float or np.ndarray
                Lower bound of the DV.
            
            ub : float or np.ndarray
                Upper bound of the DV.
        """

        # List of possible DVs
        possibleDVs = ["upper", "lower"]

        # Validating name of the DV
        if not isinstance(name, str):
            self._error("Name argument is not a string.")

        # Checking if the DV is allowed
        if name not in possibleDVs:
            self._error("{} argument is not a valid DV.".format(name))
        
        # Checking if the DV is already added
        if name in self.DV:
            self._error("{} already exists.".format(name))

        # Validating the bounds for "upper" variable
        if name == "upper":
            if not isinstance(lb, np.ndarray) or lb.ndim != 1:
                self._error("Lower bound for \"upper\" variable should be a 1D numpy array.")

            if not isinstance(ub, np.ndarray) or ub.ndim != 1:
                self._error("Upper bound for \"upper\" variable should be a 1D numpy array.")

            if len(lb) != self.options["numCST"][0]:
                self._error("Length of lower bound for \"upper\" variable is not equal to number of CST coeff for upper surface.")

            if len(ub) != self.options["numCST"][0]:
                self._error("Length of upper bound for \"upper\" variable is not equal to number of CST coeff for upper surface.")

        # Validating the bounds for "lower" variable
        elif name == "lower":
            if not isinstance(lb, np.ndarray) or lb.ndim != 1:
                self._error("Lower bound for \"lower\" variable should be a 1D numpy array.")

            if not isinstance(ub, np.ndarray) or ub.ndim != 1:
                self._error("Upper bound for \"lower\" variable should be a 1D numpy array.")

            if len(lb) != self.options["numCST"][1]:
                self._error("Length of lower bound for \"lower\" variable is not equal to number of CST coeff for lower surface.")

            if len(ub) != self.options["numCST"][1]:
                self._error("Length of upper bound for \"lower\" variable is not equal to number of CST coeff for lower surface.")

        # Validating lb and ub of the scalar DVs
        else:
            if not isinstance(lb, float):
                self._error("Lower Bound argument is not a float.")

            if not isinstance(ub, float):
                self._error("Upper Bound argument is not a float.")

            if lb >= ub:
                self._error("Lower bound is greater than or equal to upper bound.")
