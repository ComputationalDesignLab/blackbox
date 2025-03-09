# Imports
import os
import numpy as np
from mpi4py import MPI
from pygeo import DVGeometry
from prefoil import Airfoil
from prefoil.utils import readCoordFile
from smt.sampling_methods import LHS
from ...base.airfoilbaseclass import AirfoilBaseClass, DefaultOptions

# Trying to import pyvista
try:
    import pyvista
except ImportError:
    msg_pyvista = "pyVista is not installed"
else:
    msg_pyvista = None

# Trying to import matplotlib
try:
    import matplotlib.pyplot as plt
except ImportError:
    msg_matplotlib = "Matplotlib is not installed"
else:
    msg_matplotlib = None

comm = MPI.COMM_WORLD

class AirfoilFFD(AirfoilBaseClass):
    """
        This class provides methods for generating samples for a general airfoil
        using FFD parameterization.
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
        requiredOptions = ["airfoilFile", "nffd"]

        # Validating user provided options
        self._checkOptions(defaultOptions, requiredOptions, options)

        # Updating/Appending the default option list with user provided options
        self._setOptions(options)

        # Overiding/set some solver options
        self.options["solverOptions"]["printAllOptions"] = False
        self.options["solverOptions"]["printIntro"] = False
        self.options["solverOptions"]["outputDirectory"] = "."
        self.options["solverOptions"]["numberSolutions"] = False

        # Raise an error if pyvista is not installed
        if self.options["getFlowFieldData"]:
            if msg_pyvista != None:
                self._error(msg_pyvista)
            else:
                self.pyvista = pyvista

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

        # Read the coordinate file
        self.origCoords = readCoordFile(self.options["airfoilFile"])
        airfoil = Airfoil(self.origCoords)

        # Creating FFD box
        airfoil.generateFFD(nffd=int(self.options["nffd"]/2), filename=directory + "/ffd", fitted=self.options["fitted"], 
                            xmargin=self.options["xmargin"], ymarginu=self.options["ymarginu"], 
                            ymarginl=self.options["ymarginl"], coords=self.origCoords)

        # Creating DVGeometry object
        self.DVGeo = DVGeometry(directory + "/ffd.xyz")
        self.parametrization = "FFD"

        # Adding pointset to the parametrization
        self.origCoords = np.hstack(( self.origCoords, np.zeros((self.origCoords.shape[0], 1)) ))
        self.DVGeo.addPointSet(self.origCoords, "airfoil")

        # Some initializations which will be used later
        self.DV = []
        self.genSamples = 0

    # ----------------------------------------------------------------------------
    #                       Design Variable related methods
    # ----------------------------------------------------------------------------

    def addDV(self, name: str, lowerBound: list, upperBound: list) -> None:
        """
            Method for adding a DV for CST parameterization.

            Parameters
            ----------
            name: str
                Name of the DV. It can be "shape", "alpha", "mach" or "altitude".

            lowerBound: float or 1D numpy array
                Lower bound of the DV

            upperBound: float or 1D numpy array
                Upper bound of the DV
        """

        # Checking
        self._checkDV(name, lowerBound, upperBound)

        if name == "shape":
            locator = np.array(["{}".format(name)]*len(lowerBound))

            if len(self.DV) == 0:
                self.upperBound = upperBound
                self.lowerBound = lowerBound
                self.locator = locator
            else:
                self.upperBound = np.append(self.upperBound, upperBound)
                self.lowerBound = np.append(self.lowerBound, lowerBound)
                self.locator = np.append(self.locator, locator)

            # Adding FFD points as a DV
            self.DVGeo.addSpanwiseLocalDV("shape", spanIndex="k", axis="y")
            
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

        # Adding the DV to the list
        self.DV.append(name)

        # Limits for sampler
        xlimits = np.hstack((self.lowerBound.reshape(-1,1), self.upperBound.reshape(-1,1)))

        # Creating the sampler
        self.sampler = LHS(xlimits=xlimits, criterion=self.options["samplingCriterion"], random_state=self.options["randomState"])

    def removeDV(self, name: str) -> None:
        """
            Method to remove a DV. 

            Parameters
            ----------
            name: str
                Name of the DV to be removed. It can be "shape", "alpha", "mach" or "altitude".
        """

        if name not in self.DV:
            self._error("{} doesn't exists as a DV.".format(name))

        # Finding the indices to be removed
        loc = self.locator == name
        loc = loc.reshape(-1,)

        # Removing the entry in the bounds and locator
        self.lowerBound = np.delete(self.lowerBound, loc)
        self.upperBound = np.delete(self.upperBound, loc)
        self.locator = np.delete(self.locator, loc)

        # Removing the entry from DV list
        self.DV.remove(name)

        if len(self.DV) == 0:
            delattr(self, "sampler")
        else:
            # Limits for sampler
            xlimits = np.hstack((self.lowerBound.reshape(-1,1), self.upperBound.reshape(-1,1)))

            # Creating the sampler
            self.sampler = LHS(xlimits=xlimits, criterion=self.options["samplingCriterion"], random_state=self.options["randomState"])

    # ----------------------------------------------------------------------------
    #                       Methods related to validation
    # ----------------------------------------------------------------------------

    def _checkDV(self, name: str, lb, ub) -> None:
        """
            Method for validating DV user wants to add.

            Parameters
            ----------
            name: str
                Name of the DV. It can be "shape", "alpha", "mach" or "altitude".

            lowerBound: float or 1D numpy array
                Lower bound of the DV
            
            upperBound: float or 1D numpy array
                Upper bound of the DV
        """

        # List of possible DVs
        possibleDVs = ["shape", "alpha", "mach", "altitude"]

        # Validating name of the DV
        if not isinstance(name, str):
            self._error("Name argument is not a string.")

        # Checking if the DV is allowed
        if name not in possibleDVs:
            self._error("{} argument is not a valid DV.".format(name))

        # Checking if the DV is already added
        if name in self.DV:
            self._error("{} already exists.".format(name))

        # Checking if alpha can be added as a DV
        if name == "alpha":
            if self.options["alpha"] != "explicit":
                self._error("Alpha cannot be a design variable when \"alpha\" attribute in option is \"implicit\".")

        # Checking if these variables are initialized through aero problem or not
        if name in ["mach", "altitude"]:
            if name not in self.options["aeroProblem"].inputs.keys():
                self._error("You need to initialize \"{}\" in the aero problem to set it as design variable.".format(name))

        # Validating the bounds for "shape" variable
        if name == "shape":
            if not isinstance(lb, np.ndarray) or lb.ndim != 1:
                self._error("Lower bound for \"shape\" variable should be a 1D numpy array.")

            if not isinstance(ub, np.ndarray) or ub.ndim != 1:
                self._error("Upper bound for \"shape\" variable should be a 1D numpy array.")

            if self.options["fixLETE"]:
                if len(lb) != self.options["nffd"] - 2:
                    self._error("Length of lower bound array is not equal to (nffd - 2) points.")

                if len(ub) != self.options["nffd"] - 2:
                    self._error("Length of upper bound array is not equal to (nffd - 2) points.")
            else:
                if len(lb) != self.options["nffd"]:
                    self._error("Length of lower bound array is not equal to number of FFD points.")

                if len(ub) != self.options["nffd"]:
                    self._error("Length of upper bound array is not equal to number of FFD points.")

            if np.any(lb >= ub):
                self._error("Lower bound is greater than or equal to upper bound for atleast one DV.")

            # Checking if the bounds are within the limits ########## Instead of checking for bounds, implement checking for intersecting surface
            # coeff = self.DVGeo.origFFDCoef
            # index = self.DVGeo.getLocalIndex(0)
            # dist = coeff[index[:,1,0], 1] - coeff[index[:,0,0], 1]
            # allowableLowerBound = np.zeros(self.options["nffd"])
            # allowableUpperBound = np.zeros(self.options["nffd"])

            # for i in range(dist.shape[0]):
            #     allowableLowerBound[2*i] = -0.45 * dist[i]
            #     allowableLowerBound[2*i+1] = -0.45 * dist[i]
            #     allowableUpperBound[2*i] = 0.45 * dist[i]
            #     allowableUpperBound[2*i+1] = 0.45 * dist[i]

            # if np.any(lb <= allowableLowerBound):
            #     self._error("Lower bound for some FFD points is greater than or equal to 45% of the \
            #                 local FFD thickness. Reduce the bound and try again.")
                
            # if np.any(ub >= allowableUpperBound):
            #     self._error("Upper bound for some FFD points is greater than or equal to 45% of the \
            #                 local FFD thickness. Reduce the bound and try again.")

        else:
            if not isinstance(lb, float):
                self._error("Lower Bound argument is not a float.")

            if not isinstance(ub, float):
                self._error("Upper Bound argument is not a float.")

            if lb >= ub:
                self._error("Lower bound is greater than or equal to upper bound.")

    # ----------------------------------------------------------------------------
    #                   Methods related to Laplacian Smoothing
    # ----------------------------------------------------------------------------

    def LaplacianSmoothing(self, x: np.ndarray) -> np.ndarray:
        """
            Method for performing Laplacian smoothing on the FFD points.

            Parameters
            ----------
            x: 1D numpy array
                Array containing the design variable values.

            Returns
            -------
            x_smooth: 1D numpy array
                Array containing the smoothed design variable values.
        """

        # Performing checks
        if len(self.DV) == 0:
            self._error("Add design variables before running the analysis.")

        if not isinstance(x, np.ndarray):
            self._error("Input sample is not a numpy array.")

        if x.ndim != 1:
            self._error("Input sample is a single dimensional array.")

        if len(x) != len(self.lowerBound):
            self._error("Input sample is not of correct size.")

        # If no geometric design variable is present, then return the original airfoil
        if self.DVGeo.getNDV() == 0:
            self._error("No shape design variable is present to smooth", type=1)
            return x
        
        # Copying the original x
        x_smooth = x.copy()

        # Smoothing parameters
        theta = self.options["smoothingTheta"]
        maxIter = self.options["smoothingMaxIterations"]
        tolerance = self.options["smoothingTolerance"]

        # Creating dictionary from x
        loc = self.locator == "shape"
        loc = loc.reshape(-1,)
        y = x_smooth[loc]

        # Fixing the LE and TE FFD points
        if self.options["fixLETE"]:

            midpoint = y[0]/2
            y[0] -= midpoint
            y = np.append(-midpoint, y)

            midpoint = y[-1]/2
            y[-1] -= midpoint
            y = np.append(y, -midpoint)

        # Array containing index of lower and upper surface FFD points in design variable
        lowerIndex = np.linspace(0, self.options["nffd"]-2, int(self.options["nffd"]/2), dtype=int) # lower surface ffd point index in dv array
        upperIndex = np.linspace(1, self.options["nffd"]-1, int(self.options["nffd"]/2), dtype=int) # upper surface ffd point index in dv array

        # Looping control variables
        error = 1
        itr = 0

        # Smoothing operation
        while error > tolerance and itr < maxIter:

            # Getting the airfoil coordinates for previous iteration
            if self.options["fixLETE"]:
                x_smooth[loc] = np.delete(y, [0,-1]) # Dropping first and last entry, since they are not design variables
            else:
                x_smooth[loc] = y
            airfoil_prev = self.getAirfoil(x_smooth)

            # Lower surface FFD points
            for i in lowerIndex:
                if i == lowerIndex[0]:
                    y[i] = theta * y[i] + (1-theta) * y[i+2]

                elif i == lowerIndex[-1]:
                    y[i] = theta * y[i] + (1-theta) * y[i-2]

                else:
                    y[i] = theta * y[i] + (1-theta) * (y[i-2] + y[i+2])/2

            # Upper surface FFD points
            for i in upperIndex:
                if i == upperIndex[0]:
                    y[i] = theta * y[i] + (1-theta) * y[i+2]

                elif i == upperIndex[-1]:
                    y[i] = theta * y[i] + (1-theta) * y[i-2]

                else:
                    y[i] = theta * y[i] + (1-theta) * (y[i-2] + y[i+2])/2

            # Getting the airfoil coordinates for current iteration
            if self.options["fixLETE"]:
                x_smooth[loc] = np.delete(y, [0,-1])  # Dropping first and last entry, since they are not design variables
            else:
                x_smooth[loc] = y
            airfoil = self.getAirfoil(x_smooth)

            # Calculating error
            error = np.linalg.norm(airfoil - airfoil_prev) / np.linalg.norm(airfoil_prev)
            
            itr = itr + 1

        return x_smooth
