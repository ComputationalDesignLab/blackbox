# Importing python packages
import os, psutil, shutil, sys, pickle, time
import numpy as np
from pyDOE2 import lhs, fullfact
from scipy.io import savemat
from smt.sampling_methods import LHS
from mpi4py import MPI
from baseclasses import AeroProblem
from pygeo import DVGeometry
from pygeo.geo_utils.polygon import volumeTriangulatedMesh
from pygeo.geo_utils.file_io import readPlot3DSurfFile
from cgnsutilities.cgnsutilities import readGrid
from idwarp import USMesh

class DefaultOptions():
    """
        Class creates a default option list (for solvers and other settings)
        which is later edited/appended with user provided options.
    """

    def __init__(self):

        # Solver Options
        self.aeroSolver = "adflow"
        self.openfoamDir = "."
        self.structSolver = "tacs"

        # Other options
        self.directory = "output"
        self.noOfProcessors = 4
        self.sliceLocation = [] # defines slice location
        self.writeDeformedFFD = False
        self.storeFieldData = False
        self.storeMDAHistory = False

        # Sampling options
        self.samplingCriterion = "cm"
        self.randomState = None

class AeroStructFFD():
    """
        This class provides methods for generating data for a static
        aeroelastic problem using FFD parameterization.

        Currently, only flow and oml variables can be added for analysis
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
        self._getDefaultOptions()

        # Setting up the required options list
        requiredOptions = ["solverOptions", "gridFile", "structMeshFile", "ffdFile", "liftIndex", "aeroProblem", "tacsProblemSetup", "tacsElementCallback"]

        # Validating user provided options
        self._checkOptions(options, requiredOptions)

        # Updating/Appending the default option list with user provided options
        self._setOptions(options)

        # Getting abs path for the directory/files
        self.options["directory"] = os.path.abspath(self.options["directory"])
        self.options["gridFile"] = os.path.abspath(self.options["gridFile"])
        self.options["ffdFile"] = os.path.abspath(self.options["ffdFile"])
        self.options["structMeshFile"] = os.path.abspath(self.options["structMeshFile"])
        if self.options["aeroSolver"] == "dafoam":
            self.options["openfoamDir"] = os.path.abspath(self.options["openfoamDir"])

        # Setting up the folder for saving the results
        directory = self.options["directory"]

        # Creating directory for storing the results
        if not os.path.isdir(directory):
            os.system("mkdir {}".format(directory))
        else:
            os.system("rm -r {}".format(directory))
            os.system("mkdir {}".format(directory))

        # Create mesh deformation object
        self.mesh = USMesh(options={"gridFile": self.options["gridFile"]})

        # Get the surface mesh coordinates
        surfMesh = self.mesh.getSurfaceCoordinates()

        # Creating DVGeometry object
        self.DVGeo = DVGeometry(self.options["ffdFile"])

        # Number of FFD points
        self.nffd = self.DVGeo.getLocalIndex(0).flatten().shape[0]

        # Adding surface mesh co-ordinates as a pointset
        self.DVGeo.addPointSet(surfMesh, "wing_surface_mesh")

        if self.options["liftIndex"] == 2: # y
            self.spanIndex = "k" # If y is lift index, then span is along z (k)
            self.liftIndex = "y"
        elif self.options["liftIndex"] == 3: # z
            self.spanIndex = "j" # If z is lift index, then span is along y (j)
            self.liftIndex = "z"

        # Create reference axis
        self.nRefAxPts = self.DVGeo.addRefAxis("wingAxis", xFraction=0.25, alignIndex=self.spanIndex)

        # Number of twist locations
        self.nTwist = self.nRefAxPts - 1 # Root is fixed

        # Some initializations which will be used later
        self.DV = []
        self.genSamples = 0

        # Overiding/setting some adflow solver options
        if self.options["aeroSolver"] == "adflow":
            self.options["solverOptions"]["printAllOptions"] = False
            self.options["solverOptions"]["printIntro"] = False
            self.options["solverOptions"]["outputDirectory"] = "."
            self.options["solverOptions"]["numberSolutions"] = False
            self.options["solverOptions"]["printTiming"] = False
            self.options["solverOptions"]["gridFile"] = self.options["gridFile"]
            self.options["solverOptions"]["liftIndex"] = self.options["liftIndex"]
            self.options["solverOptions"]["forcesAsTractions"] = False # bcoz of meld

    # ----------------------------------------------------------------------------
    #                       Design Variable related methods
    # ----------------------------------------------------------------------------

    def addDV(self, name: str, lowerBound, upperBound) -> None:
        """
            Method for adding a DV for analysis

            Parameters
            ----------
            name: str
                Name of the design variable, possible names: "shape", "twist", "mach", "alpha", "altitude"

            lowerBound: float or np.ndarray
                Lower bound of the design variable.

            upperBound: float or np.ndarray
                Upper bound of the design variable.
        """

        # Checking
        self._checkDV(name, lowerBound, upperBound)

        # Adding the pyGeo related DVs
        if name == "shape":
            self.DVGeo.addLocalDV("shape", lower=lowerBound, upper=upperBound, axis=self.liftIndex, scale=1.0)

        elif name == "twist":
            def twist_z(val, geo):
                for i in range(1, self.nRefAxPts):
                    geo.rot_z["wingAxis"].coef[i] = -val[i - 1]

            def twist_y(val, geo):
                for i in range(1, self.nRefAxPts):
                    geo.rot_y["wingAxis"].coef[i] = val[i - 1]

            # Adding the twist DV
            if self.spanIndex == "k":
                self.DVGeo.addGlobalDV(dvName="twist", value=[0]*self.nTwist, func=twist_z, lower=lowerBound, upper=upperBound)
            elif self.spanIndex == "j":
                self.DVGeo.addGlobalDV(dvName="twist", value=[0]*self.nTwist, func=twist_y, lower=lowerBound, upper=upperBound)

        # Converting to numpy array
        if isinstance(lowerBound, float):
            lowerBound = np.array([lowerBound])
            upperBound = np.array([upperBound])

        # Creating/Appending the design variable list
        if len(self.DV) == 0:
            self.upperBound = upperBound
            self.lowerBound = lowerBound
            self.locator = np.array(["{}".format(name)]*upperBound.shape[0])
        else:
            self.upperBound = np.append(self.upperBound, upperBound)
            self.lowerBound = np.append(self.lowerBound, lowerBound)
            self.locator = np.append(self.locator, np.array(["{}".format(name)]*upperBound.shape[0]))

        # Adding the DV to the list
        self.DV.append(name)

        # Limits for sampler
        xlimits = np.hstack((self.lowerBound.reshape(-1,1), self.upperBound.reshape(-1,1)))

        # Creating the sampler
        self.sampler = LHS(xlimits=xlimits, criterion=self.options["samplingCriterion"], random_state=self.options["randomState"])

    # ----------------------------------------------------------------------------
    #                   Methods related to sample generation
    # ----------------------------------------------------------------------------

    def generateSamples(self, numSamples: int=None, doe: np.ndarray=None) -> None:
        """
            This method generates samples using user provided number of 
            samples or user provided design of experiments (doe). It performs
            the analysis for each sample and saves the results in the output directory.

            Parameters
            ----------
            numSamples: int, optional
                Number of samples to be generated. If provided, then LHS samples are generated
                using `smt` package.

            doe: np.ndarray, optional
                User provided samples (doe) of size n x numSamples.

            **Note**: Either numSamples or doe should be provided, both cannot be provided at the same time.
        """

        # Checking if inputs to the method are correct
        if numSamples is None and doe is None:
            self._error("Provide either number of samples (int) or doe (2D numpy array)")

        if numSamples is not None and doe is not None:
            self._error("Provide either number of samples (int) or doe (2D numpy array)")

        # Checking if the appropriate options are set for analysis
        if self.options["solverOptions"] == {} or self.options["aeroProblem"] == None:
            self._error("You need to set solverOptions and aeroProblem in the options dictionary for running the analysis.")

        # Performing checks
        if len(self.DV) == 0:
            self._error("Add design variables before running the analysis")

        # Sampling plan
        if isinstance(numSamples, int):
            samples = self.sampler(numSamples)

        elif isinstance(doe, np.ndarray):
            if doe.ndim != 2:
                self._error("doe argument is not a 2D numpy array.")
            
            samples = doe
            numSamples = samples.shape[0]

        else:
            self._error("Data type of numSamples/doe is not correct.")

        # Number of analysis passed/failed
        failed =[]
        totalTime = 0

        # Creating empty dictionary for storing the data
        data = {}

        # Creating and writing a description file
        description = open("{}/description.txt".format(self.options["directory"]), "a", buffering=1)
        description.write("---------------------------------------------------")
        description.write("\Static Aerostructural data generation with FFD parametrization")
        description.write("\n--------------------------------------------------")
        description.write("\nDesign variables: {}".format(self.DV))
        description.write("\nNumber of DVs: {}".format(samples.shape[1]))
        description.write("\nTotal number of samples requested: {}".format(numSamples))
        description.write("\n-----------------------------")
        description.write("\nAnalysis specific description")
        description.write("\n-----------------------------")

        # Generate data
        for sampleNo in range(numSamples):

            description.write("\nAnalysis {}: ".format(sampleNo+1))

            # Current sample
            x = samples[sampleNo,:]

            # Starting time
            t1 = time.time()

            try:
                # Getting output for specific sample
                output = self.getObjectives(x)

            except Exception as e:
                print("Error occured during the analysis. Check analysis.log in the respective folder for more details.")
                failed.append(sampleNo + 1)
                description.write("\nAnalysis failed.")

            else:
                # Check for analysis failure
                if output["fail"] == True: # Check for analysis failure
                    failed.append(sampleNo + 1)
                    description.write("\nAnalysis failed.")

                # Creating a dictionary of data
                else:
                    if self.genSamples - len(failed) == 1:
                        data["x"] = np.array(x)
                        for value in output.keys():
                            data[value] = np.array([output[value]])

                    else:
                        # Appending data dictionary created earlier
                        data["x"] = np.vstack((data["x"], x))
                        for value in output.keys():
                            data[value] = np.vstack(( data[value], np.array([output[value]]) ))

                    # Saving the results
                    savemat("{}/data.mat".format(self.options["directory"]), data)

            finally:
                # Ending time
                t2 = time.time()

                totalTime += (t2-t1)/60

                # Writing time taken to file
                description.write("\nTime taken for analysis: {} min.".format((t2-t1)/60))

        # Making generated samples 0
        self.genSamples = 0

        # Writing final results in the description file
        description.write("\n--------------------------------------")
        description.write("\nTotal time taken for analysis: {} min.".format(totalTime))
        description.write("\nNumber of successful analysis: {}".format(numSamples - len(failed)))
        description.write("\nNumber of failed analysis: {}".format(len(failed)))
        if len(failed) != 0:
            description.write("\nFailed analysis: {}".format(failed))

        # Closing the description file
        description.close()

    def getObjectives(self, x: np.ndarray) -> dict:
        """
            Method for running a single analysis.

            Parameters
            ----------
            x: np.ndarray
                design variable for which the analysis is to be run.

            Returns
            -------
            output: dict
                A dictionary containing the outputs generated from the analysis.
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

        print("Running analysis {}".format(self.genSamples + 1))

        directory = self.options["directory"]

        # Create the folder for saving the results
        os.system("mkdir {}/{}".format(directory, self.genSamples+1))

        # Getting the directory where package is saved
        pkgdir = sys.modules["blackbox"].__path__[0]

        # Setting filepath based on the solver
        if self.options["solver"] == "adflow":
            filepath = os.path.join(pkgdir, "runscripts/runscript_aerostruct.py")

        # Copy the runscript to analysis directory
        shutil.copy(filepath, "{}/{}/runscript.py".format(directory, self.genSamples+1))

        # Creating the new design variable dict
        # If there are no shape DV, then DVGeo
        # will not update the wing surface.
        newDV = {}
        for dv in self.DV:
            loc = self.locator == dv
            loc = loc.reshape(-1,)
            newDV[dv] = x[loc]

        # Creating the new design variable dict
        self.DVGeo.setDesignVars(newDV)
        newSurfMesh = self.DVGeo.update("wing_surface_mesh")

        # Update the surface mesh in IdWarp
        self.mesh.setSurfaceCoordinates(newSurfMesh)

        # Deform the volume mesh
        self.mesh.warpMesh()

        # Changing the directory to analysis folder
        os.chdir("{}/{}".format(directory, self.genSamples+1))

        # Write deformed FFD file
        if self.options["writeDeformedFFD"]:
            self.DVGeo.writePlot3d("deformedFFD.xyz")

        # Write the new grid file
        self.mesh.writeGrid('volMesh.cgns')

        # Create input file
        self._createInputFile(x)





    # ----------------------------------------------------------------------------
    #          Other required methods
    # ----------------------------------------------------------------------------

    def _checkOptions(self, options: dict, requiredOptions: list) -> None:
        """
            This method validates user provided options.

            Parameters
            ----------
            options: dict
                User provided options.

            requiredOptions: list
                List of required options.
        """

        defaultOptions = list(self.options.keys())
        allowedUserOptions = defaultOptions
        allowedUserOptions.extend(requiredOptions)
        userProvidedOptions = list(options.keys())

        # Checking if user provided option contains only allowed attributes
        if not set(userProvidedOptions).issubset(allowedUserOptions):
            self._error("Option dictionary contains unrecognized attribute(s): {}"\
                        .format(set(userProvidedOptions) - set(allowedUserOptions)))

        # Checking if user has mentioned all the requried attributes
        if not set(requiredOptions).issubset(userProvidedOptions):
            self._error("Option dictionary doesn't contain following attribute(s): {}"\
                        .format(set(requiredOptions) - set(userProvidedOptions)))

        ############ Validating solverOptions
        if not isinstance(options["solverOptions"], dict):
            self._error("\"solverOptions\" attribute is not a dictionary.")

        ############ Validating gridFile
        if not os.path.exists(os.path.abspath(options["gridFile"])):
            self._error("Provided grid file doesn't exists.")

        ############ Validating structMeshFile
        if not os.path.exists(os.path.abspath(options["structMeshFile"])):
            self._error("Provided struct mesh file doesn't exists.")

        ############ Validating ffdFile
        if not os.path.exists(os.path.abspath(options["ffdFile"])):
            self._error("Provided FFD file doesn't exists.")

        ############ Validating liftIndex
        if not isinstance(options["liftIndex"], int):
            self._error("\"liftIndex\" attribute in options is not an integer.")

        if options["liftIndex"] not in [2, 3]:
            self._error("\"liftIndex\" attribute in options is not recognized. It can be either 2 or 3.")

        ############ Validating aeroProblem
        if not isinstance(options["aeroProblem"], AeroProblem):
            self._error("\"aeroProblem\" attribute is not an aeroproblem.")

        ############ Validating tacs problem setup
        if not callable(options["tacsProblemSetup"]):
            self._error("\"tacsProblemSetup\" attribute is not a callable.")

        ############ Validating tacs element callback
        if not callable(options["tacsElementCallback"]):
            self._error("\"tacsElementCallback\" attribute is not a callable.")

        ############ Validating aerosolver
        if "aeroSolver" in userProvidedOptions:
            if not isinstance(options["aeroSolver"], str):
                self._error("\"aeroSolver\" attribute is not a string.")

            if options["aeroSolver"] not in ["adflow", "dafoam"]:
                self._error("\"aeroSolver\" should be either adflow or dafoam")

            if options["aeroSolver"] == "dafoam":

                if "openfoamDir" in userProvidedOptions:
                    if not isinstance(options["openfoamDir"], str):
                        self._error("\"openfoamDir\" attribute is not a string.")

                    for fname in ["0", "constant", "system"]:
                        if not os.path.isdir(os.path.join(os.path.abspath(options["openfoamDir"]), f"{fname}")):
                            self._error(f"The folder \"{fname}\" requried for openfoam simulation is not available at given location")

                else:
                    for fname in ["0", "constant", "system"]:
                        if not os.path.isdir(os.path.join(os.path.abspath(self.options["openfoamDir"]), f"{fname}")):
                            self._error(f"The folder \"{fname}\" requried for openfoam simulation is not available at given location")

        ############ Validating noOfProcessors
        if "noOfProcessors" in userProvidedOptions:
            if not isinstance(options["noOfProcessors"], int):
                self._error("\"noOfProcessors\" attribute is not an integer.")

            if psutil.cpu_count(False) < options["noOfProcessors"] + 1:
                self._error("\"noOfProcessors\" requested is more than available processors.")

        ############ Validating sliceLocation
        if "sliceLocation" in userProvidedOptions:
            if not isinstance(options["sliceLocation"], list):
                self._error("\"sliceLocation\" attribute is not a list of relative slice locations on wing.")

        ############ Validating writeDeformedFFD
        if "writeDeformedFFD" in userProvidedOptions:
            if not isinstance(options["writeDeformedFFD"], bool):
                self._error("\"writeDeformedFFD\" attribute is not a boolean value.")

        ############ Validating storeFieldData
        if "storeFieldData" in userProvidedOptions:
            if not isinstance(options["storeFieldData"], bool):
                self._error("\"storeFieldData\" attribute is not a boolean value.")

        ############ Validating storeMDAHistory
        if "storeMDAHistory" in userProvidedOptions:
            if not isinstance(options["storeMDAHistory"], bool):
                self._error("\"storeMDAHistory\" attribute is not a boolean value.")

        ############ Validating directory attribute
        if "directory" in userProvidedOptions:
            if not isinstance(options["directory"], str):
                self._error("\"directory\" attribute is not string.")

    def _checkDV(self, name: str, lb, ub) -> None:
        """
            Method for validating a DV user wants to add

            Parameters
            ----------
            name: str
                Name of the design variable

            lb: float or np.ndarray
                Lower bound of the design variable

            ub: float or np.ndarray
                Upper bound of the design variable
        """

        # List of possible DVs
        possibleDVs = ["shape", "twist", "alpha", "mach", "altitude"]

        # Validating name of the DV
        if not isinstance(name, str):
            self._error("Name argument is not a string.")

        # Checking if the DV is allowed
        if name not in possibleDVs:
            self._error("{} argument is not a valid DV.".format(name))

        # Checking if the DV is already added
        if name in self.DV:
            self._error("{} already exists.".format(name))

        # Checking if these variables are initialized through aero problem or not
        if name in ["mach", "altitude"]:
            if name not in self.options["aeroProblem"].inputs.keys():
                self._error("You need to initialize \"{}\" in the aero problem to set it as design variable.".format(name))

        # Validating the bounds for "shape" variable
        if name == "shape" or name == "twist":
            if not isinstance(lb, np.ndarray) or lb.ndim != 1:
                self._error("Lower bound for \"shape\" variable should be a 1D numpy array.")

            if not isinstance(ub, np.ndarray) or ub.ndim != 1:
                self._error("Upper bound for \"shape\" variable should be a 1D numpy array.")
            
            if name == "shape":
                if len(lb) != self.nffd:
                    self._error("Length of lower bound array is not equal to number of FFD points.")

                if len(ub) != self.nffd:
                    self._error("Length of upper bound array is not equal to number of FFD points.")

            elif name == "twist":
                if len(lb) != self.nTwist:
                    self._error("Length of lower bound array is not equal to number of twist variables.")

                if len(ub) != self.nTwist:
                    self._error("Length of upper bound array is not equal to number of twist variables.")

            if np.any(lb >= ub):
                self._error("Lower bound is greater than or equal to upper bound for atleast one DV.")

        else:
            if not isinstance(lb, float):
                self._error("Lower Bound argument is not a float.")

            if not isinstance(ub, float):
                self._error("Upper Bound argument is not a float.")

            if lb >= ub:
                self._error("Lower bound is greater than or equal to upper bound.")
    
    def _getDefaultOptions(self):
        """
            Setting up the initial values of options which are common across all functions.
        """
        
        defaultOptions = DefaultOptions()

        for key in vars(defaultOptions):
            value = getattr(defaultOptions, key)
            self.options[key] = value

    def _setOptions(self, options):
        """
            Method for assigning user provided options.
        """

        for key in options.keys():
            if isinstance(options[key], dict):
                if key in self.options.keys():
                    # If the value is dictionary, update the default dictionary.
                    # Otherwise, assign values.
                    self.options[key].update(options[key]) 
                else:
                    self.options[key] = options[key]
            else:
                self.options[key] = options[key]

    def _createInputFile(self):
        """
            Method to create an input file for analysis
        """

        directory = self.options["directory"]

        input = {}
        sample = {}
        parameters = self.options["fixedParameters"]
        objectives = self.options["objectives"]

        if self.options["type"] == "multi":
            for sampleNo in range(self.options["numberOfSamples"]):
                os.chdir("{}/{}".format(directory,sampleNo))

                for key in self.samples:
                    sample[key] = self.samples[key][sampleNo]

                input = {
                    "aeroSolverOptions" : self.options["aeroSolverOptions"],
                    "sample" : sample,
                    "parameters" : parameters,
                    "objectives" : objectives
                }

                filehandler = open("input.pickle", "xb")
                pickle.dump(input, filehandler)
                filehandler.close()
                os.chdir("../..")

        elif self.options["type"] == "single":
            os.chdir("{}/{}".format(directory, self.sampleNo))

            input = {
                "aeroSolverOptions" : self.options["aeroSolverOptions"],
                "sample" : self.samples,
                "parameters" : parameters,
                "objectives" : objectives
            }

            filehandler = open("input.pickle", "xb")
            pickle.dump(input, filehandler)
            filehandler.close()

            os.chdir("../..")

    def _error(self, message: str) -> None:
        """
            Method for printing errors in nice format

            Parameters
            ----------
            message: str
                Error message to be printed
        """

        # Initial message - total len is 80 characters
        msg = "\n+" + "-" * 78 + "+" + "\n" + "| Blackbox Error: "

        # Initial number of characters
        i = 16

        for word in message.split():
            if len(word) + i + 1 > 76:  # Finish line and start new one
                msg += " " * (76 - i) + " |\n| " + word + " " # Adding space and word in new line
                i = len(word) + 1 # Setting i value for new line
            else:
                msg += word + " " # Adding the word with a space
                i += len(word) + 1 # Increase the number of characters
        msg += " " * (76 - i) + " |\n" + "+" + "-" * 78 + "+" + "\n" # Adding last line
 
        print(msg, flush=True)

        exit()
