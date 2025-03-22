# This script contains parent classes for different airfoil use cases in blackbox

# General imports
import os, pickle, psutil, time, shutil, sys
import numpy as np
from baseclasses import AeroProblem
from scipy.io import savemat
from mpi4py import MPI

comm = MPI.COMM_WORLD

class DefaultOptions():
    """
        Class creates a default option list which is later 
        edited/appended with user provided options. This class
        can be used as a parent class for application specific
        default class. No need to define instance variables 
        in child which are already defined here.
    """

    def __init__(self):

        # Aero solver / Meshing Options
        self.solver = "adflow"
        self.openfoamDir = "."
        self.solverOptions = {}
        self.meshingOptions = {}
        self.aeroProblem = None

        # Other options
        self.directory = "output"
        self.noOfProcessors = 4
        self.refine = 0
        self.writeSliceFile = False
        self.writeAirfoilCoordinates = False
        self.plotAirfoil = False
        self.writeDeformedFFD = False

        # Flow-field related options
        # if getFlowFieldData is false, then all other options are useless
        self.getFlowFieldData = False
        self.region = "surface"

        # Alpha implicit related options
        self.alpha = "explicit"
        self.targetCL = 0.824
        self.targetCLTol = 1e-4
        self.startingAlpha = 2.5

        # Sampling options
        self.samplingCriterion = "cm"
        self.randomState = None

        # FFD Box related options
        self.fitted = False
        self.xmargin = 0.001
        self.ymarginu = 0.02
        self.ymarginl = 0.02

        # LE/TE related options
        self.fixLETE = True

        # Smoothing options
        self.smoothing = False
        self.smoothingTheta = 0.75
        self.smoothingMaxIterations = 100
        self.smoothingTolerance = 5e-4

class AirfoilBaseClass():
    """
        Base class for airfoil related classes.

        This class needs to be inherited by all the airfoil related classes.

        Wherever this class is used, child class needs to initialize various
        variables for proper working of the class.
    """

    # ----------------------------------------------------------------------------
    #                       Methods related to analysis
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
        if self.options["solverOptions"] == {} or self.options["meshingOptions"] == {} or self.options["aeroProblem"] == None:
            self._error("You need to set solverOptions, meshingOptions and aeroProblem in the options dictionary for running the analysis.")

        # Performing checks
        if len(self.DV) == 0:
            self._error("Add design variables before running the analysis")
        
        # Sampling plan
        if isinstance(numSamples, int):
            samples = self.sampler(numSamples)
        elif isinstance(doe, np.ndarray):
            samples = doe
            numSamples = samples.shape[0]
        else:
            self._error("Data type of numSamples/doe is not correct.")

        # Number of analysis passed/failed
        failed = []
        totalTime = 0

        # Creating empty dictionary for storing the data
        data = {}
        fieldData = {}

        # Creating and writing a description file
        description = open("{}/description.txt".format(self.options["directory"]), "a", buffering=1)
        description.write("---------------------------------------------------")
        description.write("\nAirfoil sample generation with {} parametrization".format(self.parametrization))
        description.write("\n--------------------------------------------------")
        description.write("\nDesign variables: {}".format(self.DV))
        description.write("\nLower bound for design variables:\n{}".format(self.lowerBound))
        description.write("\nUpper bound for design variables:\n{}".format(self.upperBound))
        description.write("\nTotal number of samples requested: {}".format(numSamples))
        description.write("\n-----------------------------")
        description.write("\nAnalysis specific description")
        description.write("\n-----------------------------")

        # Generate data
        for sampleNo in range(numSamples):

            description.write("\nAnalysis {}: ".format(sampleNo+1))

            # Current sample
            x = samples[sampleNo,:]

            # Laplacian smoothing
            if self.parametrization == "FFD" and self.options["smoothing"]:
                x = self.LaplacianSmoothing(x)

            description.write("\nDesign Variable: {}".format(x))

            # Starting time
            t1 = time.time()

            try:
                # Getting output for specific sample
                output, field = self.getObjectives(x)

            except:
                print("Error occured during the analysis. Check analysis.log in the respective folder for more details.")
                failed.append(sampleNo + 1)
                description.write("\nAnalysis failed.")

            else:
                # Check for analysis failure
                if output["fail"] == True: # Check for analysis failure
                    failed.append(sampleNo + 1)
                    description.write("\nAnalysis failed due solver not converging or due to some other reason.")

                # Creating a dictionary of data
                else:
                    if self.genSamples - len(failed) == 1:
                        data["x"] = np.array(x)
                        for value in output.keys():
                            data[value] = np.array([output[value]])

                        # Creating a dictionary of field data
                        if self.options["getFlowFieldData"]:
                            fieldData["x"] = np.array(x)
                            for value in field.keys():
                                if field[value].ndim == 2:
                                    fieldData[value] = field[value].reshape(1,-1,3)
                                else:
                                    fieldData[value] = field[value].reshape(1,-1)

                    else:
                        # Appending data dictionary created earlier
                        data["x"] = np.vstack((data["x"], x))
                        for value in output.keys():
                            data[value] = np.vstack(( data[value], np.array([output[value]]) ))

                        # Appending field data dictionary created earlier
                        if self.options["getFlowFieldData"]:
                            fieldData["x"] = np.vstack((fieldData["x"], x))
                            for value in field.keys():
                                if field[value].ndim == 2:
                                    fieldData[value] = np.vstack(( fieldData[value], field[value].reshape(1,-1,3) ))
                                else:
                                    fieldData[value] = np.vstack(( fieldData[value], field[value].reshape(1,-1) ))

                    # Saving the results
                    savemat("{}/data.mat".format(self.options["directory"]), data)

                    # Saving the field data
                    if self.options["getFlowFieldData"]:
                        savemat("{}/fieldData.mat".format(self.options["directory"]), fieldData)

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

    def getObjectives(self, x: np.ndarray) -> tuple:
        """
            Method for running a single analysis for a given sample.

            Parameters
            ----------

            x: 1D numpy array
                Design variable for the analysis.

            Returns
            -------
            output: tuple
                A tuple containing the outputs generated from the analysis.
                First entry of the tuple is a dictionary containing scalar outputs
                from the analysis. Second entry is a dictionary containing the
                field data, if requested. Otherwise, it is None.
        """

        # Checking if the appropriate options are set for analysis
        if self.options["solverOptions"] == {} or self.options["meshingOptions"] == {} or self.options["aeroProblem"] == None:
            self._error("You need to set solverOptions, meshingOptions and aeroProblem in the options dictionary for running the analysis.")

        # Overiding/set some solver options
        if self.options["solver"] == "adflow":
            self.options["solverOptions"]["printAllOptions"] = False
            self.options["solverOptions"]["printIntro"] = False
            self.options["solverOptions"]["outputDirectory"] = "."
            self.options["solverOptions"]["numberSolutions"] = False
            self.options["solverOptions"]["printTiming"] = False

        # Raise an error if pyvista is not installed
        if self.options["getFlowFieldData"]:
            self.options["solverOptions"]["writeSurfaceSolution"] = True

        # Getting the deformed airfoil
        points = self.getAirfoil(x)

        print("Running analysis {}".format(self.genSamples + 1))

        directory = self.options["directory"]

        # Create the folder for saving the results
        os.system("mkdir {}/{}".format(directory, self.genSamples+1))

        # Getting the directory where package is saved
        pkgdir = sys.modules["blackbox"].__path__[0]

        # Setting filepath based on the how alpha is treated
        if self.options["solver"] == "adflow":

            if type(self).__name__ == "AirfoilCSTMultipoint":
                filepath = os.path.join(pkgdir, "runscripts/airfoil/runscript_airfoil_cst_mp.py")
            else:
                if self.options["alpha"] == "explicit":
                    filepath = os.path.join(pkgdir, "runscripts/airfoil/runscript_airfoil.py")
                else:
                    # filepath = os.path.join(pkgdir, "runscripts/airfoil/runscipt_airfoil_cst_opt.py")
                    filepath = os.path.join(pkgdir, "runscripts/airfoil/runscript_airfoil_rf.py")

        elif self.options["solver"] == "dafoam":

            filepath = os.path.join(pkgdir, "runscripts/airfoil/runscript_airfoil_dafoam.py")

        # Copy the runscript to analysis directory
        shutil.copy(filepath, "{}/{}/runscript.py".format(directory, self.genSamples+1))

        # Copying openfoam standard folders
        if self.options["solver"] == "dafoam":
            ofdir = self.options["openfoamDir"]
            os.system(f"cp -r {ofdir}/0 {directory}/{self.genSamples+1}")
            os.system(f"cp -r {ofdir}/constant {directory}/{self.genSamples+1}")
            os.system(f"cp -r {ofdir}/system {directory}/{self.genSamples+1}")

        # Changing the directory to analysis folder
        os.chdir("{}/{}".format(directory, self.genSamples+1))

        if self.options["writeAirfoilCoordinates"]:
            self._writeCoords(coords=points, filename="deformedAirfoil.dat")

        if self.options["plotAirfoil"]:
            self._plotAirfoil(self.plt, self.origCoords, points)

        if self.parametrization == "FFD" and self.options["writeDeformedFFD"]:
            self.DVGeo.writePlot3d("deformedFFD.xyz")

        # Create input file
        self._creatInputFile(x)

        # Writing the surface mesh
        self._writeSurfMesh(coords=points, filename="surfMesh.xyz")

        try:
            # Spawning the runscript on desired number of processors
            child_comm = MPI.COMM_SELF.Spawn(sys.executable, args=["runscript.py"], maxprocs=self.options["noOfProcessors"])

            # Creating empty process id list
            pid_list = []

            # Getting each spawned process
            for processor in range(self.options["noOfProcessors"]):
                pid = child_comm.recv(source=MPI.ANY_SOURCE, tag=processor)
                pid_list.append(psutil.Process(pid))

            # Disconnecting from intercommunicator
            child_comm.Disconnect()

            # Waiting till all the child processors are finished
            while len(pid_list) != 0:
                for pid in pid_list:
                    if not pid.is_running():
                        pid_list.remove(pid)

            # Reading the output file containing results
            filehandler = open("output.pickle", 'rb')

        except:
            raise Exception

        else:
            # Read the output
            output = pickle.load(filehandler)
            filehandler.close()

            # Calculate the area
            x = points[:,0]
            y = points[:,1]

            area = 0.0
            N = len(x)
            j = N - 1
            for i in range(0,N):
                area += (x[j] + x[i]) * (y[j] - y[i])
                j = i
            output["area"] = abs(area)/2.0

            if self.options["getFlowFieldData"]:
                # Reading the cgns file
                filename = self.options["aeroProblem"].name + "_surf.cgns"
                reader = self.pyvista.CGNSReader(filename)
                reader.load_boundary_patch = False

                # Reading the mesh
                mesh = reader.read()

                # Setting region for extraction
                if self.options["region"] == "surface":
                    mesh = mesh[0][0]
                else:
                    mesh = mesh[0][2]

                # Get the values
                fieldData = {}

                for var in mesh.array_names:
                    # Skipping the first entry in the array
                    if var != "Base/Zone":
                        # set_active_scalars returns a tuple, and second
                        # entry contains the pyvista numpy array.
                        fieldData[var] = np.asarray(mesh.set_active_scalars(var, "cell")[1])

            else:
                fieldData = None

            return output, fieldData

        finally:
            # Cleaning the directory
            files = ["volMesh.cgns", "input.pickle", "runscript.py",
                    "output.pickle", "fort.6", "opt.hst", "surfMesh.xyz"] 
            for file in files:
                if os.path.exists(file):
                    os.system("rm {}".format(file))

            # Changing the directory back to root
            os.chdir("../..")

            # Increase the number of generated samples
            self.genSamples += 1

    def getAirfoil(self, x: np.ndarray) -> np.ndarray:
        """
            Method for getting the airfoil coordinates from a design variable
            using parameterization within pyGeo.

            Parameters
            ----------
            x: 1D numpy array
                design variable

            Returns
            -------
            points: 2D numpy array
                Airfoil coordinates based on the design variable.
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
            return self.origCoords[:,0:2]

        # Creating dictionary from x
        newDV = {}
        for dv in self.DV:
            loc = self.locator == dv
            loc = loc.reshape(-1,)
            newDV[dv] = x[loc]

        if type(self).__name__ == "AirfoilCSTMultipoint":
            for ap in self.options["aeroProblem"]:
                ap.setDesignVars(newDV)

        if self.parametrization == "FFD" and self.options["fixLETE"]:

            # Adjusting LE FFD points
            midpoint = newDV["shape"][0]/2
            newDV["shape"][0] -= midpoint
            newDV["shape"] = np.append(-midpoint, newDV["shape"])

            # Adjusting TE FFD points
            midpoint = newDV["shape"][-1]/2
            newDV["shape"][-1] -= midpoint
            newDV["shape"] = np.append(newDV["shape"], -midpoint)

        # Updating the airfoil pointset based on new DV
        self.DVGeo.setDesignVars(newDV)

        # Getting the updated airfoil points
        points = self.DVGeo.update("airfoil")[:,0:2]

        return points

    def calculateArea(self, x: np.ndarray) -> float:
        """
            Note: This function should not be called in the middle of analysis
            It should ONLY be used from outside. Do not use this method within 
            getObjectives. That method has its own implementation
            of area calculation.

            Method to calculate the area of the airfoil
            based on the value of design variable.

            Parameters
            ----------
            x: 1D numpy array
                design variable

            Returns
            -------
            area: float
                Area of the airfoil.

            **Note**: To use this method, atleast lower or upper surface CST
            coefficient should be added as a DV.
        """

        # Getting the updated airfoil points
        points = self.getAirfoil(x)
        x = points[:,0]
        y = points[:,1]

        area = 0.0
        N = len(x)
        j = N - 1
        for i in range(0,N):
            area += (x[j] + x[i]) * (y[j] - y[i])
            j = i
        area = abs(area)/2.0

        return area
    
    # ----------------------------------------------------------------------------
    #                       Methods related to validation
    # ----------------------------------------------------------------------------

    def _checkOptions(self, defaultOptions: list, requiredOptions: list, options: dict) -> None:
        """
            This is a general method for validating options dictionary provided by
            user based on default and mandatory options.

            Parameters
            ----------
            defaultOptions: list
                List of default options.

            requiredOptions: list
                List of required options.

            options: dict
                User provided options.

            **Note**: This method is for internal use only
        """

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

        ############ Validating airfoilFile
        if "airfoilFile" in userProvidedOptions:
            if not isinstance(options["airfoilFile"], str):
                self._error("\"airfoilFile\" attribute is not string.")

            if not os.path.exists(os.path.abspath(options["airfoilFile"])):
                self._error("\"airfoilFile\" doesn't exists.")
            else:
                options["airfoilFile"] = os.path.abspath(options["airfoilFile"])

        ############ Validating nffd
        if "nffd" in userProvidedOptions:
            if not isinstance(options["nffd"], int):
                self._error("\"nffd\" attribute is not an integer.")

        ############ Validating numCST
        if "numCST" in userProvidedOptions:
            if not isinstance(options["numCST"], list):
                self._error("\"numCST\" is not a list")
            else:
                if len(options["numCST"]) != 2:
                    self._error("\"numCST\" should have only two entries.")

                if not isinstance(options["numCST"][0], int):
                    self._error("First entry in \"numCST\" is not an integer.")
                elif options["numCST"][0] <= 0:
                    self._error("First entry in \"numCST\" is less than 1.")
                        
                if not isinstance(options["numCST"][1], int):
                    self._error("Second entry in \"numCST\" is not an integer.")
                elif options["numCST"][1] <= 0:
                    self._error("Second entry in \"numCST\" is less than 1.")

        ############ Validating openfoamDir option
        if "openfoamDir" in userProvidedOptions:
            if not isinstance(options["openfoamDir"], str):
                self._error("\"openfoamDir\" attribute is not a string.")

            for fname in ["0", "constant", "system"]:
                if not os.path.isdir(os.path.join(os.path.abspath(options["openfoamDir"]), f"{fname}")):
                    self._error(f"The folder \"{fname}\" requried for openfoam simulation is not available at given location")

        ############ Validating solver
        if "solver" in userProvidedOptions:
            if not isinstance(options["solver"], str):
                self._error("\"solver\" attribute is not a string.")

            if options["solver"] not in ["adflow", "dafoam"]:
                self._error("\"solver\" should be either adflow or dafoam")

            if options["solver"] == "dafoam":

                if "openfoamDir" in userProvidedOptions:
                    if not isinstance(options["openfoamDir"], str):
                        self._error("\"openfoamDir\" attribute is not a string.")

                    for fname in ["0", "constant", "system"]:
                        if not os.path.isdir(os.path.join(os.path.abspath(options["openfoamDir"]), f"{fname}")):
                            self._error(f"The folder \"{fname}\" requried for openfoam simulation is not available at given location")

                else:
                    for fname in ["0", "constant", "system"]:
                        if not os.path.isdir(os.path.join(os.path.abspath("."), f"{fname}")):
                            self._error(f"The folder \"{fname}\" requried for openfoam simulation is not available at given location")

        ############ Validating aeroProblem
        if type(self).__name__ == "AirfoilCSTMultipoint":
            if "aeroProblem" in userProvidedOptions:
                if not isinstance(options["aeroProblem"], list):
                    self._error("\"aeroProblem\" attribute is not a list.")

                for ap in options["aeroProblem"]:
                    if not isinstance(ap, AeroProblem):
                        self._error("\"aeroProblem\" attribute doesn't contain AeroProblem object.")
        else:
            if "aeroProblem" in userProvidedOptions:
                if not isinstance(options["aeroProblem"], AeroProblem):
                    self._error("\"aeroProblem\" attribute is not an aeroproblem.")

        ############ Validating solverOptions
        if "solverOptions" in userProvidedOptions:
            if not isinstance(options["solverOptions"], dict):
                self._error("\"solverOptions\" attribute is not a dictionary.")

            # if "gridFile" in options["solverOptions"].keys():
            #     self._error("\"gridFile\" attribute in solver options is not required.", type=1)

        ############ Validating meshingOptions
        if "meshingOptions" in userProvidedOptions:
            if not isinstance(options["meshingOptions"], dict):
                self._error("\"meshingOptions\" attribute is not a dictionary.")

            if "inputFile" in options["meshingOptions"].keys():
                self._error("\"inputFile\" attribute in meshing options is not required.", type=1)

        ############ Validating noOfProcessors
        if "noOfProcessors" in userProvidedOptions:
            if not isinstance(options["noOfProcessors"], int):
                self._error("\"noOfProcessors\" attribute is not an integer.")

            if psutil.cpu_count(False) < options["noOfProcessors"] + 1:
                self._error("\"noOfProcessors\" requested is more than available processors.")

        ############ Validating refine
        if "refine" in userProvidedOptions:
            if not isinstance(options["refine"], int):
                self._error("\"refine\" attribute is not an integer.")

        ############ Validating alpha
        if "alpha" in userProvidedOptions:
            if not isinstance(options["alpha"], str):
                self._error("\"alpha\" attribute is not string.")

            if options["alpha"] not in ["explicit", "implicit"]:
                self._error("\"alpha\" attribute is not recognized. It can be either \"explicit\" or \"implicit\".")

            if options["alpha"] == "implicit":
                if "targetCL" in userProvidedOptions:
                    if not isinstance(options["targetCL"], float):
                        self._error("\"targetCL\" option is not float.")

            if "targetCL" in userProvidedOptions:
                if not isinstance(options["targetCL"], float):
                    self._error("\"targetCL\" option is not float.")

            if "targetCLTol" in userProvidedOptions:
                    if not isinstance(options["targetCLTol"], float):
                        self._error("\"targetCLTol\" option is not float.")
                
            if "startingAlpha" in userProvidedOptions:
                if not isinstance(options["startingAlpha"], float):
                    self._error("\"startingAlpha\" option is not float.")

        ############ Validating writeSliceFile
        if "writeSliceFile" in userProvidedOptions:
            if not isinstance(options["writeSliceFile"], bool):
                self._error("\"writeSliceFile\" attribute is not a boolean value.")

        ############ Validating directory attribute
        if "directory" in userProvidedOptions:
            if not isinstance(options["directory"], str):
                self._error("\"directory\" attribute is not string.")

        ############ Validating writeAirfoilCoordinates attribute
        if "writeAirfoilCoordinates" in userProvidedOptions:
            if not isinstance(options["writeAirfoilCoordinates"], bool):
                self._error("\"writeAirfoilCoordinates\" attribute is not a boolean value.")

        ############ Validating plotAirfoil attribute
        if "plotAirfoil" in userProvidedOptions:
            if not isinstance(options["plotAirfoil"], bool):
                self._error("\"plotAirfoil\" attribute is not a boolean value.")

        ############ Validating getFlowFieldData attribute
        if "getFlowFieldData" in userProvidedOptions:
            if not isinstance(options["getFlowFieldData"], bool):
                self._error("\"getFlowFieldData\" attribute is not a boolean value.")

            # Checking the other related options
            if "region" in userProvidedOptions:
                if not isinstance(options["region"], str):
                    self._error("\"region\" attribute is not a string.")

                if options["region"] not in ["surface", "field"]:
                    self._error("\"region\" attribute is not recognized. It can be either \"surface\" or \"field\".")
        
        ############ Validating FFD options
        if "fitted" in userProvidedOptions:
            if not isinstance(options["fitted"], bool):
                self._error("\"fitted\" attribute is not a boolean value.")

        if "xmargin" in userProvidedOptions:
            if not isinstance(options["xmargin"], float):
                self._error("\"xmargin\" attribute is not a float.")

        if "ymarginu" in userProvidedOptions:
            if not isinstance(options["ymarginu"], float):
                self._error("\"ymarginu\" attribute is not a float.")

        if "ymarginl" in userProvidedOptions:
            if not isinstance(options["ymarginl"], float):
                self._error("\"ymarginl\" attribute is not a float.")

        if "writeDeformedFFD" in userProvidedOptions:
            if not isinstance(options["writeDeformedFFD"], bool):
                self._error("\"writeDeformedFFD\" attribute is not a boolean.")

        ############ Validating sampling options
        if "sampling" in userProvidedOptions:
            if not isinstance(options["sampling"], str):
                self._error("\"sampling\" attribute is not a string.")

            if options["sampling"] not in ["internal", "external"]:
                self._error("\"sampling\" attribute is not recognized. It can be either \"internal\" or \"external\".")

        if "samplingCriterion" in userProvidedOptions:
            if not isinstance(options["samplingCriterion"], str):
                self._error("\"samplingCriterion\" attribute is not a string.")

            if options["samplingCriterion"] not in ["c", "m", "cm", "ese"]:
                self._error("\"samplingCriterion\" attribute is not recognized. It can be \"c\", \"m\", \"cm\", or \"ese\".")

        if "randomState" in userProvidedOptions:
            if not isinstance(options["randomState"], int):
                self._error("\"randomState\" attribute is not an integer.")

        ############ Validating LE/TE options
        if "fixLETE" in userProvidedOptions:
            if not isinstance(options["fixLETE"], bool):
                self._error("\"fixLETE\" attribute is not a boolean.")

        ############ Validating smoothing options
        if "smoothing" in userProvidedOptions:
            if not isinstance(options["smoothing"], bool):
                self._error("\"smoothing\" attribute is not a boolean.")

        if "smoothingTheta" in userProvidedOptions:
            if not isinstance(options["smoothingTheta"], float):
                self._error("\"smoothingTheta\" attribute is not a float.")

        if "smoothingMaxIterations" in userProvidedOptions:
            if not isinstance(options["smoothingMaxIterations"], int):
                self._error("\"smoothingMaxIterations\" attribute is not an integer.")

        if "smoothingTolerance" in userProvidedOptions:
            if not isinstance(options["smoothingTolerance"], float):
                self._error("\"smoothingTolerance\" attribute is not a float.")

    # ----------------------------------------------------------------------------
    #                               Other methods
    # ----------------------------------------------------------------------------

    def _getDefaultOptions(self, defaultOptions) -> None:
        """
            Setting up the initial values of options.

            Input is default options class.
        """

        for key in vars(defaultOptions):
            value = getattr(defaultOptions, key)
            self.options[key] = value

    def _setOptions(self, options: dict) -> None:
        """
            Method for assigning user provided options.
        """

        for key in options.keys():
            # update strategy is different for a dictionary
            # and other type of variables
            if isinstance(options[key], dict):
                # if the option is already present in default dictionary, 
                # then append the user provided key-value pairs to the 
                # default dictionary. For example: solverOptions.
                if key in self.options.keys():
                    self.options[key].update(options[key]) 
                else:
                    self.options[key] = options[key]
            else:
                self.options[key] = options[key]

    def _creatInputFile(self, x:np.ndarray) -> None:
        """
            Method to create an input file for a specific analysis.

            Parameters
            ----------
            x: 1D numpy array
                design variable

            **Note**: This method is for internal use only
        """

        # Creating input dict
        input = {
            "solverOptions": self.options["solverOptions"],
            "aeroProblem": self.options["aeroProblem"],
            "meshingOptions": self.options["meshingOptions"],
            "refine": self.options["refine"],
            "writeSliceFile": self.options["writeSliceFile"]
        }

        # Adding non-shape DV
        if "alpha" in self.DV:
            loc = self.locator == "alpha"
            loc = loc.reshape(-1,)
            input["alpha"] = x[loc]

        if "mach" in self.DV:
            loc = self.locator == "mach"
            loc = loc.reshape(-1,)
            input["mach"] = x[loc]

        if "altitude" in self.DV:
            loc = self.locator == "altitude"
            loc = loc.reshape(-1,)
            input["altitude"] = x[loc]

        # Adding target Cl if alpha is implicit
        if self.options["alpha"] == "implicit":
            input["targetCL"] = self.options["targetCL"]
            input["targetCLTol"] = self.options["targetCLTol"]
            input["startingAlpha"] = self.options["startingAlpha"]

        # Saving the input file
        filehandler = open("input.pickle", "xb")
        pickle.dump(input, filehandler)
        filehandler.close()

    def _writeCoords(self, coords: np.ndarray, filename: str) -> None:
        """
            Writes out a set of airfoil coordinates in selig format in dat format

            Parameters
            ----------
            coords: 2D numpy array
                Airfoil coordinates.

            filename: str
                Name of the file to write the coordinates.

            **Note**: This method is for internal use only
        """

        # X and Y ccordinates of the airfoil
        x = coords[:, 0]
        y = coords[:, 1]

        with open(filename, "w") as f:
            for i in range(len(x)):
                f.write(str(round(x[i], 12)) + "\t\t" + str(round(y[i], 12)) + "\n")

        f.close()

    def _writeSurfMesh(self, coords, filename, zHeight=1.0):
        """
            Writes out surface mesh in Plot 3D format (only one element in z direction)

            Parameters
            ----------
            coords: 2D numpy array
                Airfoil coordinates.

            filename: str
                Name of the file to write the coordinates.

            zHeigth: float
                Distance to travel in z direction, default=1.0
            
            **Note**: This method is for internal use only
        """

        # X and Y ccordinates of the airfoil
        x = coords[:, 0]
        y = coords[:, 1]
        z = [0.0, zHeight]

        # Writing the file
        with open(filename, "w") as f:
            f.write("1\n")
            f.write("%d %d %d\n" % (len(x), 2, 1))
            for iDim in range(3):
                for j in range(2):
                    for i in range(len(x)):
                        if iDim == 0:
                            f.write("%g\n" % x[i])
                        elif iDim == 1:
                            f.write("%g\n" % y[i])
                        else:
                            # f.write("%g\n" % (float(j)))
                            f.write("%g\n" % z[j])

        f.close()

    def _plotAirfoil(self, plt, orig_airfoil: np.ndarray, def_airfoil: np.ndarray) -> None:
        """
            Method for plotting the base airfoil
            and the deformed airfoil.

            Parameters
            ----------
            plt: matplotlib.pyplot
                Matplotlib pyplot object.
            
            orig_airfoil: 2D numpy array
                Original airfoil coordinates.

            def_airfoil: 2D numpy array
                Deformed airfoil coordinates.
        """

        _, ax = plt.subplots()

        ax.plot(orig_airfoil[:,0], orig_airfoil[:,1], label="Original airfoil")
        ax.plot(def_airfoil[:,0], def_airfoil[:,1], label="Deformed airfoil")
        ax.set_xlabel("x/c", fontsize=14)
        ax.set_ylabel("y/c", fontsize=14)
        ax.legend(fontsize=12)

        plt.savefig("airfoil.png", dpi=400)

        plt.close()

    def _error(self, message: str, type=0) -> None:
        """
            Method for printing errors in nice manner.

            Parameters
            ----------
            message: str
                Message to be displayed.

            type: int
                Type of the message. 0 for error and 1 for warning.
                Default is 0.
        """

        # Initial message - total len is 80 characters
        msg = "\n+" + "-" * 78 + "+" + "\n" + "| Blackbox {}: ".format("Error" if type == 0 else "Warning")

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

        if type == 0:
            exit()
