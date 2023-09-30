# This script contains parent classes for different airfoil use cases in blackbox

# General imports
import os, pickle, psutil, time
import numpy as np
from pyDOE2 import lhs
from baseclasses import AeroProblem
from scipy.io import savemat
from scipy import integrate

class DefaultOptions():
    """
        Class creates a default option list which is later 
        edited/appended with user provided options. This class
        can be used as a parent class for application specific
        default class. No need to define instance variables 
        in child which are already defined here.
    """

    def __init__(self):

        # Aero solver Options
        self.aeroSolver = "adflow"

        # Other options
        self.directory = "output"
        self.noOfProcessors = 4
        self.refine = 0
        self.writeSliceFile = False
        self.writeAirfoilCoordinates = False
        self.plotAirfoil = False

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
        self.sampling = "internal"
        self.samplingCriterion = "cm"
        self.randomState = None

class AirfoilBaseClass():
    """
        Base class for airfoil related classes.

        This class needs to be inherited by all the airfoil related classes.
    """

    # ----------------------------------------------------------------------------
    #                       Methods related to analysis
    # ----------------------------------------------------------------------------

    def generateSamples(self, numSamples: int=None, doe: np.ndarray=None) -> None:
        """
            Method for generating samples.

            Inputs (only one is rquired, depending on user preference):
            numSamples (int): Number of samples
            doe (np.ndarray): User provided samples of size n x numSamples. 
                            If not provided, then LHS samples are generated.
        """

        # Performing checks
        if len(self.DV) == 0:
            self._error("Add design variables before running the analysis.")

        # Sampling plan
        if self.options["sampling"] == "internal":

            # Validation
            if numSamples is None:
                self._error("Number of samples is not provided.")
            if not isinstance(numSamples, int):
                self._error("Number of samples argument is not an integer.")

            # Generate sampling plan
            samples = self.sampler(numSamples)
    
        elif self.options["sampling"] == "external":

            # Validation
            if doe is None:
                self._error("External samples are not provided.")
            if not isinstance(doe, np.ndarray):
                self._error("Provided external samples are not a numpy array.")
            if doe.ndim != 2:
                self._error("Provided external samples are not a 2D numpy array.")
            if doe.shape[1] != len(self.DV):
                self._error("Provided external samples are not of correct size.")

            # Assign user provided samples
            samples = doe
            numSamples = samples.shape[0]

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

            description.write("\nDesign Variable: {}".format(x))

            # Starting time
            t1 = time.time()

            try:
                # Getting output for specific sample
                output, field = self.getObjectives(x)

            except Exception as e:
                print(e)
                print("Error occured during the analysis. Check analysis.log in the respective folder for more details.")
                failed.append(sampleNo + 1)
                description.write("\nAnalysis failed.")

            else:
                # Check for analysis failure
                if output["fail"] == True: # Check for analysis failure
                    failed.append(sampleNo + 1)
                    description.write("\nAnalysis failed.")
                    
                # Check for implicit alpha
                elif self.options["alpha"] == "implicit" and abs(output["cl"] - self.options["targetCL"]) > self.options["targetCLTol"]: 
                    failed.append(sampleNo + 1)
                    description.write("\nAnalysis failed.")

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

    def calculateArea(self, x: np.ndarray) -> float:
        """
            Note: This function should not be called in the middle of analysis
            It should ONLY be used from outside. Do not use this method within 
            getObjectives. That method has its own implementation
            of area calculation.

            Function to calculate the area of the airfoil
            based on the value of design variable.

            Input:
            x - 1D numpy array (value of dv).

            Ouput:
            area: area of the airfoil.

            Note: To use this method, atleast lower or upper surface CST
            coefficient should be added as a DV.
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

        if "upper" not in self.DV and "lower" not in self.DV:
            self._error("\"upper\" or \"lower\" surface is not added as design variable.")

        # Creating dictionary from x
        newDV = {}
        for dv in self.DV:
            loc = self.locator == dv
            loc = loc.reshape(-1,)
            newDV[dv] = x[loc]

        # Updating the airfoil pointset based on new DV
        self.DVGeo.setDesignVars(newDV)

        # Getting the updated airfoil points
        points = self.DVGeo.update("airfoil")[:,0:2]
        x = points[:,0]
        y = points[:,1]

        # Calculate the area using simpson's rule
        # Note: x and y are both flipped here
        area = integrate.simpson(x, y, even='avg')

        return area
    
    # ----------------------------------------------------------------------------
    #                       Methods related to validation
    # ----------------------------------------------------------------------------

    def _checkOptions(self, defaultOptions: list, requiredOptions: list, options: dict) -> None:
        """
            This is a general method for validating options dictionary provided by
            user based on default and mandatory options.
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

        ############ Validating aeroProblem
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

    def _lhs(self, numSamples) -> np.ndarray:
        """
            Method to generate the lhs samples.
        """

        # Number of dimensions
        dim = len(self.lowerBound)

        # Generating normalized lhs samples
        samples = lhs(dim, samples=numSamples, criterion='cm', iterations=100*dim)

        # Scaling the samples
        x = self.lowerBound + (self.upperBound - self.lowerBound) * samples

        return x

    def _creatInputFile(self, x:np.ndarray) -> None:
        """
            Method to create an input file for analysis.
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

    def _writeCoords(self, coords, filename) -> None:
        """
            Writes out a set of airfoil coordinates in dat format.
        """

        # X and Y ccordinates of the airfoil
        x = coords[:, 0]
        y = coords[:, 1]

        with open(filename, "w") as f:
            for i in range(len(x)):
                f.write(str(round(x[i], 12)) + "\t\t" + str(round(y[i], 12)) + "\n")

        f.close()

    def _writeSurfMesh(self, coords, filename):
        """
            Writes out surface mesh in Plot 3D format (one element in z direction)
        """

        # X and Y ccordinates of the airfoil
        x = coords[:, 0]
        y = coords[:, 1]

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
                            f.write("%g\n" % (float(j)))

        f.close()

    def _plotAirfoil(self, plt, orig_airfoil, def_airfoil) -> None:
        """
            Method for plotting the base airfoil
            and the deformed airfoil.
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

            Inputs:
            message: a string containing the message to be displayed.
            type: an integer value to specify whether it is warning or error.
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

        exit()