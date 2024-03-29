# Importing python packages
import os
import numpy as np
from scipy.io import savemat
from ..base import BaseClass

class Forrester(BaseClass):
    """
        Class contains essential methods for generating data
        of Forrester Function:

            y(x) = (6x - 2)**2 * sin(12x - 4)

        There are two values possible for type: "single" and "multi".
        Options argument is not needed for "single" type analysis, while it
        is needed for "multi" type analysis. Provide an appropriate options
        dict which contains necessary attributes. There are only four possible
        attributes:

        "directory" : Folder name where the data.mat file will be saved (string, optional).
        "numberOfSamples" : number of samples to be generated (integer).

        Note: There is no sampling method option for forrester function since
        there is only one dimension.
    """

    def __init__(self, type="multi", options=None):

        self._initialization(type, options)

    # ----------------------------------------------------------------------------
    #               All the methods for multi analysis
    # ----------------------------------------------------------------------------

    def _setupMultiAnalysis(self, options):
        """
            Method to perform initialization when the type is "multi".
        """

        # Creating an empty options dictionary and assign value of type
        self.options = {}
        self.options["type"] = "multi"

        # Setting up default options
        self._getDefaultOptions()

        # Creating list of required options
        requiredOptions = ["numberOfSamples"]

        # Validating user provided options
        self._checkOptionsForMultiAnalysis(options, requiredOptions)

        self.options["lowerBound"] = 0
        self.options["upperBound"] = 1

        # Updating/Appending the default option list with user provided options
        self._setOptions(options)

        # Setting up the folder for saving the results
        self._setDirectory()

    # def _checkOptionsForMultiAnalysis(self, options):
    #     """
    #         This method validates user provided options for type = "multi".
    #     """

    #     # Creating list of various different options
    #     defaultOptions = list(self.options.keys())
    #     requiredOptions = ["numberOfSamples", "lowerBound", "upperBound"]
    #     allowedUserOptions = defaultOptions
    #     allowedUserOptions.extend(requiredOptions)

    #     userProvidedOptions = list(options.keys())

    #     # Checking if the user provided option contains only allowed attributes
    #     if not set(userProvidedOptions).issubset(allowedUserOptions):
    #         self._error("Option dictionary contains unrecognized attribute(s).")

    #     # Checking if user has mentioned all the requried attributes
    #     if not set(requiredOptions).issubset(userProvidedOptions):
    #         self._error("Option dictionary doesn't contain all the requried options. \
    #                     {} attribute(s) is/are missing.".format(set(requiredOptions) - set(userProvidedOptions)))

    #     # Checking type of required options
    #     for attribute in requiredOptions:
    #         if type(options[attribute]) is not int:
    #             self._error("\"{}\" attribute is not an integer".format(attribute))

    #     # Setting minimum limit on number of samples
    #     if options["numberOfSamples"] < 2:
    #         self._error("Number of samples need to least 2.")

    #     # Checking correctness of bound
    #     if options["lowerBound"] >= options["upperBound"]:
    #         self._error("Lower bound is greater than upper bound.")

    #     # Validating directory attribute
    #     if "directory" in userProvidedOptions:
    #         if type(options["directory"]) is not str:
    #             self._error("\"directory\" attribute is not string.")

    def generateSamples(self):
        """
            Method to generate samples and save the data for further use.
        """

        if self.options["type"] != "multi":
            self._error("You cannot call generateSamples() method when type is not \"multi\".")

        self.samples = np.linspace(self.options["lowerBound"], self.options["upperBound"], self.options["numberOfSamples"])
        self.samples = np.reshape(self.samples, (-1,1))

        self.y = self._function(self.samples)

        data = {"x" : self.samples, "y" : self.y }

        print("{} samples generated.".format(self.options["numberOfSamples"]))

        # Saving data file in the specified folder
        os.chdir(self.options["directory"])
        savemat("data.mat", data)
        os.chdir("../")

    # ----------------------------------------------------------------------------
    #           All the methods for single analysis
    # ----------------------------------------------------------------------------

    def _setupSingleAnalysis(self):
        """
            Method to setup object for single analysis
        """
        
        self.options = {}
        self.options["type"]  = "single"

        # There is no other special initialization needed for 
        # Forrester function when the type is "single".

    def getObjectives(self, x):
        """
            Method to generate a single output y, based on input x.
            Input x should be an integer. Output y will also be an integer.
        """

        if self.options["type"] != "single":
            self._error("You cannot call getObjectives() method when type is not \"single\".")

        # if type(x) != int and type(x) != float:
        #     self._error("Provided x is not an integer.")

        return self._function(x)

    # ----------------------------------------------------------------------------
    #          Other required methods, irrespective of type of analysis.
    # ----------------------------------------------------------------------------

    def _function(self, x):
        """
            Forrester function. Note: Output of function should always be
            of size num_samples X 1. Here, input x is already of size num_samples X 1,
            so no need to use reshape.
        """

        y = (6*x - 2)**2 * np.sin(12*x - 4)

        return y
