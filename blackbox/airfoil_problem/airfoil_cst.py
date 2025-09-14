import os
from ..base_classes.airfoil_base_problem import *
from .cst import CST
from .airfoil_options import AirfoilOptions
from ..msg import print_msg

class AirfoilCST(AirfoilBaseProblem):
    """
        Class for generating airfoil data using CST parameterization
    """

    def __init__(self, options: AirfoilOptions):

        try:
            assert isinstance(options, AirfoilOptions), "options argument should be an object of AirfoilOptions class"

        except Exception as e:
            print_msg(str(e))

        self.options = options

        # Getting abs path for the storage directory
        self.options.directory = os.path.abspath(self.options.directory)

        # Creating directory for storing the results
        if not os.path.isdir(self.options.directory):
            os.system("mkdir {}".format(self.options.directory))
        else:
            os.system("rm -r {}".format(self.options.directory))
            os.system("mkdir {}".format(self.options.directory))

        # Read the coordinate file
        self.orig_coords = CST.read_coord_file(self.options.airfoil_file)

        # Initializing the parametrization object
        self.parametrization = CST(self.options.airfoil_file, num_cst=[self.options.num_cst_upper, self.options.num_cst_lower])

        # Some initializations which will be used later
        self.design_variables = []
        self.samples_generated = 0

    def add_design_varaible(self, name: str, lower_bound: list, upper_bound: list) -> None:
        """
            Method for adding a design variable within CST parameterization

            Parameters
            ----------
            name : str
                Name of the design variable, valid names are: "upper", "lower", "alpha", "mach" or "altitude"

            lower_bound : float or numpy.ndarray
                Lower bound of the design variable

            upper_bound : float or numpy.ndarray
                Upper bound of the design variable
        """

        # Checking
        self._check_design_variable(name, lower_bound, upper_bound)

        if name == "upper" or name == "lower":
            locator = np.array(["{}".format(name)]*len(lower_bound))

            if len(self.design_variables) == 0:
                self.upper_bound = upper_bound
                self.lower_bound = lower_bound
                self.locator = locator
            else:
                self.upper_bound = np.append(self.upper_bound, upper_bound)
                self.lower_bound = np.append(self.lower_bound, lower_bound)
                self.locator = np.append(self.locator, locator)
        else:
            locator = np.array(["{}".format(name)])

            if len(self.design_variables) == 0:
                self.upper_bound = np.array([upper_bound])
                self.lower_bound = np.array([lower_bound])
                self.locator = np.array([locator])
            else:
                self.upper_bound = np.append(self.upper_bound, upper_bound)
                self.lower_bound = np.append(self.lower_bound, lower_bound)
                self.locator = np.append(self.locator, locator)    

        # Adding design varialbe within CST parametrization
        if name not in ["alpha", "mach", "altitude"]:
            self.parametrization.add_design_variable(name)

        # Adding the DV to the list
        self.design_variables.append(name)

    def remove_design_variable(self, name: str) -> None:
        """
            Method to remove a design variable 

            Parameters
            ----------
            name : str
                Name of the design variable. It can be "upper", "lower", "alpha", "mach" or "altitude".
        """

        if name not in self.design_variables:
            print_msg("{} doesn't exists as a design variable.".format(name))

        # Finding the indices to be removed
        loc = self.locator == name
        loc = loc.reshape(-1,)

        # Removing the entry in the bounds and locator
        self.lower_bound = np.delete(self.lower_bound, loc)
        self.upper_bound = np.delete(self.upper_bound, loc)
        self.locator = np.delete(self.locator, loc)

        # Remove design variable from CST parameterization
        if name not in ["alpha", "mach", "altitude"]:
            self.parametrization.remove_design_variable(name)

        # Removing the entry from design variable list
        self.design_variables.remove(name)

    def _check_design_variable(self, name: str, lb, ub) -> None:
        """
            Method for validating a design variable.

            Parameters
            ----------

            name : str
                Name of the design variable, it can be "upper", "lower", "alpha", "mach" or "altitude".

            lb : float or numpy.ndarray
                Lower bound of the design variable

            ub : float or numpy.ndarray
                Upper bound of the design variable
        """

        # List of possible design variables
        possible_design_variables = ["upper", "lower", "alpha", "mach", "altitude"]

        try:
            assert isinstance(name, str), "name argument should be a string"
            assert name in possible_design_variables, f"{name} is not a valid design variable"
            assert name not in self.design_variables, f"{name} already exists as a design variable"

            if name == "alpha":
                assert self.options.alpha == "explicit", "Alpha cannot be a design variable when \"alpha\" attribute in option is \"implicit\""

            if name in ["mach", "altitude"]:
                assert self.options.aero_problem is not None, f"Aero problem needs to be defined in options to set \"{name}\" as design variable"
                assert name in self.options.aero_problem.inputs.keys(), f"You need to initialize \"{name}\" in the aero problem to set it as design variable"

            # Validating the bounds for upper and lower variable
            if name in ["upper", "lower"]:
                assert isinstance(lb, np.ndarray) and lb.ndim == 1, f"Lower bound for \"{name}\" variable should be a 1D numpy array"
                assert isinstance(ub, np.ndarray) and ub.ndim == 1, f"Upper bound for \"{name}\" variable should be a 1D numpy array"
                assert len(lb) == self.options.num_cst_upper and len(ub) == self.options.num_cst_upper, f"Length of lower and upper bound for \"{name}\" variable is not equal to number of CST coeff"
                assert np.all(lb < ub), f"Lower bound is greater than or equal to upper bound for \"{name}\" variable"

            # Validating lb and ub of the scalar DVs
            else:
                assert isinstance(lb, float), f"Lower bound for \"{name}\" variable should be a float"
                assert isinstance(ub, float), f"Upper bound for \"{name}\" variable should be a float"
                assert lb < ub, f"Lower bound is greater than or equal to upper bound for \"{name}\" variable"

        except Exception as e:
            print_msg(str(e))
