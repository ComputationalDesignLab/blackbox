import os, sys, pickle, psutil, h5py, json
import numpy as np
from mpi4py import MPI
from time import time

from .cst import CST
from .utils import AirfoilADflowOptions
from ..base_classes.base_problem import BaseProblem

comm = MPI.COMM_WORLD

class AirfoilADflow(BaseProblem):

    def __init__(self, options: AirfoilADflowOptions) -> None:
        """
            Class for performing airfoil analysis using ADflow solver
            
            CST is used for parametrizing the airfoil in this problem

            pyHyp is used for creating mesh around the airfoil

            Refer to the documentation for more details about the analysis pipeline

            Parameters
            ----------
            options: AirfoilOptions
                airfoil option dataclass object
        """

        assert isinstance(options, AirfoilADflowOptions), "options argument should be an object of AirfoilOptions class"

        self.options = options

        # Creating directory for storing the results
        if not os.path.isdir(self.options.directory):
            os.system("mkdir {}".format(self.options.directory))
        else:
            os.system("rm -r {}".format(self.options.directory))
            os.system("mkdir {}".format(self.options.directory))

        # Initializing the parametrization object
        self.parametrization = CST(self.options.airfoil_file, num_cst=[self.options.num_cst_upper, self.options.num_cst_lower])

        # Some initializations which will be used later
        self.parameters = []
        self.mask = np.array([])
        self.bounds = (np.array([]), np.array([]))
        self.samples_generated = 0

    def add_parameter(self, name: str, lower_bound: float | np.ndarray, upper_bound: float | np.ndarray) -> None:
        """
            Method for adding a parameter for the airfoil problem

            Parameters
            ----------
            name: str
                name of the parameter to be added. It can only be:
                    * `lower_cst`: lower surface cst coefficients
                    * `upper_cst`: upper surface cst coefficients
                    * `alpha`: angle of attack of the flow
                    * `mach`: mach number of the flow
                    * `reynolds`: reynolds number for the flow

            lower_bound: float or np.ndarray
                lower bound for the parameter. It should be a 1D numpy array
                if the parameter is `lower` or `upper`, otherwise a float

            upper_bound: float or np.ndarray
                upper bound for the parameter. It should be a 1D numpy array
                if the parameter is `lower` or `upper`, otherwise a float
        """

        self._check_variable(name, lower_bound, upper_bound)

        if name.lower() in ["upper_cst", "lower_cst"]:
            mask = [f"{name.lower()}"] * len(lower_bound)
            lb = np.append(self.bounds[0], lower_bound)
            ub = np.append(self.bounds[1], upper_bound)

        else:
            mask = [f"{name.lower()}"]
            lb = np.append(self.bounds[0], np.array([lower_bound]))
            ub = np.append(self.bounds[1], np.array([upper_bound]))
            
        self.bounds = (lb, ub)
        self.mask = np.append(self.mask, mask)
        self.parameters.append(name.lower())

    def __call__(self, x: np.ndarray, return_results: bool = False) -> None | dict:
        """
            Method to evaluate given x. It can be a single or multiple samples of
            size (d,) or (N,d) where N is the number of samples and d is the number of
            parameters added

            Parameters
            ----------
            x: np.ndarray
                a numpy array representing values for added parameters

            return_results: bool
                flag to determine if the results should be returned or not.
                This should be set to True only when you want this function to
                return the results.
                
                `NOTE`: This is only useful when you are doing sequential
                data generation
        """

        assert len(self.parameters) > 0, "add some parameters before running analysis"
        assert isinstance(x, np.ndarray), "given sample 'x' should be a numpy array"

        x = np.atleast_2d(x)

        assert x.shape[1] == self.bounds[0].shape[0], "size of given sample 'x' is not same as the number of parameters"

        # Creating and writing a description file
        description = open("{}/description.txt".format(self.options.directory), "a", buffering=1)

        if self.samples_generated == 0:

            description.write("---------------------------------------------------")
            description.write("\nAirfoil analysis using ADflow")
            description.write("\n--------------------------------------------------")
            description.write(f"\nVariables: {self.parameters}")
            description.write(f"\nLower bound for design variables:\n{self.bounds[0]}")
            description.write(f"\nUpper bound for design variables:\n{self.bounds[1]}")
            description.write("\n-----------------------------")
            description.write("\nAnalysis specific description")
            description.write("\n-----------------------------")

        if return_results:
            output = {}

        for i in range(x.shape[0]):

            description.write(f"\nAnalysis {self.samples_generated+1}:")

            t1 = time()

            try:
                self._run_analysis(x[i,:])

                if return_results:
                    output[self.samples_generated+1] = self.read_results(f"{self.options.directory}/{self.samples_generated+1}")
                    output[self.samples_generated+1]["parameters"] = x[i,:]

            except Exception as e:
                description.write(f"\n Error: {e}")
                description.write(f"\n -------- Analysis failed ----------")

                if return_results:
                    output[f"{i}"] = {}
                    output[f"{i}"]["scalar"] = {"fail": True}
                
            finally:
                # Write time taken for analysis to desc file
                description.write(f"\nTime taken for analysis: {(time()-t1)/60} min.\n")

                self.samples_generated += 1

        description.close() # close the description file

        if return_results:
            return output

    def calculate_area(self, x: np.ndarray) -> float:
        """
            Method to calculate the area of an airfoil shape based on given parameters

            Parameters
            ----------
            x: 1D numpy array
                design variable

            Returns
            -------
            area: float
                set of parameters based on which airfoil area will be computed
        """

        # Getting the updated airfoil points
        points = self._get_airfoil(x)

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
    
    def read_results(self, dir_name: str) -> dict:
        """
            Method to read scalars.json and field_outputs.hdf5 files containing results from an
            airfoil analysis and convert it to a dictionary which can be used later

            Parameters
            ----------
            dir_name: str
                name of directory where the results from the analysis are saved

            Returns
            -------
                out: a dictionary containing data stored in the folder
        """

        assert isinstance(dir_name, str), "file name should be a string"

        output = {}

        with open(f"{dir_name}/scalar_outputs.json", "r") as fp:
            scalar_data = json.load(fp)
        fp.close()

        output["scalar"] = scalar_data

        # read the field outputs hdf5 file
        f = h5py.File(f"{dir_name}/field_outputs.hdf5", "r")
        for key in list(f.keys()): # only for reading a group
            output[key] = {}
            for k in list(f[key].keys()):
                output[key][k] = f[key][k][()]
        f.close()

        return output
    
    ###########################################################
    ######### Below methods are for internal use only #########
    ###########################################################
    
    def _run_analysis(self, x: np.ndarray) -> None:
        """
            Method to evaluate a given set of parameters

            Parameters
            ----------
            x: np.ndarray
                a 1D numpy array representing values for added parameters
        """

        assert len(self.parameters) > 0, "add some parameters before running analysis"
        assert isinstance(x, np.ndarray) and x.ndim == 1, "given sample 'x' should be a 1D numpy array"
        assert x.shape[0] 

        print("Running analysis {}".format(self.samples_generated + 1))

        directory = self.options.directory

        # Create the folder for saving the results
        os.system("mkdir {}/{}".format(directory, self.samples_generated + 1))

        # Getting the directory where package is saved
        pkgdir = sys.modules["blackbox"].__path__[0]
        filepath = os.path.join(pkgdir, "airfoil_adflow_problem/runscript.py")

        # Copy the runscript to analysis directory
        os.system(f"cp {filepath} {directory}/{self.samples_generated + 1}/runscript.py")

        # Changing the directory to analysis folder
        os.chdir("{}/{}".format(directory, self.samples_generated + 1))

        points = self._get_airfoil(x)

        if self.options.write_airfoil_coordinates:
            self._write_coords(coords=points, filename="deformed_airfoil.dat")

        if self.options.plot_airfoil:
            self._plot_airfoil(self.parametrization.orig_coords, points)

        # write surface mesh
        self._write_surf_mesh(coords=points, filename="surf_mesh.xyz")

        # write parameters
        self._write_parameters(x)

        # Create input file
        self._create_input_file(x)

        try:
            # Spawning the runscript on desired number of processors
            child_comm = MPI.COMM_SELF.Spawn(sys.executable, args=["runscript.py"], maxprocs=self.options.num_processors)

            # Creating empty process id list
            pid_list = []

            # Getting each spawned process
            for processor in range(self.options.num_processors):
                pid = child_comm.recv(source=MPI.ANY_SOURCE, tag=processor)
                pid_list.append(psutil.Process(pid))

            # Disconnecting from intercommunicator
            child_comm.Disconnect()

            # Waiting till all the child processors are finished
            while len(pid_list) != 0:
                for pid in pid_list:
                    if not pid.is_running():
                        pid_list.remove(pid)

        except Exception as e: 
            print(e)

        finally:

            # Cleaning the directory
            files = ["vol_mesh.cgns", "input.pickle", "runscript.py", "surf_mesh.xyz"] 
            
            for file in files:
                if os.path.exists(file):
                    os.system(f"rm {file}")

            # Changing the directory back to root
            os.chdir("../..")

    def _get_airfoil(self, x: np.ndarray) -> np.ndarray:
        """
            Method for getting the airfoil coordinates for a set of parameters

            Parameters
            ----------
            x: 1D numpy array
                set of parameters based on which airfoil shape will be computed

            Returns
            -------
            points: 2D numpy array
                airfoil coordinates representing deformed airfoil shape
        """

        # if no cst variables, then return original coordinates
        if "upper_cst" not in self.parameters and "lower_cst" not in self.parameters:
            return self.parametrization.orig_coords
        
        # upper surface cst coeff
        if "upper_cst" in self.parameters:
            upper_cst_coeff = x[self.mask == "upper_cst"]
        else:
            upper_cst_coeff = self.parametrization.upper_cst
        
        # lower surface cst coeff
        if "lower_cst" in self.parameters:
            lower_cst_coeff = x[self.mask == "lower_cst"]
        else:
            lower_cst_coeff = self.parametrization.lower_cst

        coords = self.parametrization.compute_airfoil_coordinates(upper_cst_coeff, lower_cst_coeff)

        return coords
    
    def _write_coords(self, coords: np.ndarray, filename: str) -> None:
        """
            Method to write airfoil coordinates in selig format in a dat file

            Parameters
            ----------
            coords: 2D numpy array
                airfoil coordinates that will be written

            filename: str
                name of the file to write the coordinates
        """

        # X and Y ccordinates of the airfoil
        x = coords[:, 0]
        y = coords[:, 1]

        with open(filename, "w") as f:
            for i in range(len(x)):
                f.write(str(round(x[i], 12)) + "\t\t" + str(round(y[i], 12)) + "\n")

        f.close()

    def _plot_airfoil(self, orig_airfoil: np.ndarray, def_airfoil: np.ndarray) -> None:
        """
            Method for plotting the base airfoil and the deformed airfoil

            Parameters
            ----------
            orig_airfoil: 2D numpy array
                original airfoil coordinates

            def_airfoil: 2D numpy array
                deformed airfoil coordinates
        """

        import matplotlib.pyplot as plt

        _, ax = plt.subplots()

        ax.plot(orig_airfoil[:,0], orig_airfoil[:,1], label="Original airfoil")
        ax.plot(def_airfoil[:,0], def_airfoil[:,1], label="Deformed airfoil")
        ax.set_xlabel("x/c", fontsize=14)
        ax.set_ylabel("y/c", fontsize=14)
        ax.legend(fontsize=12)

        plt.savefig("airfoil.png", dpi=400)

        plt.close()

    def _write_surf_mesh(self, coords: np.ndarray, filename: str) -> None:
        """
            Method to write surface mesh in Plot 3D format (only one element in z direction)

            Parameters
            ----------
            coords: np.ndarray
                a 2D numpy array representing x and y coordinates of the airfoil

            filename: str
                name of the file to write the coordinates
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

    def _create_input_file(self, x: np.ndarray) -> None:
        """
            Method to create an input file for a specific analysis
            for a given parameter set

            Parameters
            ----------
            x: np.ndarray
                1D numpy array representing a single set of parameters
        """

        # Creating input dict
        input = {
            "solver_options": self.options.solver_options,
            "aero_problem": self.options.aero_problem,
            "meshing_options": self.options.meshing_options,
            "refine": self.options.refine,
            "write_slice_file": self.options.write_slice_file,
            "scalar_outputs": self.options.scalar_outputs,
            "alpha_type": self.options.alpha,
            "target_CL": self.options.target_CL,
            "target_CL_tol": self.options.target_CL_tol,
            "starting_alpha": self.options.starting_alpha
        }

        # Adding non-shape DV
        if "alpha" in self.parameters:
            mask = self.mask == "alpha"
            mask = mask.reshape(-1,)
            input["alpha"] = x[mask]

        if "mach" in self.parameters:
            mask = self.mask == "mach"
            mask = mask.reshape(-1,)
            input["mach"] = x[mask]

        if "reynolds" in self.parameters:
            mask = self.mask == "reynolds"
            mask = mask.reshape(-1,)
            input["reynolds"] = x[mask]

        # Saving the input file
        filehandler = open("input.pickle", "xb")
        pickle.dump(input, filehandler)
        filehandler.close()

    def _write_parameters(self, x):
        """
            Method to write a set of prameters `x` to a json file
        """

        assert len(self.parameters) > 0, "add some parameters before calling this method"

        parameters = {}
        for name in self.parameters:
            mask = self.mask == name
            if name == "lower_cst" or name == "upper_cst":
                parameters[name] = x[mask].tolist()
            else:
                parameters[name] = x[mask].item()

        with open("parameters.json", "w") as fp:
            json.dump(parameters, fp, indent=4)
        fp.close()

    def _check_variable(self, name: str, lower_bound: float | np.ndarray, upper_bound: float | np.ndarray):
        """
            Method for validating a given parameter before adding it.

            Parameters
            ----------
            name: name of the parameter. It can be "upper_cst", "lower_cst", 
                "alpha", "mach" or "reynolds"

            lb: lower bound of the parameter

            ub: upper bound of the parameter
        """

        # List of possible parameters
        valid_parameters = ["upper_cst", "lower_cst", "alpha", "mach", "reynolds"]

        # Validating name of the parameter
        assert isinstance(name, str), "'name' argument must be a string"

        assert name.lower() in valid_parameters, f"'{name}' is not a valid parameter"

        assert name.lower() not in self.parameters, f"'{name}' is already added as a parameter"

        if name.lower() == "alpha":
            assert self.options.alpha == "explicit", "'alpha' cannot be a parameter when \"alpha\" attribute in options is 'implicit'"

        if name.lower() in ["mach", "reynolds"]:
            assert name.lower() in self.options.aero_problem.inputs.keys(), f"initialize '{name}' in the aero problem to set it as parameter"

        # Validating bounds
        for bound in [lower_bound, upper_bound]:
            if name.lower() in ["upper_cst", "lower_cst"]:
                assert isinstance(bound, np.ndarray) and bound.ndim == 1, f"lower and upper bound should be a 1D numpy array if name is 'lower_cst' or 'upper_cst'"
            else:
                assert isinstance(bound, float), f"lower and upper bound should be float if name is 'alpha', 'mach' or 'reynolds'"
            
        # Validating shape of bounds
        if name.lower() == "upper_cst":
            assert upper_bound.shape[0] == lower_bound.shape[0] == self.options.num_cst_upper, "length of upper and/or lower bound is not same as the number of CST coefficiets for upper surface"

        if name.lower() == "lower_cst":
            assert upper_bound.shape[0] == lower_bound.shape[0] == self.options.num_cst_lower, "length of upper and/or lower bound is not same as the number of CST coefficiets for lower surface"

        assert np.all(upper_bound > lower_bound), "upper bound should be greater than lower bound"
