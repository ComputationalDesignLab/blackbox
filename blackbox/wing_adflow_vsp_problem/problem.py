import os, sys, pickle, json, h5py
import numpy as np
import pyvista as pv
from time import time
from packaging import version
from .options import WingADflowVSPOptions
from ..base_classes.base_problem import BaseProblem

class WingADflowVSP(BaseProblem):

    def __init__(self, options: WingADflowVSPOptions):
        """
            Class for performing wing analysis using ADflow solver
            
            OpenVSP is used for parametrizing the airfoil in this problem

            idwarp is used for deforming the given volume mesh based on the wing geometry changes

            Refer to the documentation for more details about the analysis pipeline

            Parameters
            ----------
            options: WingADflowVSPOptions
                wing vsp option dataclass object
        """

        assert isinstance(options, WingADflowVSPOptions), "options argument should be an object of AirfoilOptions class"

        self.options = options

        # check if required packages are available with specific versions
        try:
            import adflow
        except:
            raise RuntimeError(
                "ADFLOW solver is not installed or can not be imported\n\n"
                "You can follow the installation guide: https://mdolab-adflow.readthedocs-hosted.com/en/latest/install.html"
            )
        else:
            MIN_VERSION = version.parse("2.11.0")
            CURR_VERSION = version.parse(adflow.__version__)
            assert CURR_VERSION >= MIN_VERSION, f"ADFlow version is {CURR_VERSION} but minimum {MIN_VERSION} is required"

        if self.options.write_surface_output or self.options.write_volume_output:
            try:
                import pyvista, h5py
            except ImportError as e:
                raise ValueError(
                    "`pyvista` and `h5py` are not installed, it is required when surface or volume output is set to `True`"
                    "Run this command to install both these pacages: pip install pyvista h5py"
                ) from e

        # Creating directory for storing the results
        if not os.path.isdir(self.options.directory):
            os.system("mkdir {}".format(self.options.directory))
        else:
            os.system("rm -r {}".format(self.options.directory))
            os.system("mkdir {}".format(self.options.directory))

        # ensure alpha is not added as a DV
        if self.options.alpha == "implicit":
            for key, val in self.options.aero_problem.DVs.items():
                assert val.key != "alpha"

        # Some initializations which will be used later
        # self.parameters = []
        # self.mask = np.array([])
        # self.bounds = (np.array([]), np.array([]))
        self.samples_generated = 0

    @property
    def bounds(self):
        """
            Method to get upper and lower bound based on added parameters
        """

        lb = np.array([])
        ub = np.array([])

        for key, val in self.options.aero_problem.DVs.items():
            lb = np.append(lb, val.lower)
            ub = np.append(ub, val.upper)

        if self.options.vsp_file is not None:

            lb = np.append(lb, self.options.vsp_file.bounds[0])
            ub = np.append(ub, self.options.vsp_file.bounds[1])

        return (lb, ub)

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

        assert len(self.bounds[0].shape[0]) > 0, "add some parameters before running analysis"
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

    def read_results(self):

        pass

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

        import psutil
        from mpi4py import MPI

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
            "wing_vsp": self.options.vsp_file,
            "aero_problem": self.options.aero_problem,
            "write_slice_file": self.options.write_slice_file,
            "scalar_outputs": self.options.scalar_outputs,
            "alpha_type": self.options.alpha,
            "target_CL": self.options.target_CL,
            "target_CL_tol": self.options.target_CL_tol,
            "starting_alpha": self.options.starting_alpha
        }

        # Saving the input file
        filehandler = open("input.pickle", "xb")
        pickle.dump(input, filehandler)
        filehandler.close()

    def _write_parameters(self, x):
        """
            Method to write a set of prameters `x` to a json file
        """

        assert len(self.bounds[0].shape[0]) > 0, "add some parameters before calling this method"
        assert isinstance(x, np.ndarray) and x.ndim == 1, "x must be a 1D numpy array"
        assert x.shape[0] == self.bounds[0].shape[0], "number of entries in x is not same as the number of parameters"

        # separate flow and shape variables
        x_flow = x[:len(self.options.aero_problem.DVs)]
        x_shape = x[len(self.options.aero_problem.DVs):]

        for val in x_flow:

            print

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