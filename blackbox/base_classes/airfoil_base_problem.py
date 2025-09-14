from abc import abstractmethod
from .base_problem import BaseProblem
from dataclasses import dataclass
import os, pickle, psutil, time, shutil, sys
import numpy as np
from baseclasses import AeroProblem
from scipy.io import savemat
from mpi4py import MPI

from ..msg import print_msg
from pygeo import DVGeometry

class AirfoilBaseProblem(BaseProblem):
    """
        Base class for all airfoil problems.

        This class needs to be inherited by all the airfoil problems.

        Wherever this class is used, child class needs to initialize various variables for proper working of the class.
    """

    def __call__(self, x):

        # self.check_input(x)

        x = np.atleast_2d(x)

        for i in range(x.shape[0]):



            # Getting the deformed airfoil
            points = self.get_airfoil_coordinates(x)

            directory = self.options.directory



        print("Running analysis {}".format(self.genSamples + 1))

    def _run_analysis(self, x: np.ndarray) -> tuple:
        """
            Method for running an analysis for a given sample.

            Parameters
            ----------
            x: 1D numpy array
                design variable for the analysis

            Returns
            -------
            output: tuple
                A tuple containing the outputs generated from the analysis.
                First entry of the tuple is a dictionary containing scalar outputs
                from the analysis. Second entry is a dictionary containing the
                field data, if requested. Otherwise, it is None
        """

        # Checking if the appropriate options are set for analysis

        # assert 
        
        if self.options["solverOptions"] == {} or self.options["meshingOptions"] == {} or self.options["aeroProblem"] == None:
            self._error("You need to set solverOptions, meshingOptions and aeroProblem in the options dictionary for running the analysis.")

        # Overiding/set some solver options
        if self.options.solver == "adflow":
            self.options.solverOptions["printAllOptions"] = False
            self.options.solverOptions["printIntro"] = False
            self.options.solverOptions["outputDirectory"] = "."
            self.options.solverOptions["numberSolutions"] = False
            self.options.solverOptions["printTiming"] = False

        # Raise an error if pyvista is not installed
        if self.options.get_flowfield_data:
            self.options.solverOptions["writeSurfaceSolution"] = True

        # Getting the deformed airfoil
        points = self.get_airfoil_coordinates(x)

        print("Running analysis {}".format(self.samples_generated + 1))

        directory = self.options.directory

        # Create the folder for saving the results
        os.system("mkdir {}/{}".format(directory, self.samples_generated+1))

        # Getting the directory where package is saved
        pkgdir = sys.modules["blackbox"].__path__[0]

        # Setting filepath based on the solver
        if self.options.solver == "adflow":

            if type(self).__name__ == "AirfoilCSTMultipoint":
                filepath = os.path.join(pkgdir, "airfoil_problem/runscripts/runscript_airfoil_cst_mp.py")
            else:
                if self.options["alpha"] == "explicit":
                    filepath = os.path.join(pkgdir, "airfoil_problem/runscripts/runscript_airfoil.py")
                else:
                    filepath = os.path.join(pkgdir, "airfoil_problem/runscripts/runscript_airfoil_rf.py")

        elif self.options["solver"] == "dafoam":

            filepath = os.path.join(pkgdir, "airfoil_problem/runscripts/runscript_airfoil_dafoam.py")

        # Copy the runscript to analysis directory
        shutil.copy(filepath, "{}/{}/runscript.py".format(directory, self.samples_generated+1))

        # Copying openfoam standard folders
        if self.options.solver == "dafoam":
            ofdir = self.options.openfoam_directory
            os.system(f"cp -r {ofdir}/0 {directory}/{self.samples_generated+1}")
            os.system(f"cp -r {ofdir}/constant {directory}/{self.samples_generated+1}")
            os.system(f"cp -r {ofdir}/system {directory}/{self.samples_generated+1}")

        # Changing the directory to analysis folder
        os.chdir("{}/{}".format(directory, self.samples_generated+1))

        if self.options.write_airfoil_coordinates:
            self._write_coords(coords=points, filename="deformed_airfoil.dat")

        if self.options.plot_airfoil:
            self._plot_airfoil(self.plt, self.parametrization.orig_coords, points)

        if isinstance(self.parametrization, DVGeometry) and self.options.write_deformed_ffd:
            self.parametrization.writePlot3d("deformed_ffd.xyz")

        # Create input file
        self._creat_input_file(x)

        # Writing the surface mesh
        if self.options.solver == "dafoam":
            zSpan = 0.1
        elif self.options.solver == "adflow":
            zSpan = 1.0

        self._write_surf_mesh(coords=points, filename="surf_mesh.xyz", zSpan=zSpan)

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

            # Reading the output file containing results
            filehandler = open("output.pickle", 'rb')

        except Exception as e:
            print_msg(str(e), type=1)

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

            if self.options.get_flowfield_data:
                # Reading the cgns file
                filename = self.options.aero_problem.name + "_surf.cgns"
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
                field_data = {}

                for var in mesh.array_names:
                    # Skipping the first entry in the array
                    if var != "Base/Zone":
                        # set_active_scalars returns a tuple, and second
                        # entry contains the pyvista numpy array.
                        field_data[var] = np.asarray(mesh.set_active_scalars(var, "cell")[1])

            else:
                field_data = None

            return output, field_data
        
        finally:

            # Cleaning the directory
            files = ["vol_mesh.cgns", "input.pickle", "runscript.py",
                    "output.pickle", "fort.6", "opt.hst", "surf_mesh.xyz"]
            
            for file in files:
                if os.path.exists(file):
                    os.system("rm {}".format(file))

            if self.options.solver == "dafoam":
                os.system("rm -r 0")
                os.system("rm -r constant")
                os.system("rm -r system")

            # Changing the directory back to root
            os.chdir("../..")

            # Increase the number of generated samples
            self.samples_generated += 1

    def get_airfoil_coordinates(self, x: np.ndarray) -> np.ndarray:
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

        # Some validations
        try:
            assert len(self.design_variables) != 0, "Add design variables before running the analysis"
            assert isinstance(x, np.ndarray), "Input sample is not a numpy array"
            assert x.ndim == 1, "Input sample is a multi-dimensional array"
            assert len(x) == len(self.lower_bound), "Input sample is not of correct size"

        except AssertionError as e:
            print_msg(str(e))

        # If no geometric design variable is present, then return the original airfoil
        if self.parameterization.getNDV() == 0:
            return self.origCoords[:,0:2]

        # Creating dictionary from x
        new_design_variable = {}

        for dv in self.design_variables:
            loc = self.locator == dv
            loc = loc.reshape(-1,)
            new_design_variable[dv] = x[loc]

        if type(self).__name__ == "AirfoilCSTMultipoint":
            for ap in self.options.aero_problem:
                ap.setDesignVars(new_design_variable)

        if isinstance(self.parametrization,DVGeometry) and self.options.fix_LETE:

            # Adjusting LE FFD points
            midpoint = new_design_variable["shape"][0]/2
            new_design_variable["shape"][0] -= midpoint
            new_design_variable["shape"] = np.append(-midpoint, new_design_variable["shape"])

            # Adjusting TE FFD points
            midpoint = new_design_variable["shape"][-1]/2
            new_design_variable["shape"][-1] -= midpoint
            new_design_variable["shape"] = np.append(new_design_variable["shape"], -midpoint)

        # Updating the airfoil pointset based on new DV
        self.DVGeo.setDesignVars(new_design_variable)

        # Getting the updated airfoil points
        points = self.DVGeo.update("airfoil")[:,0:2]

        return points

    def calculate_area(self, x: np.ndarray) -> float:
        """
            Note: This function should not be called in the middle of analysis
            It should ONLY be used from outside. Do not use this method within 
            run_analysis. That method has its own implementation
            of area calculation

            Method to calculate the area of the airfoil
            based on the value of design variable

            Parameters
            ----------
            x: 1D numpy array
                design variable

            Returns
            -------
            area: float
                Area of the airfoil
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
    
    def _creat_input_file(self, x:np.ndarray) -> None:
        """
            Method to create an input file for a specific analysis

            Parameters
            ----------
            x: 1D numpy array
                design variable

            **Note**: This method is for internal use only
        """

        # Creating input dict
        input = {
            "solver_options": self.options.solver_options,
            "aero_problem": self.options.aero_problem,
            "meshing_options": self.options.meshing_options,
            "refine": self.options.refine,
            "write_slice_file": self.options.write_slice_file
        }

        # Adding non-shape design_variable
        for var_name in ["alpha", "mach", "altitude"]:
            if var_name in self.design_variables:
                loc = self.locator == f"{var_name}"
                loc = loc.reshape(-1,)
                input[f"{var_name}"] = x[loc]

        # Adding target Cl if alpha is implicit
        if self.options.alpha == "implicit":
            input["targetCL"] = self.options["targetCL"]
            input["targetCLTol"] = self.options["targetCLTol"]
            input["starting_alpha"] = self.options.starting_alpha

        # Saving the input file
        filehandler = open("input.pickle", "xb")
        pickle.dump(input, filehandler)
        filehandler.close()
        
    def _write_coords(self, coords: np.ndarray, filename: str) -> None:
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

    def _write_surf_mesh(self, coords, filename, zSpan=1.0):
        """
            Writes out surface mesh in Plot 3D format (only one element in z direction)

            Parameters
            ----------
            coords: 2D numpy array
                Airfoil coordinates.

            filename: str
                Name of the file to write the coordinates.

            zSpan: float
                Width in the z direction, default=1.0
            
            **Note**: This method is for internal use only
        """

        # X and Y ccordinates of the airfoil
        x = coords[:, 0]
        y = coords[:, 1]
        z = [0.0, zSpan]

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

    def _plot_airfoil(self, plt, orig_airfoil: np.ndarray, def_airfoil: np.ndarray) -> None:
        """
            Method for plotting the base airfoil
            and the deformed airfoil.

            Parameters
            ----------
            plt: matplotlib.pyplot
                matplotlib pyplot object.
            
            orig_airfoil: 2D numpy array
                original airfoil coordinates.

            def_airfoil: 2D numpy array
                deformed airfoil coordinates.
        """

        _, ax = plt.subplots()

        ax.plot(orig_airfoil[:,0], orig_airfoil[:,1], label="Original")
        ax.plot(def_airfoil[:,0], def_airfoil[:,1], label="Deformed")
        ax.set_xlabel("x/c", fontsize=14)
        ax.set_ylabel("y/c", fontsize=14)
        ax.legend(fontsize=12)
        plt.tight_layout()
        plt.savefig("airfoil.png", dpi=400)
        plt.close()
