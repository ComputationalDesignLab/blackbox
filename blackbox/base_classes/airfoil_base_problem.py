from abc import abstractmethod
from .base_problem import BaseProblem
from dataclasses import dataclass
import os, pickle, psutil, time, shutil, sys
import numpy as np
from baseclasses import AeroProblem
from scipy.io import savemat
from mpi4py import MPI
from typing import Optional

@dataclass
class AirfoilOptions:

    # CFD solver and Meshing options
    solver: str = "adflow"
    solver_options: dict = {}
    meshing_options: dict = {}
    aero_problem: Optional[AeroProblem] = None
    open_foam_directory: str = "."

    # Other options
    directory: str = "output"


class AirfoilBaseProblem(BaseProblem):
    """
        Base class for all airfoil problems.

        This class needs to be inherited by all the airfoil related classes.

        Wherever this class is used, child class needs to initialize various
        variables for proper working of the class.
    """

    def __call__(self):

        pass

    @abstractmethod
    def _evaluate(self):

        pass

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
            assert len(self.dv) != 0, "Add design variables before running the analysis"
            assert isinstance(x, np.ndarray), "Input sample is not a numpy array"
            assert x.ndim == 1, "Input sample is a multi-dimensional array"
            assert len(x) == len(self.lower_bound), "Input sample is not of correct size"

        except AssertionError as e:
            self._error(str(e))
            
        # If no geometric design variable is present, then return the original airfoil
        if self.DVGeo.getNDV() == 0:
            return self.origCoords[:,0:2]

        # Creating dictionary from x
        new_dv = {}
        for dv in self.DV:
            loc = self.locator == dv
            loc = loc.reshape(-1,)
            new_dv[dv] = x[loc]

        if type(self).__name__ == "AirfoilCSTMultipoint":
            for ap in self.options["aeroProblem"]:
                ap.setDesignVars(new_dv)

        if self.parametrization == "FFD" and self.options["fixLETE"]:

            # Adjusting LE FFD points
            midpoint = new_dv["shape"][0]/2
            new_dv["shape"][0] -= midpoint
            new_dv["shape"] = np.append(-midpoint, new_dv["shape"])

            # Adjusting TE FFD points
            midpoint = new_dv["shape"][-1]/2
            new_dv["shape"][-1] -= midpoint
            new_dv["shape"] = np.append(new_dv["shape"], -midpoint)

        # Updating the airfoil pointset based on new DV
        self.DVGeo.setDesignVars(new_dv)

        # Getting the updated airfoil points
        points = self.DVGeo.update("airfoil")[:,0:2]

        return points
