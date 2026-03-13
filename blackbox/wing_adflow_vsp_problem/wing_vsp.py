import os
from importlib.metadata import version
import numpy as np

def parse_version(v):
    return tuple(int(x) for x in v.split("."))

class WingVSP():

    def __init__(self, vsp_file: str) -> None:
        """
            Class for parsing OpenVSP file and adding various wing shape parameters

            Parameters
            ----------
            vsp_file: str
                vsp3 file that contains the wing model

            `NOTE`: The vsp model should only contain one component of type `WING`
        """

        # check openvsp versions
        try:
            import openvsp
        except:
            raise RuntimeError(
                "OpenVSP Python API is not installed or can not be imported\n\n"
                "You can follow the installation guide: https://github.com/OpenVSP/OpenVSP/tree/main/src/python_api/packages"
            )
        else:
            MIN_VERSION = parse_version("3.28.0")
            CURR_VERSION = parse_version(version("openvsp"))
            assert CURR_VERSION >= MIN_VERSION, f"openvsp version is {CURR_VERSION} but minimum {MIN_VERSION} is required"

        # check pygeo versions
        try:
            from pygeo import DVGeometryVSP
        except:
            raise RuntimeError(
                "pyGeo is not installed or can not be imported\n\n"
                "You can follow the installation guide: https://mdolab-pygeo.readthedocs-hosted.com/en/latest/install.html"
            )
        else:
            MIN_VERSION = parse_version("1.17.0")
            CURR_VERSION = parse_version(version("pygeo"))
            assert CURR_VERSION >= MIN_VERSION, f"pygeo version is {CURR_VERSION} but minimum {MIN_VERSION} is required"

        assert isinstance(vsp_file, str), "`vsp_file` argument must be string"

        self.vsp_file = os.path.abspath(vsp_file) # set abs path

        # check if vsp file path exists
        if not os.path.exists(self.vsp_file):
            raise ValueError("Provided vsp file path does not exist")
        
        if hasattr(openvsp, "VSPVehicle"):
            self.vsp_model = openvsp.VSPVehicle()
        else:
            self.vsp_model = openvsp

        # clear the vsp model
        self.vsp_model.ClearVSPModel()

        # read the model file
        self.vsp_model.ReadVSPFile(self.vsp_file)

        # Some checks
        geoms = self.vsp_model.FindGeoms()
        assert(len(geoms) == 1) and self.vsp_model.GetGeomTypeName(geoms[0]) == "Wing", "Your OpenVSP model should contain only one component of type `WING`"

        # Storing some VSP variables
        self.component_id = geoms[0]
        self.component_name = self.vsp_model.GetGeomName(geoms[0])
        self.xsec_surf_id = self.vsp_model.GetXSecSurf(self.component_id, 0)
        self.number_of_xsec = self.vsp_model.GetNumXSec(self.xsec_surf_id)
        self.number_of_sections = self.number_of_xsec - 1

        # intialize some empty variables
        self.parameters = {}
        self.mask = np.array([])
        self.bounds = (np.array([]), np.array([]))

    def add_global_parameter(
        self,
        name: str,
        lower: float,
        upper: float,
    ) -> None:
        """
            Method for adding global wing parameters for analysis

            Some of the examples are span, dihedral, sweep

            Parameters
            ----------
            name: str
                name of the parameter to be added. It can only be:
                    * `Span`: span of the section
                    * `Dihedral`: dihedral of the section
                    * `Sweep`: sweep of the section

                `NOTE`: You can add a single parameter for all sections, refer to
                        `section_id` argument for more detail

            lower: float
                lower bound for the parameter

            upper: float
                upper bound for the parameter

            section_id: int
                an integer denoting to which section the parameters belongs to.
                The default is -1, which indicates that a single value controls
                this parameter across all the sections. `NOTE`: This parameter must
                not be more than number of sections in the VSP model
        """

        possible_global_paramaters = ["Span", "Dihedral", "Sweep"]

        assert isinstance(name, str), "`name` argument must be string"
        assert name in possible_global_paramaters, f"`{name}` is not a valid parameter"
        assert name not in self.parameters.keys(), f"`{name}` is already added as a parameter"
        assert isinstance(lower, float) and isinstance(upper, float), "`lower` and `upper` must be float values"
        assert lower < upper, "`upper` must be greater than `lower`"

        # empty list for storing pygeo openvsp parameter names
        dvgeo_name = []

        # store openvsp parameter names
        for i in range(self.number_of_sections):
            dvgeo_name.append(f"{self.component_name}:XSec_{i+1}:{name}")

        # update bounds
        lb = np.append(self.bounds[0], lower)
        ub = np.append(self.bounds[1], upper)

        self.mask = np.append(self.mask, name)
        self.bounds = (lb, ub)
        self.parameters[name] = dvgeo_name

    def add_local_parameter(
        self,
        name: str,
        lower: float,
        upper: float,
        section_id: int
    ) -> None:
        """
            Method for adding local/sectional wing parameters for analysis

            Some of the examples are twist, chord

            Parameters
            ----------
            name: str
                name of the parameter to be added. It can only be:
                    * `Span`: span of the section
                    * `Dihedral`: dihedral of the section
                    * `Sweep`: sweep of the section

                `NOTE`: You can add a single parameter for all sections, refer to
                        `section_id` argument for more detail

            lower: float
                lower bound for the parameter

            upper: float
                upper bound for the parameter

            section_id: int
                an integer denoting to which section the parameters belongs to

                `NOTE`: `section_id` must not be more than number of sections in the VSP model
        """

        possible_section_paramaters = ["Twist"]

        assert isinstance(name, str), "`name` argument must be string"
        assert name in possible_section_paramaters, f"`{name}` is not a valid parameter"
        assert name not in self.parameters.keys(), f"`{name}` is already added as a parameter"
        assert isinstance(lower, float) and isinstance(upper, float), "`lower` and `upper` must be float values"
        assert lower < upper, "`upper` must be greater than `lower`"
        assert isinstance(section_id, int) and section_id >=0 and section_id <= self.number_of_xsec, f"`section_id` must be an integer between 0 and {self.number_of_xsec}"

        # update bounds
        lb = np.append(self.bounds[0], lower)
        ub = np.append(self.bounds[1], upper)

        # set updated variables
        self.mask = np.append(self.mask, name)
        self.bounds = (lb, ub)
        self.parameters[name] = [f"{self.component_name}:XSec_{section_id}:{name}"]

    def get_pygeo_parameter_dict(self, x: np.ndarray) -> dict:
        """
            Method to get pygeo openvsp parameter dictionary for a given
            set of parameters

            This method is useful when you want to create a pygeo openvsp
            object from scratch and add design variables to it

            Parameters
            ----------
            x: np.ndarray
                1D numpy array containing a single set of parameters
        """

        assert isinstance(x, np.ndarray) and x.ndim == 1, "`x` should be a 1D numpy array"
        assert x.shape[0] == self.bounds[0].shape[0], "provided `x` does not have correct number of parameters"

        dvgeo_params = {}

        for key, vals in self.parameters.items():
            for val in vals:
                if key == "Span":
                    dvgeo_params[val] = x[self.mask == key].item()/self.number_of_sections
                else:
                    dvgeo_params[val] = x[self.mask == key].item()

        return dvgeo_params
    
    ###########################################################
    ######### Below methods are for internal use only #########
    ###########################################################
