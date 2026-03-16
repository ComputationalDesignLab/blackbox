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
        assert isinstance(lower, float) and isinstance(upper, float), "`lower` and `upper` must be float values"
        assert lower < upper, "`upper` must be greater than `lower`"
        for key in self.parameters.keys():
            assert name not in key, f"`{name}` is already added as a global/local parameter"

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
                an integer denoting which section the parameters belongs to

                `NOTE`: `section_id` must not be more than number of sections in the VSP model
        """

        possible_section_paramaters = ["Span", "Root_Chord", "Tip_Chord", "Dihedral", "Sweep", "Twist"]

        assert isinstance(name, str), "`name` argument must be string"
        assert name in possible_section_paramaters, f"`{name}` is not a valid parameter"
        assert isinstance(lower, float) and isinstance(upper, float), "`lower` and `upper` must be float values"
        assert lower < upper, "`upper` must be greater than `lower`"
        assert isinstance(section_id, int) and 1 <= section_id <= self.number_of_sections, f"`section_id` must be an integer between 1 and {self.number_of_sections}"

        for key in self.parameters.keys():
            if name == key:
                raise AssertionError(f"`{name}` is already added as a global parameter")
            else:
                assert f"{name}_{section_id}" not in key, f"`{name}_{section_id}` is already added as a local parameter"

        # update bounds
        lb = np.append(self.bounds[0], lower)
        ub = np.append(self.bounds[1], upper)

        # set updated variables
        self.mask = np.append(self.mask, name+f"{section_id}")
        self.bounds = (lb, ub)
        self.parameters[name+f"{section_id}"] = [f"{self.component_name}:XSec_{section_id}:{name}"]

    def add_airfoil_cst_parameters(
        self,
        surface: str,
        lower: np.ndarray,
        upper: np.ndarray,
        sections: list,
    ) -> None:
        """
            Method to add entire upper or lower airfoil as a parameter for analysis

            Currently, airfoil is parameterized using CST method only

            Parameters
            ----------
            surface: str
                a string denoting whether to add `upper` or `lower` airfoil surface

            lower: np.ndarray
                a 1D numpy array containing lower bound for all the CST coefficients

            upper: np.ndarray
                a 1D numpy array containing upper bound for all the CST coefficients

            sections: list
                a list containing section indices whose CST coefficients are to be 
                added as a parameter. `NOTE`: If there is more than one entry in
                this list, then all those sections will have uniform airfoil shape
                parametrization
        """

        # some initial checks
        assert surface.lower() in ["upper", "lower"], "`surface` must be a either 'upper' or 'lower'"
        assert isinstance(lower, np.ndarray) and isinstance(upper, np.ndarray) and lower.ndim == upper.ndim == 1 and lower.shape[0] == upper.shape[0], "`lower` and `upper` arugment must be a 1D numpy array of same shape"
        assert np.all(lower < upper), "`upper` bound must be greater than `lower` bound"
        assert isinstance(sections, list), "`sections` argument must be a list"

        surface = surface.lower()

        # select CST function/variables dynamically
        if surface == "upper":
            degree_func = self.vsp_model.GetUpperCSTDegree
            group_prefix = "UpperCoeff"
            param_prefix = "Au"
        else:
            degree_func = self.vsp_model.GetLowerCSTDegree
            group_prefix = "LowerCoeff"
            param_prefix = "Al"

        # individual section checks
        for sec in sections:

            # check list entries
            assert isinstance(sec, int) and 0 <= sec <= self.number_of_sections, f"The entires in `sections` list must be an integer between 0 and {self.number_of_sections}"

            # get xsec id
            xsec_id = self.vsp_model.GetXSec(self.xsec_surf_id, sec)

            # check number xsec type
            assert self.vsp_model.GetXSecShape(xsec_id) == self.vsp_model.XS_CST_AIRFOIL, f"airfoil at section {sec} is not CST"

            num_coeff = degree_func(xsec_id) + 1 # get number of cst coeffs

            # check number of cst coeffs
            assert len(lower) == num_coeff and len(upper) == num_coeff, (
                f"lower/upper size must match number of CST coefficients ({num_coeff}) at section {sec}"
            )

        # empty array for storing
        mask = np.array([])

        # loop through each coefficient
        for i in range(num_coeff):

            dvgeo_name = []

            for sec in sections:

                dvgeo_name.append(f"{self.component_name}:{group_prefix}_{sec}:{param_prefix}_{i}")

                # check if this parameter is already added
                for key in self.parameters.keys():
                    assert f"airfoil_xsec_{sec}_{param_prefix}{i}" not in key, f"CST coefficient {i} for {surface} airfoil at section {sec} is already added as a parameter"

            # parameter name
            param_name = f"airfoil_xsec_{'_'.join(map(str,sections))}_{param_prefix}{i}"

            # update parameter list
            self.parameters[param_name] = dvgeo_name

            # mask for the parameter
            mask = np.append(mask, param_name)

        # update bounds
        lb = np.append(self.bounds[0], lower)
        ub = np.append(self.bounds[1], upper)

        # update variables
        self.mask = np.append(self.mask, mask)
        self.bounds = (lb, ub)

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
