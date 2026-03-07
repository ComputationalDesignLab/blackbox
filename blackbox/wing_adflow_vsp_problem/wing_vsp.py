import os
from pygeo import DVGeometryVSP
import numpy as np

class WingVSP():

    def __init__(self, vsp_file:str):

        assert isinstance(vsp_file, str), "`vsp_file` argument must be string"

        vsp_file = os.path.abspath(vsp_file) # set abs path

        # check if vsp file path exists
        if not os.path.exists(vsp_file):
            raise ValueError("Provided vsp file path does not exist")

        # Initializing the dvgeo_vsp object
        self.dvgeo_vsp = DVGeometryVSP(vsp_file)

        # Some checks
        geoms = self.dvgeo_vsp.vspModel.FindGeoms()
        assert(len(geoms) == 1) and self.dvgeo_vsp.vspModel.GetGeomTypeName(geoms[0]) == "Wing", "Your OpenVSP model should contain only one component of type `WING`"

        # Storing some VSP variables
        self.component_id = geoms[0]
        self.component_name = self.dvgeo_vsp.vspModel.GetGeomName(geoms[0])
        self.xsec_surf_id = self.dvgeo_vsp.vspModel.GetXSecSurf(self.component_id, 0)
        self.number_of_xsec = self.dvgeo_vsp.vspModel.GetNumXSec(self.xsec_surf_id)
        self.number_of_sections = self.number_of_xsec - 1

        self.parameters = []
        self.mask = np.array([])
        self.bounds = (np.array([]), np.array([]))

        self.possible_section_paramaters = ["Span", "Dihedral", "Sweep", "Twist"]

    def add_parameter(
        self,
        name: str,
        lower: float,
        upper: float,
        section_id: int = -1
    ):
        """
            Method for adding wing section parameters for analysis
            
            Some of the examples are span, dihedral, sweep, twist, chord length

            Parameters
            ----------
            name: str
                name of the parameter to be added. It can only be:
                    * `span`: span of the section
                    * `dihedral`: dihedral of the section
                    * `sweep`: sweep of the section
                    * `twist`: twist of the section

                `NOTE`: You can add a single parameter for all sections, refer to
                        `section_id` argument for more detail

            lower: float or np.ndarray
                lower bound for the parameter

            upper: float or np.ndarray
                upper bound for the parameter

            section_id: int
                an integer denoting to which section the parameters belongs to.
                The default is -1, which indicates that a single value controls
                this parameter across all the sections. `NOTE`: This parameter must
                not be more than number of sections in the VSP model
        """

        assert isinstance(name, str)
        assert isinstance(lower, float) and isinstance(upper, float)
        assert lower < upper
        assert isinstance(section_id, int) and section_id > -2

        assert name in self.possible_section_paramaters
        assert section_id <= self.number_of_sections and section_id >= -1

        dvgeo_name = []

        if section_id == -1:

            for i in range(self.number_of_sections):
                
                self.dvgeo_vsp.addVariable(
                    component=self.component_name,
                    group=f"XSec_{i+1}",
                    parm=name,
                    scaledStep=False
                )

                dvgeo_name.append(f"{self.component_name}:XSec_{i+1}:{name}")

        else:
            
            self.dvgeo_vsp.addVariable(
                component=self.component_name,
                group=f"XSec_{section_id}",
                parm=name,
                scaledStep=False
            )

            self.parameters.append(f"{name}_{section_id}")

            dvgeo_name.append(f"{self.component_name}:XSec_{i+1}:{name}")

        self.parameters.append({
                "name": name,
                "dvgeo_name": dvgeo_name
            })

        lb = np.append(self.bounds[0], lower)
        ub = np.append(self.bounds[1], upper)
        self.bounds = (lb, ub)

    def set_parameter(self, params: dict):

        pass

        # dvgeo_params = {}

        # for key, value in params.items():

        #     if key in 
        
            

        #     else:

        #         dvgeo_params[key] = value

