import numpy as np
from scipy.special import factorial
from ..msg import print_msg

class CST():

    def __init__(self, airfoil_file, num_cst):
        """
            Class implementing a basic 2D Class Shape Transformation method for airfoil parameterization

            Following are the assumptions regarding airfoil dat file:

                - The airfoil dat file should be in selig format i.e. coordinates should 
                  start from (1.0,0.0), move in counter-clockwise direction and end at (1.0,0.0)
                - The airfoil can have sharp or blunt trailing edge

            This is a simplified version of the `DVGeometryCST` module in the pyGeo library

            Parameters
            ----------
            airfoil_file: str
                name of the dat file that contains airfoil coordinates
            num_cst: int or list of int
                number of CST coefficients for parameterizing the airfoil. If ``num_cst`` is an int,
                the value will be used for both upper and lower. If it is a list of two integers, the 
                first value defines the number of CST coefficients for upper surface and the second is 
                the number of CST coefficients for lower surface
        """

        # Checks
        if not isinstance(airfoil_file, str):
            raise ValueError("`airfoil_file` should be a str")

        if isinstance(num_cst, int):
            self.nCST_upper = num_cst
            self.nCST_lower = num_cst
        elif isinstance(num_cst, list):
            self.nCST_upper = num_cst[0]
            self.nCST_lower = num_cst[1]
        else:
            raise ValueError("`num_cst` should be an int or a list of int with two entries")
        
        coords = self.read_coord_file(airfoil_file)

        self.orig_coords = coords.copy()

        # Some validation for coordinate file
        if self.orig_coords[0,0] != self.orig_coords[-1,0]:
            print_msg("The X coordinate of airfoil doesn't start and end at same point.")
        
        if self.orig_coords[0,1] != self.orig_coords[-1,1]:
            print_msg("The Y coordinate of airfoil doesn't start and end at same point.")

        if self.orig_coords[0,0] != 1.0 or self.orig_coords[0,1] != 0.0:
            print_msg("The coordinates of airfoil doesn't start at (1.0, 0.0)")

        if np.min(self.orig_coords[:,0]) < 0.0 or np.max(self.orig_coords[:,0]) > 1.0:
            print_msg("The X coordinates of airfoil are not in range [0,1]")

        if self.orig_coords[np.argmin(self.orig_coords[:,0]),1] != 0.0:
            print_msg("The Y coordinate of airfoil at the LE is not 0.0")

        # LE index
        idxLE = np.argmin(coords[:,0])

        # Split coordinates using the LE index
        self.upper_coords = coords[:idxLE,:]
        self.lower_coords = coords[idxLE:,:]

        ######### Upper surface

        # find index of TE cords
        idxTE =  np.where(self.upper_coords[:,0] == self.upper_coords[0,0])[0]

        # Number of points in TE
        numTE = len(idxTE)

        if numTE > 1:
            self.upperTE_coords = self.upper_coords[idxTE,:]
            self.upper_yTE = self.upperTE_coords[-1,1] - self.upperTE_coords[0,1] # np.linalg.norm(self.upperTE_coords[0,1] - self.upperTE_coords[-1,1])
            self.upper_coords = self.upper_coords[numTE-1:,:]
        else:
            self.upperTE_coords = None
            self.upper_yTE = 0.0

        # find CST coefficents for upper surface
        self.upper_CST = self.compute_cst_coefficients(self.upper_coords[:,0], self.upper_coords[:,1], self.upper_yTE, self.nCST_upper)

        ######### Lower surface pre-processing

        # find index of TE cords
        idxTE =  np.where(self.lower_coords[:,0] == self.lower_coords[-1,0])[0]

        # Number of points in TE
        numTE = len(idxTE)

        if numTE > 1:
            self.lowerTE_coords = self.lower_coords[idxTE,:]
            self.lower_yTE = self.lowerTE_coords[0,1] - self.lowerTE_coords[-1,1] # np.linalg.norm(self.lowerTE_coords[0,1] - self.lowerTE_coords[-1,1])
            self.lower_coords = self.lower_coords[:-numTE+1,:]
        else:
            self.lowerTE_coords = None
            self.lower_yTE = 0.0

        # find CST coefficents for lower surface
        self.lower_CST = self.compute_cst_coefficients(self.lower_coords[:,0], self.lower_coords[:,1], self.lower_yTE, self.nCST_lower)

        print(f"### CST coefficients for the coordinates in {airfoil_file} ########")
        print(f"Upper surface: {self.upper_CST}")
        print(f"Lower surface: {self.lower_CST}")

        # Default values of other DVs
        self.N1 = np.array([0.5])
        self.N2 = np.array([1.0])

        # Initializing the variables
        self.variables = {
            "upper": self.upper_CST,
            "lower": self.lower_CST,
            "N1": self.N1,
            "N2": self.N2
        }

        self.DVs = []

    def compute_cst_coefficients(self, x, y, yte, numCoeffs):
        """
            Compute the CST coefficients using the provided coordinates

            The input should be either upper or lower airfoil surface, not both

            Parameters
            ----------
            x : ndarray (# pts,)
                x coordinates of the curve
            y : ndarray (# pts,)
                y coordinates of the curve
        
            Returns
            -------
            ndarray (# coeff)
                CST coefficients fitted to the given curve
        """

        C = self.compute_class_functions(x)

        S = self.compute_shape_functions(x, np.ones(numCoeffs))

        A = C*S

        w = np.linalg.lstsq(A.transpose(), y - x*yte, rcond=None)[0]

        return w


    def compute_cst_coordinates(self, x, w, N1, N2, yte):
        """
            Compute the y coordinates for given x coordinates and CST coefficients

            The x coordinates should in range [0,1] and yte is normalized by the chord

            Parameters
            ----------
            x : ndarray (# pts,)
                x coordinates at which to compute the CST curve height
            w : ndarray (# coeff,)
                CST coefficient array
            N1 : float
                First class shape parameter, default=0.5
            N2 : float
                Second class shape parameter, default=1.0
            yte : float, default=0.0
                y coordinate of the trailing edge (used to define trailing edge thickness).
                Note that the trailing edge will be twice this thick, assuming the same ``yte``
                value is used for both the upper and lower surfaces.

            Returns
            -------
            ndarray (# pts,)
                y coordinates of the CST curve
        """
        
        C = self.compute_class_functions(x, N1, N2)

        S = self.compute_shape_functions(x, w)

        return C * S.sum(axis=0) + yte * x


    def compute_class_functions(self, x, N1=0.5, N2=1.0):
        """
            Compute the class shape of a CST curve

            Parameters
            ----------
            x : ndarray (# pts,)
                x coordinates at which to compute the CST curve height
            N1 : float
                First class shape parameter, default = 0.5
            N2 : float
                Second class shape parameter, default = 1.0

            Returns
            -------
            ndarray (# pts,)
                y coordinates of the class shape
        """

        return x ** N1 * (1.0 - x) ** N2
    
    
    def add_design_variable(self, name):#, lowerBound, upperBound):
        """
            Add design variable for modifying the airfoil shape

            Parameters
            ----------
            name: str
                name of the design variable. Possible names are:
                    - `upper`: upper surface CST coefficients
                    - `lower`: lower surface CST coefficients
                    - `N1`: first class shape parameter for both upper and lower surfaces
                    - `N2`: second class shape parameter for both upper and lower surfaces
        """

        if name not in ["upper", "lower", "N1", "N2"]:
            raise ValueError('name must be one of the following: \"upper\", \"lower\", \"N1\", \"N2\"')
        
        if name in self.DVs:
            raise ValueError(f'{name} already exists as a DV')
        
        # if not isinstance(lowerBound, np.ndarray):
        #     raise ValueError(f'lowerBound should be a numpy array')

        # if not isinstance(upperBound, np.ndarray):
        #     raise ValueError(f'upperBound should be a numpy array')

        # lowerBound: ndarray
        #     lower bound for the variable(s)
        # upperBound: ndarray
        #     lower bound for the variable(s)

        self.DVs.append(name)


    def remove_design_variable(self, name):#, lowerBound, upperBound):
        """
            Remove design variable which are added for modifying the airfoil shape

            Parameters
            ----------
            name: str
                name of the design variable. Possible names are:
                    - `upper`: upper surface CST coefficients
                    - `lower`: lower surface CST coefficients
                    - `N1`: first class shape parameter for both upper and lower surfaces
                    - `N2`: second class shape parameter for both upper and lower surfaces
        """
        
        if name not in self.DVs:
            raise ValueError(f'{name} is not a DV')
        
        # if not isinstance(lowerBound, np.ndarray):
        #     raise ValueError(f'lowerBound should be a numpy array')

        # if not isinstance(upperBound, np.ndarray):
        #     raise ValueError(f'upperBound should be a numpy array')

        # lowerBound: ndarray
        #     lower bound for the variable(s)
        # upperBound: ndarray
        #     lower bound for the variable(s)

        self.DVs.remove(name)


    def set_design_vars(self, DVs):
        """
            Set design variables based on the given DVs

            Parameters
            ----------
            DVs: dict
                A dictionary containing the value of design variables. Each
                value should be a 1D numpy array
        """

        # Performing checks
        if len(self.DVs) == 0:
            raise ValueError("No design variables are added")

        if not isinstance(DVs, dict):
            raise ValueError("The input should be a dictionary")
    
        for name in DVs.keys():

            if name in self.DVs: # only update the CST variables, ignore other variables

                if not isinstance(DVs[name], np.ndarray):
                    raise ValueError(f"Value of {name} should be a 1D numpy array")
                
                if DVs[name].shape != self.variables[name].shape:
                    raise ValueError(f"The input shape {DVs[name].shape} is different from \
                                    the expected shape {self.variables[name].shape} for the {name} DV")

                self.variables[name] = DVs[name]


    def update(self, dv_dict):
        """
            Compute the airfoil coordinates based on the given design variables and
            return the coordinates in the selig format

            Parameters
            ----------
            dv_dict: dict
                a dictionary containing the name and value of the variable(s)
        """

        # Set design variable
        # self.set_design_vars(dv_dict)

        # Get upper surface coordinates
        upperYCoords = self.compute_cst_coordinates(self.upper_coords[:,0], self.variables["upper"], 
                                   self.variables["N1"], self.variables["N2"], self.upper_yTE)
        
        upper = np.hstack(( self.upper_coords[:,0].reshape(-1,1), upperYCoords.reshape(-1,1) ))

        # Append upper TE points
        if self.upperTE_coords is not None:
            upper = np.vstack(( self.upperTE_coords[:-1,:] , upper ))
        
        # Get lower surface coordinates
        lowerYCoords = self.compute_cst_coordinates(self.lower_coords[:,0], self.variables["lower"], 
                                   self.variables["N1"], self.variables["N2"], self.lower_yTE)
    
        lower = np.hstack(( self.lower_coords[:,0].reshape(-1,1), lowerYCoords.reshape(-1,1) ))

        # Append lower TE points
        if self.lowerTE_coords is not None:
            lower = np.vstack(( lower, self.lowerTE_coords[1:,:] ))

        coords = np.vstack(( upper, lower ))

        # Reset the DV values
        self.variables = {
            "upper": self.upper_CST,
            "lower": self.lower_CST,
            "N1": self.N1,
            "N2": self.N2
        }

        return coords


    def getNDV(self):
        """
            Method to compute the total number of variables added for parameterization
        """

        return len(self.DVs)


    def compute_shape_functions(self, x, w, dtype=float):
        """
            Compute the Bernstein polynomial shape function of a CST curve

            This function assumes x has been normalized to the range [0,1].

            Parameters
            ----------
            x : ndarray (# pts,)
                x coordinates at which to compute the CST curve height
            w : ndarray (# coeff,)
                CST coefficient array
            dtype : type, optional
                Type for instantiated arrays, by default float

            Returns
            -------
            ndarray (# coeff, # pts)
                Bernstein polynomials for each CST coefficient
        """

        numCoeffs = len(w)

        order = numCoeffs - 1

        S = np.zeros((numCoeffs, len(x)), dtype=dtype)

        facts = factorial(np.arange(0, order + 1))

        for i in range(numCoeffs):
            binom = facts[-1] / (facts[i] * facts[order - i])
            S[i] = w[i] * binom * x ** (i) * (1.0 - x) ** (order - i)

        return S

    @staticmethod
    def read_coord_file(filename, headerlines=0):
        """
            Method to read an airfoil dat file in selig format.

            Each xy coordinate should be on a new line

            Parameters
            ----------
            filename : str
                the file to read from, including the '.dat' extension

            headerlines : int
                the number of lines to skip at the beginning of the file

            Returns
            -------
            x : Ndarray [N,2]
                The coordinates read from the file
        """

        with open(filename, "r") as f:

            for _i in range(headerlines):
                f.readline()
                
            r = []

            while True:
                line = f.readline()

                if not line:
                    break  # end of file

                if line.isspace():
                    break  # blank line

                r.append([float(s) for s in line.split()])

                X = np.array(r)

        return X