import numpy as np
from scipy.special import factorial
from typing import Union, List

class CST:

    def __init__(self, airfoil_file: str, num_cst: Union[int, List[int]]):
        """
            Class implementing a basic 2D Class Shape Transformation method for airfoil parameterization

            Following are the assumptions regarding airfoil dat file:

                - The airfoil dat file should be in selig format i.e. coordinates should 
                  start from (1.0,0.0), move in counter-clockwise direction and end at (1.0,0.0)
                - The airfoil can have sharp or blunt trailing edge

            This is a simplified version of the `DVGeometryCST` module in the pyGeo library:
            https://mdolab-pygeo.readthedocs-hosted.com/en/latest/DVGeometryCST.html

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
        assert isinstance(airfoil_file, str), "`airfoil_file` should be a str"

        if isinstance(num_cst, int):
            self.num_cst_upper = num_cst
            self.num_cst_lower = num_cst
        elif isinstance(num_cst, list):
            self.num_cst_upper = num_cst[0]
            self.num_cst_lower = num_cst[1]
        else:
            raise ValueError("`num_cst` should be an int or a list of int with two entries")
        
        coords = self.read_coord_file(airfoil_file)

        self.orig_coords = coords.copy()

        # Some validation for coordinate file
        assert self.orig_coords[0,0] == self.orig_coords[-1,0], "The X coordinate of airfoil doesn't start and end at same point."
        assert self.orig_coords[0,0] == 1.0 and self.orig_coords[0,1] == 0.0, "The coordinates of airfoil doesn't start at (1.0, 0.0)"
        assert np.min(self.orig_coords[:,0]) >= 0.0 and np.max(self.orig_coords[:,0]) <= 1.0, "The X coordinates of airfoil are not in range [0,1]"
        assert self.orig_coords[np.argmin(self.orig_coords[:,0]),1] == 0.0, "The Y coordinate of airfoil at the LE is not 0.0"

        # LE index
        idx_le = np.argmin(coords[:,0])

        # Split coordinates using the LE index
        self.upper_coords = coords[:idx_le,:]
        self.lower_coords = coords[idx_le:,:]

        ######### Upper surface pre-processing

        # find index of TE cords
        idx_te = np.where(self.upper_coords[:,0] == self.upper_coords[0,0])[0]

        # Number of points in TE
        num_te_pts = len(idx_te)

        if num_te_pts > 1:
            self.upper_te_coords = self.upper_coords[idx_te,:]
            self.upper_yte = self.upper_te_coords[-1,1] - self.upper_te_coords[0,1]
            self.upper_coords = self.upper_coords[num_te_pts-1:,:]
        else:
            self.upper_te_coords = None
            self.upper_yte = 0.0

        # find CST coefficents for upper surface
        self.upper_cst = self.compute_cst_coefficients(self.upper_coords[:,0], self.upper_coords[:,1], self.upper_yte, self.num_cst_upper)

        ######### Lower surface pre-processing

        # find index of TE cords
        idx_te =  np.where(self.lower_coords[:,0] == self.lower_coords[-1,0])[0]

        # Number of points in TE
        num_te_pts = len(idx_te)

        if num_te_pts > 1:
            self.lower_te_coords = self.lower_coords[idx_te,:]
            self.lower_yte = self.lower_te_coords[0,1] - self.lower_te_coords[-1,1] # np.linalg.norm(self.lower_te_coords[0,1] - self.lower_te_coords[-1,1])
            self.lower_coords = self.lower_coords[:-num_te_pts+1,:]
        else:
            self.lower_te_coords = None
            self.lower_yte = 0.0

        # find CST coefficents for lower surface
        self.lower_cst = self.compute_cst_coefficients(self.lower_coords[:,0], self.lower_coords[:,1], self.lower_yte, self.num_cst_lower)

        print(f"### CST coefficients for the coordinates in {airfoil_file} ########")
        print(f"Upper surface: {self.upper_cst}")
        print(f"Lower surface: {self.lower_cst}")

        # Other CST parameters
        self.n1 = 0.5
        self.n2 = 1.0

    def compute_airfoil_coordinates(self, upper_cst_coeffs: np.ndarray, lower_cst_coeffs: np.ndarray) -> np.ndarray:
        """
            Compute the airfoil coordinates based on the given CST coefficients and
            return the coordinates in the selig format

            Parameters
            ----------
            upper_cst_coeffs

            dv_dict: dict
                dictionary containing the name and value of the variable(s)

            Returns
            -------
            coords: np.ndarray
                a 2D numpy array containing x and y coordinates of the airfoil
        """

        assert isinstance(upper_cst_coeffs, np.ndarray) and isinstance(lower_cst_coeffs, np.ndarray), "given upper and lower CST must be numpy array"
        assert upper_cst_coeffs.ndim == 1 and lower_cst_coeffs.ndim == 1, "given upper and lower CST arrays must be 1D"
        assert upper_cst_coeffs.shape[0] == self.num_cst_upper, "given number of upper CST coefficients is not correct"
        assert lower_cst_coeffs.shape[0] == self.num_cst_lower, "given number of lower CST coefficients is not correct"

        # Get upper surface coordinates
        upper_ycoords = self.compute_cst_coordinates(self.upper_coords[:,0], upper_cst_coeffs, self.n1, self.n2, self.upper_yte)
        
        upper = np.hstack(( self.upper_coords[:,0].reshape(-1,1), upper_ycoords.reshape(-1,1) ))

        # Append upper TE points
        if self.upper_te_coords is not None:
            upper = np.vstack(( self.upper_te_coords[:-1,:] , upper ))
        
        # Get lower surface coordinates
        lower_ycoords = self.compute_cst_coordinates(self.lower_coords[:,0], lower_cst_coeffs, self.n1, self.n2, self.lower_yte)
    
        lower = np.hstack(( self.lower_coords[:,0].reshape(-1,1), lower_ycoords.reshape(-1,1) ))

        # Append lower TE points
        if self.lower_te_coords is not None:
            lower = np.vstack(( lower, self.lower_te_coords[1:,:] ))

        coords = np.vstack(( upper, lower ))

        return coords
    
    def compute_cst_coefficients(self, x: np.ndarray, y: np.ndarray, yte: float, num_coeffs: int) -> np.ndarray:
        """
            Compute the CST coefficients for given the airfoil upper/lower surface coordinates

            The input should be either upper or lower airfoil surface, not both

            Parameters
            ----------
            x: ndarray (# pts,)
                x coordinates of the curve
            y: ndarray (# pts,)
                y coordinates of the curve
            yte: float
                trailing edge thickness of the airfoil
            num_coeffs: int
                number of CST coefficients to use for fitting
        
            Returns
            -------
            ndarray (# coeff)
                CST coefficients fitted to the given curve
        """

        C = self.compute_class_functions(x)

        S = self.compute_shape_functions(x, np.ones(num_coeffs))

        A = C*S

        w = np.linalg.lstsq(A.transpose(), y - x*yte, rcond=None)[0]

        return w

    def compute_cst_coordinates(self, x: np.ndarray, w: np.ndarray, n1: float, n2: float, yte: float) -> np.ndarray:
        """
            Compute the y coordinates for given x coordinates and CST coefficients

            The x coordinates should in range [0,1] and yte is normalized by the chord

            Parameters
            ----------
            x : ndarray (# pts,)
                x coordinates at which to compute the CST curve height
            w : ndarray (# coeff,)
                CST coefficient array
            n1 : float
                First class shape parameter, default=0.5
            n2 : float
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
        
        C = self.compute_class_functions(x, n1, n2)

        S = self.compute_shape_functions(x, w)

        return C * S.sum(axis=0) + yte * x

    def compute_class_functions(self, x: np.ndarray, n1: float = 0.5, n2: float = 1.0) -> np.ndarray:
        """
            Compute the class shape of a CST curve

            Parameters
            ----------
            x : ndarray (# pts,)
                x coordinates at which to compute the CST curve height
            n1 : float
                First class shape parameter, default = 0.5
            n2 : float
                Second class shape parameter, default = 1.0

            Returns
            -------
            ndarray (# pts,)
                y coordinates of the class shape
        """

        return x ** n1 * (1.0 - x) ** n2

    def compute_shape_functions(self, x: np.ndarray, w: np.ndarray):
        """
            Compute the Bernstein polynomial shape function of a CST curve

            This function assumes x has been normalized to the range [0,1].

            Parameters
            ----------
            x : ndarray (# pts,)
                x coordinates at which to compute the CST curve height
            w : ndarray (# coeff,)
                CST coefficient array

            Returns
            -------
            ndarray (# coeff, # pts)
                Bernstein polynomials for each CST coefficient
        """

        num_coeffs = len(w) # number of coefficients

        order = num_coeffs - 1 # order of the polynomial

        S = np.zeros((num_coeffs, len(x)))

        facts = factorial(np.arange(0, order + 1))

        for i in range(num_coeffs):
            binom = facts[-1] / (facts[i] * facts[order - i])
            S[i] = w[i] * binom * x ** (i) * (1.0 - x) ** (order - i)

        return S

    @staticmethod
    def read_coord_file(filename: str, headerlines: int = 0):
        """
            Method to read an airfoil dat file in selig format.

            Each (x,y) coordinate should be on a new line

            Parameters
            ----------
            filename: str
                the file to read from, including the '.dat' extension

            headerlines: int
                the number of lines to skip at the beginning of the file

            Returns
            -------
            x: ndarray [N,2]
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

        return np.array(r)
