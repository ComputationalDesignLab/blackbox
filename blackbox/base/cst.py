import numpy as np
from scipy.special import factorial

class CST():

    def __init__(self, airfoilFile, numCST):
        """
            Class implementing a basic 2D Class Shape Transformation method for airfoil parameterization

            Following are the assumptions regarding airfoil dat file:

                - The airfoil dat file should be in selig format i.e. coordinates should 
                  start from (1.0,0.0), move in counter-clockwise direction and end at (1.0,0.0)
                - The airfoil can have sharp or blunt trailing edge

            This is a simplified version of the `DVGeometryCST` module in the pyGeo library

            Parameters
            ----------
            airfoilFile: str
                name of the dat file that contains airfoil coordinates
            numCST: int or list of int
                number of CST coefficients for parameterizing the airfoil. If ``numCST`` is an int,
                the value will be used for both upper and lower. If it is a list of two integers, the 
                first value defines the number of CST coefficients for upper surface and the second is 
                the number of CST coefficients for lower surface
        """

        # Checks
        if not isinstance(airfoilFile, str):
            raise ValueError("`airfoilFile` should be a str")

        if isinstance(numCST, int):
            self.nCSTUpper = numCST
            self.nCSTLower = numCST
        elif isinstance(numCST, list):
            self.nCSTUpper = numCST[0]
            self.nCSTLower = numCST[1]
        else:
            raise ValueError("`numCST` should be an int or a list of int with two entries")
        
        coords = readCoordFile(airfoilFile)

        self.origCoords = coords.copy()

        # LE index
        idxLE = np.argmin(coords[:,0])

        # Split coordinates using the LE index
        self.upperCoords = coords[:idxLE,:]
        self.lowerCoords = coords[idxLE:,:]

        ######### Upper surface

        # find index of TE cords
        idxTE =  np.where(self.upperCoords[:,0] == self.upperCoords[0,0])[0]

        # Number of points in TE
        numTE = len(idxTE)

        if numTE > 1:
            self.upperTECoords = self.upperCoords[idxTE,:]
            self.upperYTE = self.upperTECoords[-1,1] - self.upperTECoords[0,1] # np.linalg.norm(self.upperTECoords[0,1] - self.upperTECoords[-1,1])
            self.upperCoords = self.upperCoords[numTE-1:,:]
        else:
            self.upperTECoords = None
            self.upperYTE = 0.0

        # find CST coefficents for upper surface
        self.upperCST = self.computeCSTCoefficients(self.upperCoords[:,0], self.upperCoords[:,1], self.upperYTE, self.nCSTUpper)

        ######### Lower surface pre-processing

        # find index of TE cords
        idxTE =  np.where(self.lowerCoords[:,0] == self.lowerCoords[-1,0])[0]

        # Number of points in TE
        numTE = len(idxTE)

        if numTE > 1:
            self.lowerTECoords = self.lowerCoords[idxTE,:]
            self.lowerYTE = self.lowerTECoords[0,1] - self.lowerTECoords[-1,1] # np.linalg.norm(self.lowerTECoords[0,1] - self.lowerTECoords[-1,1])
            self.lowerCoords = self.lowerCoords[:-numTE+1,:]
        else:
            self.lowerTECoords = None
            self.lowerYTE = 0.0

        # find CST coefficents for lower surface
        self.lowerCST = self.computeCSTCoefficients(self.lowerCoords[:,0], self.lowerCoords[:,1], self.lowerYTE, self.nCSTLower)

        print(f"### CST coefficients for the coordinates in {airfoilFile} ########")
        print(f"Upper surface: {self.upperCST}")
        print(f"Lower surface: {self.lowerCST}")

        # Default values of other DVs
        self.N1 = np.array([0.5])
        self.N2 = np.array([1.0])

        # Initializing the variables
        self.variables = {
            "upper": self.upperCST,
            "lower": self.lowerCST,
            "N1": self.N1,
            "N2": self.N2
        }

        self.DVs = []

    def computeCSTCoefficients(self, x, y, yte, numCoeffs):
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

        C = self.computeClassFunctions(x)

        S = self.computeShapeFunctions(x, np.ones(numCoeffs))

        A = C*S

        w = np.linalg.lstsq(A.transpose(), y - x*yte, rcond=None)[0]

        return w


    def computeCSTCoordinates(self, x, w, N1, N2, yte):
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
        
        C = self.computeClassFunctions(x, N1, N2)

        S = self.computeShapeFunctions(x, w)

        return C * S.sum(axis=0) + yte * x


    def computeClassFunctions(self, x, N1=0.5, N2=1.0):
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
    
    
    def addDV(self, name):#, lowerBound, upperBound):
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


    def removeDV(self, name):#, lowerBound, upperBound):
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


    def setDesignVars(self, DVs):
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


    def update(self, dvDict):
        """
            Compute the airfoil coordinates based on the given design variables and
            return the coordinates in the selig format

            Parameters
            ----------
            dvDict: dict
                a dictionary containing the name and value of the variable(s)
        """

        # Set design variable
        # self.setDesignVars(dvDict)

        # Get upper surface coordinates
        upperYCoords = self.computeCSTCoordinates(self.upperCoords[:,0], self.variables["upper"], 
                                   self.variables["N1"], self.variables["N2"], self.upperYTE)
        
        upper = np.hstack(( self.upperCoords[:,0].reshape(-1,1), upperYCoords.reshape(-1,1) ))

        # Append upper TE points
        if self.upperTECoords is not None:
            upper = np.vstack(( self.upperTECoords[:-1,:] , upper ))
        
        # Get lower surface coordinates
        lowerYCoords = self.computeCSTCoordinates(self.lowerCoords[:,0], self.variables["lower"], 
                                   self.variables["N1"], self.variables["N2"], self.lowerYTE)
    
        lower = np.hstack(( self.lowerCoords[:,0].reshape(-1,1), lowerYCoords.reshape(-1,1) ))

        # Append lower TE points
        if self.lowerTECoords is not None:
            lower = np.vstack(( lower, self.lowerTECoords[1:,:] ))

        coords = np.vstack(( upper, lower ))

        # Reset the DV values
        self.variables = {
            "upper": self.upperCST,
            "lower": self.lowerCST,
            "N1": self.N1,
            "N2": self.N2
        }

        return coords


    def getNDV(self):
        """
            Method to compute the total number of variables added for parameterization
        """

        return len(self.DVs)


    def computeShapeFunctions(self, x, w, dtype=float):
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


def readCoordFile(filename, headerlines=0):
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
