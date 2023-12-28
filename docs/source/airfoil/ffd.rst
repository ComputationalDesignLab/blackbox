***************************
Sample generation with FFD
***************************

This section explains how to use ``AirfoilFFD`` module for generating samples. There are typically three
main steps involved in the process: setting up options and initializing the module, adding design variables 
and generating samples.

Setting up options
------------------

First step involves creating options dictionary which is used for initializating the module. The ``airfoilFile``
and ``nffd`` are the two mandatory options, rest all are optional, please refer :ref:`options<options>` 
section for more details. Following snippet of the code shows an example::

    solverOptions = {
        # Common Parameters
        "monitorvariables": ["cl", "cd", "cmz", "yplus"],
        "writeTecplotSurfaceSolution": True,
        "writeSurfaceSolution": False,
        "writeVolumeSolution": False,
        # Physics Parameters
        "equationType": "RANS",
        "smoother": "DADI",
        "MGCycle": "sg",
        "nsubiterturb": 10,
        "nCycles": 7000,
        # ANK Solver Parameters
        "useANKSolver": True,
        "ANKSubspaceSize": 400,
        "ANKASMOverlap": 3,
        "ANKPCILUFill": 4,
        "ANKJacobianLag": 5,
        "ANKOuterPreconIts": 3,
        "ANKInnerPreconIts": 3,
        # NK Solver Parameters
        "useNKSolver": True,
        "NKSwitchTol": 1e-6,
        "NKSubspaceSize": 400,
        "NKASMOverlap": 3,
        "NKPCILUFill": 4,
        "NKJacobianLag": 5,
        "NKOuterPreconIts": 3,
        "NKInnerPreconIts": 3,
        # Termination Criteria
        "L2Convergence": 1e-14
    }

    meshingOptions = {
        # ---------------------------
        #        Input Parameters
        # ---------------------------
        "unattachedEdgesAreSymmetry": False,
        "outerFaceBC": "farfield",
        "autoConnect": True,
        "BC": {1: {"jLow": "zSymm", "jHigh": "zSymm"}},
        "families": "wall",
        # ---------------------------
        #        Grid Parameters
        # ---------------------------
        "N": 129,
        "s0": 1e-6,
        "marchDist": 100.0,
    }

    # Creating aeroproblem for adflow
    ap = AeroProblem(
        name="ap", alpha=2.0, mach=0.734, reynolds=6.5e6, reynoldsLength=1.0, T=288.15, 
        areaRef=1.0, chordRef=1.0, evalFuncs=["cl", "cd", "cmz"], xRef = 0.25, yRef = 0.0, zRef = 0.0
    )

    nffd = 20 # Number of FFD points

    # Options for blackbox
    options = {
        # Requried options
        "airfoilFile": "rae2822.dat",
        "nffd": 20,
        # FFD Box options
        "fitted": True,
        "ymarginl": 0.015,
        "ymarginu": 0.015,
        # Sampling options
        "samplingCriterion": "ese",
        # Fixing the LE/TE
        "fixLETE": True,
        # Smoothing options
        "smoothing": True,
        "smoothingTolerance": 5e-4,
        "smoothingTheta": 0.6,
        # Other options
        "noOfProcessors": 8,
        "aeroProblem": ap,
        "solverOptions": solverOptions,
        "meshingOptions": meshingOptions,
        "writeAirfoilCoordinates": True,
        "plotAirfoil": True,
        "writeDeformedFFD": True,
    }

    # Example for generating samples
    airfoil = AirfoilFFD(options=options)

Firstly, required packages and modules are imported. Then, ``solverOptions`` and ``meshingOptions`` are 
created which determine the solver and meshing settings. Refer `ADflow <https://mdolab-adflow.readthedocs-hosted.com/en/latest/options.html>`_
and `pyHyp <https://mdolab-pyhyp.readthedocs-hosted.com/en/latest/options.html>`_ options for more details.
Then, `AeroProblem <https://mdolab-baseclasses.readthedocs-hosted.com/en/latest/pyAero_problem.html>`_
object is created which contains details about the flow conditions and the desired output variables are 
defined using ``evalFuncs`` argument. Then, ``options`` dictionary is created, refer :ref:`options<options>` 
section for more details. Finally, the ``AirfoilCST`` module is initialized using the options dictionary.

Adding design variables
-----------------------

Next step is to add design variables based on which samples will be generated. The ``addDV`` method needs three arguments:

- ``name (str)``: name of the design variable to add. The available design variables are:

    - ``shape``: FFD control points which parameterize the airfoil shape
    - ``alpha``: Angle of attack for the analysis
    - ``mach``: Mach number for the analysis
    - ``altitude``: Altitude for the analysis
- ``lowerBound (numpy array or float)``: lower bound for the variable
- ``upperBound (numpy array or float)``: upper bound for the variable

    .. note::
        When ``shape`` variable is to be added, the lower and upper bound should be a 1D numpy array of the same size 
        as the number of FFD points mentioned in the ``options`` dictionary. For other cases, lower
        and upper bound should be float.

Following code adds ``alpha`` and ``shape`` as design variables::

    airfoil.addDV("alpha", 2.0, 3.0)

    lb = -np.ones(nffd) * 0.1
    ub = np.ones(nffd) * 0.1
    airfoil.addDV("shape", lowerBound=lb, upperBound=ub)

Here, the upper and lower bound for ``shape`` variable is set to 0.1 and -0.1, respectively.
You can also remove a design variable using ``removeDV`` method. It takes only one input which is the name of the variable.

Generating samples and accessing data
---------------------------------------

After adding design variables, generating samples is very easy. You just need to use ``generateSamples`` 
method from the initialized object. This method has two arguments:

- ``numSamples (int)``: number of samples to generate
- ``doe (numpy array)``: 2D numpy array in which each row represents a specific sample

.. note::
    You can either provide ``numSamples`` or ``doe`` i.e. both them are mutually exclusive.
    If both are provided, then an error will be raised.

Typically, ``numSamples (int)`` should be used for generating samples. This option will internally generate doe based on the 
options provided while initializating the module and run the analysis. In some cases, you might want to generate samples based on your own doe. In that
case, you use ``doe (numpy array)`` argument. Following snippet of the code will generate 10 samples::

    airfoil.generateSamples(numSamples=10)

You can see the following output upon successful completion of sample generation process:

- A folder with the name specificed in the ``directory`` option (or the default name - *output*) is created. This folder contains all the generated
  files/folders.

- Within the main output folder, there will be subfolders equal to the number of samples you requested. Each of the folder corresponds to the specific
  analysis performed. It will contain log.txt which contains the output from mesh generation and solver. There will be other files depending on the 
  options provided to solver and blackbox.

- ``data.mat`` file which contains:

    - **Input variable**: a 2D numpy array ``x`` in which each row represents a specific sample based on which analysis is performed. The number
      of rows will be usually equal to the number of samples argument in the ``generateSamples`` method. But, many times few of the analysis
      fail. It depends a lot on the solver and meshing options, so set those options after some tuning.

      .. note::
          The order of values in each row is based on how you add design variables. In this tutorial, first ``alpha`` is added as
          design variable. Then, lower and upper surface CST coefficients are added. Thus, first value in each row will be alpha, next 6
          values will be upper surface CST coefficients and last 6 will be lower surface CST coefficients.

    - **Output variables**: There are two kinds of output variables - mandatory and user specificed. The ``evalFuncs`` argument in the aero problem
      decides the user desired variables. Along with these variables, `area` of the airfoil is the mandatory objective.


  Following snippet shows how to access the data.mat file. In this tutorial, ``evalFuncs`` argument contains 
  ``cl``, ``cd``, ``cmz``. So, data.mat will contain these variables, along with ``area``::

    from scipy.io import loadmat
    data = loadmat("data.mat") # mention the location of mat file

    x = data["x"]
    cl = data["cl"]
    cd = data["cd"]
    cmz = data["cmz"]
    area = data["area"]

- ``ffd.xyz``: contains the coordinates of the FFD box in plot3D format.

- ``description.txt``: contains various informations about the sample generation such as design variables, bounds, number of failed analysis, etc.
