*******************************
Sample generation using DAFoam
*******************************

This section explains how to use ``WingFFD`` module for generating samples using DAFoam. There are typically three
main steps involved in the process: setting up options and initializing the module, adding design variables 
and generating samples. The ``WingFFD`` module is demonstrated using a CRM wing example. This script is
available in ``examples`` directory of the repository on github.

Setting up options
------------------

First step involves creating options dictionary which is used for initializating the module. There are five
mandatory options: ``solverOptions``, ``ffdFile``, ``gridFile``, ``liftIndex`` and ``aeroProblem``, rest all are optional,
please refer :ref:`options<options>` section for more details. Following snippet of the code shows an example::

    from blackbox import WingFFD
    from baseclasses import AeroProblem
    import numpy as np

    solverOptions = {
        "designSurfaces": ["wing"],
        "solverName": "DARhoSimpleCFoam",
        "primalMinResTol": 1e-6,
        "primalMinResTolDiff": 1e2,
        "checkMeshThreshold": {
            "maxAspectRatio": 1000.0,
            "maxNonOrth": 70.0,
            "maxSkewness": 4.0,
            "maxIncorrectlyOrientedFaces": 0,
    }

    # Creating aeroproblem for adflow
    # Chord ref is 1.0 since all the dimensions are scaled according to it
    ap = AeroProblem(
        name="crm", alpha=2.0, mach=0.85, reynolds=5e6, reynoldsLength=1.0, T=298.15, 
        areaRef=3.407014, chordRef=1.0, evalFuncs=["cl", "cd", "cmy"], xRef=1.2077, yRef=0.0, zRef=0.007669
    )

    # Options for blackbox
    options = {
        "solver": "dafoam",
        "solverOptions": solverOptions,
        "openfoamDir": ".",
        "gridFile": "crm_volMesh.cgns",
        "ffdFile": "crm_ffd.xyz",
        "liftIndex": 3, # Very important
        "aeroProblem": ap,
        "noOfProcessors": 8,
        "writeDeformedFFD": True,
        "getSurfaceForces": True
    }

    # Initialize the class
    wing = WingFFD(options=options)

Firstly, required packages and modules are imported. Then, ``solverOptions`` dictionary is created, refer DAFoam documentation
for more details. Then, `AeroProblem <https://mdolab-baseclasses.readthedocs-hosted.com/en/latest/pyAero_problem.html>`_
object is created which contains details about the flow conditions and the desired output variables are 
defined using ``evalFuncs`` argument. Then, ``options`` dictionary is created, refer :ref:`options<options>` 
section for more details.

  .. note::
    In solver options dictionary for DAFoam, no need to mention arguments such as ``primalBC``, ``function``, 
    ``normalizeStates``, ``inputInfo``. These options are internally set before each analysis.

    .. note::
    In the options dictionary, you can provide the location of ``0``, ``constant``, and ``system`` directories 
    (required by openfoam) by using ``openfoamDir`` argument. This is the location where blackbox will look for these directories.
    The default location is ``"."`` i.e. current working directory, refer :ref:`options<options>` section for more details.

Adding design variables
-----------------------

Next step is to add design variables based on which samples will be generated. The ``addDV`` method needs three arguments:

- ``name (str)``: name of the design variable to add. The available design variables are:

    - ``shape``: FFD control points which parameterize the airfoil shape
    - ``twist``: Twist of the airfoil sections along the wing span, number of sections will be one less than the sections defined in the FFD file since root is fixed
    - ``alpha``: Angle of attack for the analysis
    - ``mach``: Mach number for the analysis
    - ``altitude``: Altitude for the analysis

- ``lowerBound (numpy array or float)``: lower bound for the variable
- ``upperBound (numpy array or float)``: upper bound for the variable

    .. note::
        When ``shape`` variable is to be added, the lower and upper bound should be a 1D numpy array of the same size 
        as the number of FFD points. The number of FFD points can be accessed via ``nffd`` attribute of the class.

        When ``twist`` variable is to be added, the lower and upper bound should be a 1D numpy array of the same size 
        as the number of section defined in the FFD file minus one. The twist is defined in degrees. The number of twist
        sections can be accessed via ``nTwist`` attribute of the class.

        For other cases, lower and upper bound should be float.

Following code snippet adds ``alpha`` and ``mach`` as design variables::

    # Add alpha as a design variable
    wing.addDV("alpha", lowerBound=1.5, upperBound=3.5)

    # Add mach as a design variable
    wing.addDV("mach", lowerBound=0.7, upperBound=0.8)

Generating samples and accessing data
---------------------------------------

After adding design variables, generating samples is very easy. You just need to use ``generateSamples`` 
method from the initialized object. This method has two arguments:

- ``numSamples (int)``: number of samples to generate
- ``doe (numpy array)``: 2D numpy array in which each row represents a specific sample

.. note::
    You can either provide ``numSamples`` or ``doe`` i.e. both of them are mutually exclusive.
    If both are provided, then an error will be raised.

Typically, ``numSamples (int)`` should be used for generating samples. This option will internally generate doe based on the 
options provided while initializating the module. In some cases, you might want to generate samples based on your own doe. In that
case, you use ``doe (numpy array)`` argument. Following snippet of the code will generate 5 samples using internally generated doe::

    wing.generateSamples(numSamples=5)

You can see the following output upon successful completion of sample generation process:

- A folder with the name specificed in the ``directory`` option (or the default name - *output*) is created. This folder contains all the generated
  files/folders.

- Within the main output folder, there will be subfolders equal to the number of samples you requested. Each of the folder corresponds to the specific
  analysis performed. It will contain log.txt which contains the output from mesh generation and solver. There will be other files depending on the 
  options provided to solver and blackbox.

- ``data.mat`` file which contains:

  - **Input variable**: a 2D numpy array ``x`` in which each row represents a specific sample based on which analysis is performed. The number
    of rows will be usually equal to the number of samples argument in the ``generateSamples`` method. But, many times few of the analysis
    fail. It depends a lot on the solver options, so set those options after some tuning.

    .. note::
        The order of values in each row is based on how you add design variables. In this tutorial, first ``alpha`` is added as
        design variable and then shape coefficients are added. Thus, first value in each row will be alpha, next ``nffd``
        values will be FFD coefficients, and then ``nTwist`` values will be twist values.

  - **Outputs**: There are two kinds of outputs - mandatory and user specificed. The ``evalFuncs`` argument in the aero problem
    decides the user desired outputs. Along with these outputs, `volume` of the wing is the mandatory output. Following snippet 
    shows how to access the data.mat file. In this tutorial, ``evalFuncs`` argument contains ``cl``, ``cd``, ``cmy``. So, data.mat 
    will contain these variables, along with ``volume``::

        from scipy.io import loadmat
        data = loadmat("data.mat") # mention the location of mat file

        x = data["x"]
        cl = data["cl"]
        cd = data["cd"]
        cmy = data["cmy"]
        volume = data["volume"]

- ``description.txt``: contains various informations about the sample generation such as design variables, bounds, number of failed analysis, etc.

Extracting surface force data
-------------------------------

To get forces at each surface mesh coordinate, set the option ``getSurfaceForces`` to True, refer :ref:`options<options>` 
section for more details. This creates a ``surfaceForces.mat`` file in each corresponding analysis folder. This mat file contains two
variables: ``forces`` and ``surfaceCoords``. The ``forces`` is a 2d numpy of size Ncoords x 3, where Ncoords is the number of surface 
mesh points and 3 corresponds to x, y, and z component of the force vector. The ``surfaceCoords`` is also a 2d numpy array of size 
Ncoords x 3, containing surface mesh coordinates. These forces and surface coordinates can be extracted as shown in previous section 
using ``loadmat`` function from scipy.